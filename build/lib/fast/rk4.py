# -*- coding: utf-8 -*-

#************************************************************************
#       Copyright (C) 2014 - 2017 Oscar Gerardo Lazo Arjona             *
#              <oscar.lazo@correo.nucleares.unam.mx>                    *
#                                                                       *
#  This file is part of FAST.                                           *
#                                                                       *
#  FAST is free software: you can redistribute it and/or modify         *
#  it under the terms of the GNU General Public License as published by *
#  the Free Software Foundation, either version 3 of the License, or    *
#  (at your option) any later version.                                  *
#                                                                       *
#  FAST is distributed in the hope that it will be useful,              *
#  but WITHOUT ANY WARRANTY; without even the implied warranty of       *
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
#  GNU General Public License for more details.                         *
#                                                                       *
#  You should have received a copy of the GNU General Public License    *
#  along with FAST.  If not, see <http://www.gnu.org/licenses/>.        *
#                                                                       *
#************************************************************************

sage_included = 'sage' in globals().keys()
if not sage_included:
	from math import sqrt
	from misc import Mu,IJ,find_phase_transformation, format_double
	from misc import calculate_iI_correspondence
	from misc import Theta,dot_product,find_omega_min,laser_detunings
	from misc import detuning_combinations
	import os

from time import time

def add_line(Ne,mu,musign,laser,l,epsilonsign,r,ir,jr,irho,jrho,conj=False):
	dp=dot_product(laser[l-1],epsilonsign,r,ir,jr)
	if dp==0: return ''
	line ='    y('+str(mu)+')= y('+str(mu)+')'+musign
	line+='E0('+str(l)+')*'
	if sage_included:
		line+='('+format_double(real(dp))+','+format_double(imag(dp))+')*'
	else:
		line+='('+format_double(dp.real)+','+format_double(dp.imag)+')*'
	
	nu=Mu(irho,jrho,s=1,N=Ne)
	if nu==0:
		#We have a term proportional to rho11
		line+='rho11'
	elif conj:
		line+='conjg(x('+str(nu)+'))'
	else:
		line+='x('+str(nu)+')'
	
	return line+'\n'

def write_rk4(path,name,laser,omega,gamma,r,Lij,states=None,verbose=1):
	r"""
    This function writes the Fortran code needed to calculate the time evolution of the density matrix elements
    `\rho_{ij}` using the Runge-Kutta method of order 4.

    INPUT:

    - ``path`` - A string with the working directory where all files will be stored. It must end with ``/``.

    - ``name`` - A string with the name of the experiment. All files produced will begin with this name.
    
    - ``laser`` - A list of Laser objects (see the Laser class).
    
    - ``omega`` - A matrix or list of lists containing the frequency differences `\omega_{ij}`.
    
    - ``gamma`` - A matrix or list of lists containing the spontaneous decay frequencies `\gamma_{ij}`.
    
    - ``r`` - A list of three matrices or lists of lists containing the components of the position operator `r_{-1ij},r_{0ij},r_{1ij}`.
    
    - ``Lij`` - A list with elements of the form ``[i,j,[l1,l2,...]]`` representing the sets `L_{ij}` of which lasers excite wich transitions. It does not need to contain an element for all ``i,j`` pairs, but only those which have a laser that excites them.
    
    - ``Omega`` - A floating point number indicating the frequency scale for the equations. The frequencies ``omega`` and ``gamma`` are divided by this number. If ``None`` the equations and the input are taken in SI units.

    OUTPUT:
    
    - A file ``name.f90`` is created in ``path``.
    """

	global omega_rescaled
	
	t0=time()
	Ne=len(omega[0])
	Nl=len(laser)

	if states==None: states=range(1,Ne+1)

	#We make some checks
	for i in range(Ne):
		for j in range(Ne):
			b1=not ('.' in str(omega[i][j]) or 'e' in str(omega[i][j]))
			if b1: raise ValueError,'omega must be composed of floating point numbers.'
			b2=not ('.' in str(gamma[i][j]) or 'e' in str(gamma[i][j]))
			if b2: raise ValueError,'gamma must be composed of floating point numbers.'

	#We rescale the frequencies as requested.
	#~ if Omega != None:
		#~ omega_rescaled=[[omega[i][j]/Omega for j in range(Ne)] for i in range(Ne)]
		#~ #gamma=[[gamma[i][j]/Omega for j in range(Ne)] for i in range(Ne)]
	#~ else:
		#~ omega_rescaled=omega[:]
	omega_rescaled=omega[:]
		

	#We determine wether it is possible to eliminate explicit time-dependance
	theta=find_phase_transformation(Ne,Nl,r,Lij)
	#We find the detunings if required

	#We construct the correspondence i <-> I between degenerate and non-degenerate indices.
	i_d,I_nd,Nnd=calculate_iI_correspondence(omega)
	#We get wich transitions each laser induces
	detunings,detuningsij=laser_detunings(Lij,Nl,i_d,I_nd,Nnd)
	
	#We get how many transitions each  laser induces
	detuning_indices=[len(detunings[i]) for i in range(Nl)]
	
	#The number of detunings
	Nd=sum([len(detunings[l]) for l in range(Nl)])
	combinations=detuning_combinations(detuning_indices)	

	code0='''program evolution_rk4
	implicit none
	complex*16, dimension('''+str(Ne*(Ne+1)/2-1)+') :: x,k1,k2,k3,k4\n'
	code0+='''    real*8 :: dt,t,ddelta,delta,delta0
	integer :: i,j,n,ldelta,ndelta,detuning_index,n_aprox,n_mod

	logical :: print_steps,run_spectrum\n'''
	code0+='    real*8, dimension('+str(Nl)+') :: E0,detuning_knob\n'
	code0+='    real*8, dimension('+str(Nd)+') :: detuning\n\n'

	code0+="    open(unit=1,file='"+path+name+".dat',status='unknown')\n\n"
	code0+='    n_aprox=1500\n'
	code0+='    !We load the parameters\n'
	code0+="    open(unit=2,file='"+path+name+"_params.dat',status='unknown')\n" 
	code0+='''    read(2,*) n
    read(2,*) dt
    read(2,*) print_steps
    read(2,*) x
    read(2,*) E0\n'''

	code0+='    read(2,*) detuning_knob\n'
	code0+='    read(2,*) run_spectrum\n\n'
	code0+='''    if (run_spectrum) then
		read(2,*) ldelta
		read(2,*) ndelta
		read(2,*) ddelta
		close(2)
		delta0=detuning_knob(ldelta)
		n_mod=ndelta/n_aprox
    else
		ldelta=1; ndelta=1; ddelta=0; delta=0
		close(2)
		n_mod=n/n_aprox
    end if
    if (n_mod==0) n_mod=1\n\n\n'''


	#We add the code to caculate all the initial detunings for each laser.
	code0+='	!We calculate the initial detunings.\n'
	#We find the minimal frequency corresponding to each laser.
	omega_min,omega_min_indices=find_omega_min(omega_rescaled,Nl,detuningsij,i_d,I_nd)
	
	det_index=1
	for l in range(Nl):
		omega0=omega_min[l]
		for p in detuningsij[l]:
			code0+='	detuning('+str(det_index)+')='
			code0+=format_double(omega0-omega_rescaled[p[0]][p[1]])+'+detuning_knob('+str(l+1)+')\n'
			det_index+=1
	code0+='\n'

	code0+='''	t=0
	if (.not. run_spectrum) WRITE(1,*) t,real(x),imag(x('''+str(Ne)+''':))\n'''
	code0+='''	!We start the detuning variation\n'''
	code0+='	delta=detuning_knob(ldelta)\n'
	code0+='''    do j=1,ndelta
		!We run the Runge Kutta method
		t=0.0
		do i=1,n-1\n'''

	code0+='            call f(x          , t       , k1,   E0, detuning, detuning_knob)\n'
	code0+='            call f(x+0.5*k1*dt, t+dt*0.5, k2,   E0, detuning, detuning_knob)\n'
	code0+='            call f(x+0.5*k2*dt, t+dt*0.5, k3,   E0, detuning, detuning_knob)\n'
	code0+='            call f(x    +k3*dt, t+dt    , k4,   E0, detuning, detuning_knob)\n'

	code0+='''			x= x+(k1+2*k2+2*k3+k4)*dt/6
			if (print_steps.and. .not. run_spectrum) print*,'t=',t,'delta=',delta
			t= t+ dt
			
			if (isnan(real(x(1)))) stop 1
			if (.not. run_spectrum .and. mod(i,n_mod)==0) WRITE(1,*) t,real(x),imag(x('''+str(Ne)+''':))
		end do
		if (print_steps) print*, 'delta=',delta,'percentage=',100*(delta-delta0)/(ddelta*ndelta)
		
		!We recalculate the detunings
		if (run_spectrum) then
			if (mod(j,n_mod)==0) WRITE(1,*) delta,real(x),imag(x('''+str(Ne)+''':))
			delta=delta+ddelta
			detuning_knob(ldelta)=detuning_knob(ldelta)+ddelta\n'''
	


	#We add the code to caculate all detunings for each laser
	#This way of assigining a global index ll to the detunings ammounts to
	#   ll=   number_of_previous_detunings 
	#       + number_of_detuning_ordered_by_row_and_from_left_to_right_column
	#like this
	#->
	#-> ->
	#-> -> ->
	#for each l
	
	#We find the minimal frequency corresponding to each laser		
	omega_min,omega_min_indices=find_omega_min(omega_rescaled,Nl,detuningsij,i_d,I_nd)
	
	det_index=1
	for l in range(Nl):
		omega0=omega_min[l]
		code0+='			if (ldelta=='+str(l+1)+') then\n'
		for p in detuningsij[l]:
			code0+='				detuning('+str(det_index)+')=detuning('+str(det_index)+')'
			#code0+='+('+str(omega0-omega_rescaled[p[0]][p[1]])+'+ddelta\n'
			code0+='+ddelta\n'
			det_index+=1
		code0+='			end if\n'
	
	code0+='''		end if
	
	
	end do
	
    close(1)
    
end program\n\n'''


	code0+='subroutine f(x,t,y,    E0, detuning,detuning_knob)\n'
	
	code0+='''    implicit none
    real*8, intent(in) :: t\n'''
	code0+='    complex*16, dimension('+str(Ne*(Ne+1)/2-1)+'), intent(in)  :: x\n'
	code0+='    complex*16, dimension('+str(Ne*(Ne+1)/2-1)+'), intent(out) :: y\n'
	code0+='    real*8, dimension('+str(Nl)+'), intent(in) :: E0,detuning_knob\n'
	code0+='    real*8, dimension('+str(Nd)+'), intent(in) :: detuning\n\n'

	code0+='    complex*16 :: I,fact,aux\n'
	code0+='    real*8 :: rho11\n\n'
	code0+='    I=(0,1D0)\n'
	
	#We establish the scaling of the equations
	#~ if Omega==None:
		#~ h  =1.054571726e-34; e=1.602176565e-19
		#~ code0+='    fact=I*'+str(e/h)+'\n'
	#~ else:
		#~ #code0+='    fact=I*'+str(float(Omega/sqrt(2)))+'\n'
		#~ code0+='    fact=I*'+str(float(1/sqrt(2)))+'\n'
		#~ #code0+='    fact=I*'+str(float(1/(sqrt(2)*Omega)))+'\n'
	code0+='    fact=I*'+format_double(float(1/sqrt(2)))+'\n'
	
	#We give the code to calculate rho11
	code0+='    rho11=1\n'
	for i in range(1,Ne):
		code0+='    rho11=rho11 -x('+str(i)+')\n'
	code0+='\n\n'
	
	
	####################################################################
	#We produce the code for the first order equations.
	####################################################################
	if len(theta)>0:
		code=''
		for mu in range(1,Ne*(Ne+1)/2):
			i,j,s=IJ(mu,Ne)
			#print 'ecuacion mu=',mu,',i,j=',i,j
			eqmu='    y('+str(mu)+')= 0\n'
			
			####################################################################
			#We add the terms associated with the effective hamiltonian
			#other than those associated with the phase transformation.
			
			for k in range(1,Ne+1):
				#Case 1
				if k>=j:
					for l in Lij[i-1][k-1]:
						if k>i:
							#print 'E0^',l,-1,'r',i,k,'rho',k,j,'case 1.1'
							eqmu+=add_line(Ne,mu,'+',laser,l,-1,  r,i,k,  k,j)
						elif k<i:
							#print 'E0^',l, 1,'r',i,k,'rho',k,j,'case 1.2'
							eqmu+=add_line(Ne,mu,'+',laser,l, 1,  r,i,k,  k,j)
							
				#Case 2
				elif k<j:
					for l in Lij[i-1][k-1]:
						if k>i:
							#print 'E0^',l,-1,'r',i,k,'rhoa',j,k,'case 2.1'
							eqmu+=add_line(Ne,mu,'+',laser,l,-1,  r,i,k,  j,k,True)
						elif k<i:
							#print 'E0^',l, 1,'r',i,k,'rhoa',j,k,'case 2.2'
							eqmu+=add_line(Ne,mu,'+',laser,l, 1,  r,i,k,  j,k,True)
				#Case 3
				if k<=i:
					for l in Lij[k-1][j-1]:
						if k<j:
							#print 'E0^',l,-1,'r',k,j,'rho',i,k,'case 3.1'
							eqmu+=add_line(Ne,mu,'-',laser,l,-1,  r,k,j,  i,k)
						elif k>j:
							#print 'E0^',l, 1,'r',k,j,'rho',i,k,'case 3.2'
							eqmu+=add_line(Ne,mu,'-',laser,l, 1,  r,k,j,  i,k)
				#Case 4
				elif k>i:
					for l in Lij[k-1][j-1]:
						if k<j:
							#print 'E0^',l,-1,'r',k,j,'rhoa',k,i,'case 4.1'
							eqmu+=add_line(Ne,mu,'-',laser,l,-1,  r,k,j,  k,i,True)
						elif k>j:
							#print 'E0^',l, 1,'r',k,j,'rhoa',k,i,'case 4.2'
							eqmu+=add_line(Ne,mu,'-',laser,l, 1,  r,k,j,  k,i,True)

			eqmu+='    y('+str(mu)+')=y('+str(mu)+')*fact\n'
			
			####################################################################
			#We add the terms associated with the phase transformation.

			extra=Theta(i,j,theta,omega_rescaled,omega_min,detunings,detuningsij,
						combinations,detuning_indices,Lij,i_d,I_nd,Nnd,
						verbose=verbose,states=states)

			if extra!='':
				eqmu+='    y('+str(mu)+')=y('+str(mu)+') + I*('+extra+')*x('+str(mu)+')\n'

			####################################################################
			
			#~ if i==j:
				#~ for k in range(1,Ne+1):
					#~ if k < i:
						#~ muii=Mu(i,i,s=1,N=Ne)
						#~ eqmu+='    y('+str(mu)+')= y('+str(mu)+') - ('+format_double(gamma[i-1][k-1])+')*x('+str(muii)+')\n'
					#~ elif k > i:
						#~ mukk=Mu(k,k,s=1,N=Ne)
						#~ eqmu+='    y('+str(mu)+')= y('+str(mu)+') - ('+format_double(gamma[i-1][k-1])+')*x('+str(mukk)+')\n'
				#~ eqmu+='\n'
			#~ else:
				#~ eqmu+='    y('+str(mu)+')= y('+str(mu)+') - ('+format_double(gamma[i-1][j-1]/2)+')*x('+str(mu)+')\n'
#~ 

			####################################################################
			code+=eqmu+'\n'

		#We add the terms associated with spontaneous decay.		
		#First for populations.
		for i in range(2,Ne+1):
			mu=Mu(i,i,1,Ne)
			for k in range(1,Ne+1):
				gams=0
				if k<i:
					gams+=gamma[i-1][k-1]
				elif k>i:
					nu=Mu(k,k,1,Ne)
					ga=gamma[i-1][k-1]
					if ga != 0:
						code+='    y('+str(mu)+')=y('+str(mu)+')'
						code+='-('+format_double(ga)+')*x('+str(nu)+')\n'
				if gams!=0:
					code+='    y('+str(mu)+')=y('+str(mu)+')'
					code+='-('+format_double(gams)+')*x('+str(mu)+')\n'

		#And now for coherences	
		for i in range(1,Ne+1):
			for j in range(1,i):
				gams=gamma[i-1][j-1]/2
				if gams!=0:
					for a in range(i+1,Ne+1):
						mu=Mu(a,i,+1,Ne)
						code+='    y('+str(mu)+')=y('+str(mu)+')'
						code+='-('+format_double(gams)+')*x('+str(mu)+')\n'

						#~ mu=Mu(a,i,-1,Ne)
						#~ code+='    y('+str(mu)+')=y('+str(mu)+')'
						#~ code+='-('+format_double(gams)+')*x('+str(mu)+')\n'

					for b in range(1,i):
						mu=Mu(i,b,+1,Ne)
						code+='    y('+str(mu)+')=y('+str(mu)+')'
						code+='-('+format_double(gams)+')*x('+str(mu)+')\n'

						#~ mu=Mu(i,b,-1,Ne)
						#~ code+='    y('+str(mu)+')=y('+str(mu)+')'
						#~ code+='-('+format_double(gams)+')*x('+str(mu)+')\n'




		
		####################################################################
		####################################################################
		####################################################################
		#code+='    y=y/'+str(Omega)+'\n'
		
		f=file(path+name+'.f90','w')
		code=code0+code+'end subroutine\n'
		f.write(code)
		f.close()
				
		return time()-t0
	else:
		print 'There was no phase transformation capable of eliminating explicit time dependance.'

def run_rk4(path,name,E0,laser_frequencies,  N_iter,dt,N_states,
				spectrum_of_laser=None,N_delta=None,frequency_step=None,frequency_end=None,
				rho0=None,print_steps=False,integrate=False,
				save_systems=False):
	"""This function runs the Runge-Kutta method compiled in path+name..."""

	t0=time()
	params =str(N_iter)+'\n'
	params+=str(dt)+'\n'
	#We give the flag on wether to print each time step.
	if print_steps:
		params+='.true.\n'
	else:
		params+='.false.\n'
	
	#We give the initial value of rho
	N_vars=N_states*(N_states+1)/2-1
	if rho0==None:
		params+=''.join(['(0.0,0.0) ' for i in range(N_vars)])
	elif len(rho0)==N_states-1:
		if sage_included:
			params+=''.join(['('+str(real(i))+','+str(imag(i))+') ' for i in rho0])
		else:
			params+=''.join(['('+str(i.real)+','+str(i.imag)+') ' for i in rho0])
		params+=''.join(['(0.0,0.0) ' for i in range(N_vars-N_states+1)])
	elif len(rho0)==N_vars:
		params+=''.join(['('+str(real(i))+','+str(imag(i))+') ' for i in rho0])

	params+='\n'
	#We give the amplitude of the electrical fields.
	params+=''.join([str(i)+' ' for i in E0])+'\n'
	
	#We give the detuning of each laser (taken from the lowest frequency transition).
	params+=''.join([str(i)+' ' for i in laser_frequencies])+'\n'

	#We give the flag on wether to calculate spectrums or time evolution.
	if spectrum_of_laser==None:
		params+='.false.'
	else:
		if frequency_end !=None:
			if frequency_step !=None:
				raise ValueError,'both frequency_end and frequency_step were specified.'
			if N_delta==1:
				frequency_step=0.0
			else:
				frequency_step=(frequency_end-laser_frequencies[spectrum_of_laser-1])/(N_delta-1)
			
		#frequency_step=frequency_end
		params+='.true.\n'
		params+=str(spectrum_of_laser)+'\n'
		params+=str(N_delta)+'\n'
		params+=str(frequency_step)
	#print params
	
	f=file(path+name+'_params.dat','w')
	f.write(params)
	f.close()
	
	os.system(path+name)
	
	return time()-t0



