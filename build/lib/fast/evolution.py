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
	from math import sqrt,log,pi
	from misc import Mu,IJ,find_phase_transformation, format_double
	from misc import write_equations_code
	from rk4 import write_rk4, run_rk4

	from stationary import analyze_zeros
	from time import time
	import os
else:
	from time import time

def write_evolution(path,name,laser,omega,gamma,r,Lij,states=None,
                    excluded_mu=[],rk4=False,verbose=1):
	
	if rk4:
		return write_rk4(path,name,laser,omega,gamma,r,Lij,
						 states=states,verbose=verbose)

	t0=time()
	Ne=len(omega[0])
	Nl=len(laser)
	N_excluded_mu=len(excluded_mu)
	Nrho=Ne**2-1
	print_times=False

	dummy=write_equations_code(path,name,laser,omega,gamma,r,Lij,
				states=states,excluded_mu=excluded_mu,verbose=verbose)
	
	code,Nd,row_check,col_check,rhs_check,Ne,N_excluded_mu,states,omega_min,detuningsij,omega_rescaled=dummy

	
	code0='''program evolution_diagonalization
	implicit none\n'''
	code0+='''    real*8 :: dt,ddelta,delta0
	real*8, allocatable, dimension(:) :: t,delta
	integer :: i,j,mu,n,ldelta,ndelta,detuning_index,n_aprox,n_mod,info

	logical :: print_steps,run_spectrum,save_systems,save_eigenvalues,integrate,use_netcdf\n'''
	code0+='    real*8, dimension('+str(Nl)+') :: E0,detuning_knob\n'
	code0+='    real*8, dimension('+str(Nd)+') :: detuning\n\n'
	
	code0+="	complex*16, dimension("+str(Nrho)+") :: r_amp,rho0\n"
	code0+="	complex*16, allocatable, dimension(:,:) :: rho,rho_spectrum\n"
	code0+="	real*8, dimension("+str(Nrho)+") :: rho_inf\n"
	code0+="	complex*16, dimension("+str(Nrho)+","+str(Nrho)+") :: U\n"
	code0+="	complex*16, dimension("+str(Nrho)+") :: lam\n\n"

	if print_times: code0+="	real*8 :: t7,t8\n\n"

	code0+='    !We load the parameters\n'
	code0+='    n_aprox=1500\n'
	# We break the path name into several lines if it is needed.
	long_line="    open(unit=2,file='"+path+name+"_params.dat',status='unknown')\n" 
	if len(long_line)>=72:
		long_line=long_line[:72]+"&\n&"+long_line[72:]
	code0+=long_line

	code0+='''    read(2,*) n
    read(2,*) dt
    read(2,*) print_steps
    read(2,*) rho0
    read(2,*) E0\n'''

	code0+='    read(2,*) detuning_knob\n'
	code0+='    read(2,*) run_spectrum\n'
	code0+='    read(2,*) save_systems\n'
	code0+='    read(2,*) save_eigenvalues\n'
	code0+='    read(2,*) use_netcdf\n'
	code0+='    read(2,*) integrate\n\n'
	
	code0+='''    if (run_spectrum) then
		read(2,*) ldelta
		read(2,*) ndelta
		read(2,*) ddelta
		
		
		allocate(rho(n,'''+str(Nrho)+'''),stat=info)
		allocate(rho_spectrum(ndelta,'''+str(Nrho)+'''),stat=info)
		rho(1,:)=rho0
    else
		ldelta=1; ndelta=1; ddelta=0
		
		allocate(rho(n,'''+str(Nrho)+'''),stat=info)
		rho(1,:)=rho0
    end if
    close(2)
    
    delta0=detuning_knob(ldelta)
    
    !We build the time and delta axis
    allocate(t(n),stat=info)
    allocate(delta(ndelta),stat=info)
    do i=1,n
		t(i)=(i-1)*dt
    end do

    do i=1,ndelta
		delta(i)=delta0+(i-1)*ddelta
    end do
        
    if (n_mod==0) n_mod=1\n\n'''

	#code0+='    n_aprox=1500\n'
	code0+="""    if (save_eigenvalues) then\n"""
    
	# We break the path name into several lines if it is needed.
	long_line="""open(unit=3,file='"""+path+name+"""_eigenvalues.dat',status='unknown')\n"""
	if len(long_line)>=72:
		long_line=long_line[:72]+"&\n&"+long_line[72:]
	code0+=long_line
	
	code0+="""	end if\n\n"""

	#~ #We add the code to caculate all the initial detunings for each laser.
	#~ code0+='	!We calculate the initial detunings.\n'
	#~ #We find the minimal frequency corresponding to each laser.
	#~ #omega_min,omega_min_indices=find_omega_min(omega_rescaled,Nl,detuningsij)
	#~ 
	#~ det_index=1
	#~ for l in range(Nl):
		#~ omega0=omega_min[l]
		#~ for p in detuningsij[l]:
			#~ code0+='	detuning('+str(det_index)+')='
			#~ code0+=format_double(omega0-omega_rescaled[p[0]][p[1]])+'+detuning_knob('+str(l+1)+')\n'
			#~ det_index+=1
	#~ code0+='\n'

	code0+='''	!We start the detuning variation, which might consist of a single detuning.\n'''

	code0+='''    do j=1,ndelta
		!We calculate the solution for specific E0, detuning_knob, and initial rho,
		!this returns U, r_amp, rho_inf, lam.
		detuning_knob(ldelta)=delta(j)

		call solve(E0,detuning_knob, rho0,U,r_amp,rho_inf,lam,save_systems)
		
		if (save_eigenvalues) then
			write(3,*) delta,real(lam)
			write(3,*) delta,imag(lam)
		end if\n\n'''
	
	if print_times:
		code0+="		call cpu_time(t7)\n"
	code0+='''
		!We calculate the time evolution with the solution just computed.
		do i=2,n\n'''


	code0+='''
			if (print_steps.and. .not. run_spectrum) print*,'t=',t,'delta=',delta
			
			do mu=1,'''+str(Nrho)+'''
				rho(i,mu)=r_amp(mu)*cdexp(lam(mu)*t(i))
			end do
			rho(i,:)= matmul(U , rho(i,:)) + rho_inf
			
			if (isnan(real(rho(i,1)))) stop 1
			!if (.not. run_spectrum ) WRITE(1,*) t,real(rho(i,:))
		end do\n'''
	if print_times:
		code0+='''
		call cpu_time(t8)
		print*,'time to calculate the n points of time:',t8-t7\n'''
	
	code0+='''
		if (print_steps) print*, 'delta=',delta,'percentage=',100*(delta-delta0)/(ddelta*ndelta)
		
		!We save the time evolution.
		if (.not. run_spectrum ) then'''
	#We decide whether to use netcdf.
	from config import use_netcdf
	
	if use_netcdf:
		code0+='''
			call save_matrix_and_vector("'''+path+name+'''.nc",n,'''+str(Ne**2-1)+''',real(rho),t)
		else
			rho_spectrum(j,:)=rho(n,:)
        end if
'''
	else:
		code0+="\n"
		long_line="			open(unit=1,file='"+path+name+".dat',status='unknown')\n" 
		if len(long_line)>=72:
			long_line=long_line[:72]+"&\n&"+long_line[72:]
		code0+=long_line

		code0+='''
			do i=1,n
				WRITE(1,*) t(i),real(rho(i,:))
			end do
			close(1)
		else
			rho_spectrum(j,:)=rho(n,:)
		end if
		'''
	
	code0+='''
	end do

	!We save a spectrum.
	if (run_spectrum) then\n'''
	if use_netcdf:
		code0+='''
        call save_matrix_and_vector("'''+path+name+'''.nc",ndelta,'''+str(Ne**2-1)+''',real(rho_spectrum),delta)'''
	else:
		code0+='''
		open(unit=1,file="'''+path+name+'''.dat",status='unknown')
		do i=1,ndelta
			WRITE(1,*) delta(i),real(rho_spectrum(i,:))
		end do
		close(1)'''
	
	code0+='''
	end if
	    
end program
'''
	
	if use_netcdf:
		code0+='''
subroutine save_matrix_and_vector(file_name,m,n,matrix,vector)
	use netcdf
	implicit none
	character (len = *), intent(in) :: file_name
	integer, intent(in) :: m,n
	real*8, dimension(m,n), intent(in) :: matrix
	real*8, dimension(m), intent(in) :: vector
	
	integer :: ncid, varid_matrix, varid_vector, dimids(2)
	integer :: x_dimid, y_dimid
	
	call check( nf90_create(file_name, NF90_CLOBBER, ncid) )
	call check( nf90_def_dim(ncid, "x", n, x_dimid) )
	call check( nf90_def_dim(ncid, "y", m, y_dimid) )
	dimids =  (/ y_dimid, x_dimid /)
	call check( nf90_def_var(ncid, "matrix", NF90_DOUBLE, dimids,  varid_matrix) )
	call check( nf90_def_var(ncid, "vector", NF90_DOUBLE, y_dimid, varid_vector) )
	call check( nf90_enddef(ncid) )
	
	call check( nf90_put_var(ncid, varid_matrix, matrix) )
	call check( nf90_put_var(ncid, varid_vector, vector) )
	call check( nf90_close(ncid) )
end subroutine

subroutine check(status)
	use netcdf
	implicit none
	integer, intent ( in) :: status

	if(status /= nf90_noerr) then 
		print *, trim(nf90_strerror(status))
		stop "Stopped"
	end if
end subroutine check  

\n\n'''
		
	code0+=r"""subroutine solve(E0,detuning_knob,rho0,U,r_amp,rho_inf,lam,save_systems)
	implicit none
    integer :: i,j,k,mu,nu,alpha	
	real*8, dimension("""+str(Nl)+"""), intent(in) :: E0,detuning_knob
	complex*16, dimension("""+str(Nrho)+""",1), intent(in) :: rho0
	complex*16, dimension("""+str(Nrho)+""","""+str(Nrho)+"""), intent(out) :: U
	complex*16, dimension("""+str(Nrho)+""",1), intent(out) :: r_amp
	real*8, dimension("""+str(Nrho)+""",1), intent(out) :: rho_inf
	complex*16, dimension("""+str(Nrho)+"""), intent(out) :: lam
	logical, intent(in) :: save_systems"""
	

	code0+=r"""
	real*8, dimension("""+str(Nd)+""") :: detuning
	real*8, dimension("""+str(Nrho)+""","""+str(Nrho)+""") :: A
	complex*16, dimension("""+str(Nrho)+""",1) :: B
	integer, dimension("""+str(Nrho)+""") :: IPIV

    integer :: Nrho
    real*8, dimension("""+str(Nrho)+""") :: WR,WI
    real*8, dimension(1,"""+str(Nrho)+""") :: VL
    real*8, dimension("""+str(Nrho)+""","""+str(Nrho)+""") :: VR,VRT
    real*8, dimension("""+str(6*Nrho)+""") :: WORK
    integer :: INFO
    
    complex*16 :: II
    complex*16, dimension("""+str(Nrho)+""","""+str(Nrho)+""") :: Ui,Us,Uis
    logical :: band\n"""
	
	if print_times:
		code0+="""\n    real*8 :: t1,t2,t3,t4,t5,t6,t7,t8\n"""
    
	code0+="""
	Nrho="""+str(Nrho)+"""
	A=0
	B=0	
"""

	####################################################################
	#~ #We make the program print A and B
	#~ code+="""	do j=1,"""+str(Ne**2-1)+"""
		#~ print*, A(j,:)
	#~ end do
		#~ print*,
		#~ do j=1,"""+str(Ne**2-1)+"""
			#~ print*,B(j,:)
		#~ end do
		#~ print* \n"""
	####################################################################
	
	if print_times:
		code+="""	call cpu_time(t2)
	print*,'time to form A:',t2-t1\n\n"""
	
	code+="	if (save_systems) then\n"
	
	long_line="    open(unit=4,file='"+path+name+"_AB.dat',status='unknown')\n" 
	if len(long_line)>=72:
		long_line=long_line[:72]+"&\n&"+long_line[72:]
	code+=long_line

	code+="		do j=1,"+str(Ne**2-1-N_excluded_mu)+'\n'
	code+="			write(4,*) A(j,:),B(j,1)\n"
	code+="		end do\n"
	code+="		close(4)\n"
	code+="	end if\n\n"
	
	#code+='	call dgesv('+str(Ne**2-1-N_excluded_mu)+', 1, A, '+str(Ne**2-1-N_excluded_mu)+', IPIV, B, '+str(Ne**2-1-N_excluded_mu)+', INFO)\n'
	#code+="	print*,'INFO',INFO\n"
	#code+="""	if (INFO>0) B=-11\n"""
	#code+="""	if (INFO>0) print*, 'For frequencies',detuning_knob,'The system could not be solved, exit code:',INFO\n"""
	#code+="""	if (INFO>0) then\n"""
	#code+="""		do j="""+str(Ne**2-1)+"""\n"""
	#code+="""			print*,\n"""	
	#code+="""		end do\n"""
	#code+="""	endif\n"""
	####################################################################
	#We check the rows and columns searching for rows of zeros and
	#columns of zeros.
	excluded_mu=analyze_zeros(row_check,col_check,rhs_check,Ne,N_excluded_mu,states)

	if excluded_mu!=[]:
		t0=time()-t0
		t_extra=write_stationary(path,name,laser,omega,gamma,r,Lij,
				use_symbolic_phase_transformation=True,states=states,excluded_mu=excluded_mu)
		
		return t_extra+t0
	else:		
		f=file(path+name+'.f90','w')
		
		code+="""
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!We begin the solution process.
    call DGEEV( 'N', 'V', Nrho, A, Nrho, WR, WI, VL, 1, VR, Nrho, WORK, 6*Nrho, INFO )\n"""
    
		if print_times:
			code+="""\n	call cpu_time(t3)
	print*,'time to diagonalize:',t3-t2\n\n"""
	
		code+="""    II=(0,1)
    !We build the matrix U and the eigenvalues from the information given by DGEEV.
    band=.false.
    do i=1,Nrho
		lam(i)=WR(i)+II*WI(i)
		!print*,i,lam(i),WR(i),WI(i)
		if (WI(i)==0.0d0       ) then
			U(:,i)=VR(:,i)
		elseif (band) then
			U(:,i)=VR(:,i-1)-II*VR(:,i)
			band=.false.
		elseif     (WI(i)==-WI(i+1)) then
			U(:,i)=VR(:,i)+II*VR(:,i+1)
			band=.true.
		else
			print*,'ERROR: I could not determine wether the ith eigenvalue was real or a complex pair.'
			INFO=0
			print*,1/INFO
		end if
    end do\n"""
    
		if print_times:
			code+="""\n	call cpu_time(t4)
	print*,'time to form U and lam:',t4-t3\n"""

    
		code+="""
	!We calculate the inverse of U.
    Us=U; Ui=0

	call ZGETRF( Nrho, Nrho, Us, Nrho, IPIV, INFO )
	if (INFO/=0) then
		INFO=0
		PRINT*,1/INFO
	end if
    call ZGETRI( Nrho, Us, Nrho, IPIV, WORK, Nrho, INFO )
	if (INFO/=0) then
		INFO=0
		PRINT*,1/INFO
	end if
	Ui=Us\n"""
		if print_times:
			code+="""
	call cpu_time(t5)
	print*,'time to form U and lam:',t5-t4\n"""

		code+="""	!We calculate d (which will be made by successive modifications of b)
	b=matmul(Ui,b)
	do mu=1,Nrho
		b(mu,1)=b(mu,1)/lam(mu)
	end do
"""

		if print_times:
			code+="""
	call cpu_time(t6)
	print*,'time to calculate d:',t6-t5\n"""


		code+="""	!We calculate the vector r
	r_amp=matmul(Ui,rho0) - b
"""
		if print_times:
			code+="""
	call cpu_time(t7)
	print*,'time to calculate r_amp:',t7-t6\n"""

		code+="""	!We calculate rho_inf
	rho_inf=matmul(U,b)
"""
		if print_times:
			code+="""
	call cpu_time(t8)
	print*,'time to calculate rho_fin:',t8-t7
	print*,'total time:',t8-t1\n"""
	

		#~ code+="""
	#~ !We calculate rho_inf, and rho_breve
	#~ rho_inf=0
	#~ rho_breve=0
	#~ do alpha=1,Nrho
		#~ do mu=1,Nrho
			#~ do nu=1,Nrho
				#~ rho_inf(alpha)=rho_inf(alpha) + U(alpha,mu) * Ui(mu,nu) * b(nu,1)/lam(mu)
				#~ rho_breve(alpha,mu) = rho_breve(alpha,mu) + U(alpha,mu)*Ui(mu,nu)*(rho0(nu) - b(nu,1)/lam(mu) )
			#~ end do
		#~ end do
	#~ end do\n"""
		
		code=code0+code+'end subroutine\n'
		f.write(code)
		f.close()
		
		return time()-t0

def run_evolution(path,name,E0,laser_frequencies,  N_iter,dt,N_states,
				spectrum_of_laser=None,N_delta=None,frequency_step=None,frequency_end=None,
				rho0=None,print_steps=False,
				integrate=False,
				save_systems=False,save_eigenvalues=False,rk4=False,use_netcdf=True):
	"""This function runs the Runge-Kutta method compiled in path+name..."""
	def py2f_bool(bool_var):
		if bool_var:
			return ".true.\n"
		else:
			return ".false.\n"
	
	if rk4:
		return run_rk4(path,name,E0,laser_frequencies,  N_iter,dt,N_states,
				spectrum_of_laser=spectrum_of_laser,N_delta=N_delta,
				frequency_step=frequency_step, frequency_end=frequency_end,
				rho0=rho0,print_steps=print_steps,
				integrate=integrate,save_systems=save_systems)
	
	t0=time()
	params =str(N_iter)+'\n'
	params+=str(dt)+'\n'
	#We give the flag on wether to print each time step.
	params+=py2f_bool(print_steps)
	
	#We give the initial value of rho
	N_vars=N_states*(N_states+1)/2-1
	if rho0==None:
		params+=''.join(['(0.0,0.0) ' for i in range(N_states**2-1)])
	elif len(rho0)==N_states-1:
		if sage_included:
			params+=''.join(['('+str(real(i))+','+str(imag(i))+') ' for i in rho0])
		else:
			params+=''.join(['('+str(i.real)+','+str(i.imag)+') ' for i in rho0])
		params+=''.join(['(0.0,0.0) ' for i in range( N_states**2 -N_states)])
	elif len(rho0)==N_vars:
		params+=''.join(['('+str(real(i))+','+str(imag(i))+') ' for i in rho0])
	else:
		raise ValueError,'rho0 had an invalid number of elements.'

	params+='\n'
	#We give the amplitude of the electrical fields.
	params+=''.join([str(i)+' ' for i in E0])+'\n'
	
	#We give the detuning of each laser (taken from the lowest frequency transition).
	params+=''.join([str(i)+' ' for i in laser_frequencies])+'\n'

	#We give the flag on wether to calculate spectrums or time evolution.
	if spectrum_of_laser==None:
		params+='.false.\n'
		params+=py2f_bool(save_systems)
		params+=py2f_bool(save_eigenvalues)
		params+=py2f_bool(use_netcdf)
		params+=py2f_bool(integrate)

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
		params+=py2f_bool(save_systems)
		params+=py2f_bool(save_eigenvalues)
		params+=py2f_bool(use_netcdf)
		params+=py2f_bool(integrate)

		params+=str(spectrum_of_laser)+'\n'
		params+=str(N_delta)+'\n'
		params+=str(frequency_step)
	#print params
	
	f=file(path+name+'_params.dat','w')
	f.write(params)
	f.close()
	
	os.system(path+name)
	return time()-t0

def get_eigenvalues(path,name):
	f=file(path+name+'_eigenvalues.dat')
	d=f.readlines()
	f.close()
	re=[[ float(n) for n in d[2*i][:-1].split()] for i in range(len(d)/2)]
	im=[[ float(n) for n in d[2*i+1][:-1].split()] for i in range(len(d)/2)]

	return re,im

def characteristic_times(path,name,Omega=1):
	r'''This function can be called after calling ``run_diagonalization`` if the option ``save_eigenvalues`` 
is set to ``True``. It will return the oscillation periods, and the shortest and half lives.
The results are lists ordered as:
	``[detunings, oscillation_periods_1, oscillation_periods_N**2-1, half_lives_1, half_lives_N**2-1]``.
'''
	re,im=get_eigenvalues(path,name)
	
	log12=log(0.5)
	Nr=len(re[0]); Nd=len(re)
	half_lives=[]; oscillation_periods=[]
	for i in range(Nr):
		col_half=[]
		col_osci=[]
		for j in range(Nd):
			if re[j][i]>=0.0 and i!=0:
				raise ValueError,'an eigenvalue was greater or equall to zero:'+str(re[j][i])+'.'
			else:
				col_half+=[log12/re[j][i]/Omega]
			
			if im[j][i]==0.0 and i!=0:
				col_osci+=[float('inf')]
			else:
				col_osci+=[abs(2*pi/im[j][i])/Omega]
		half_lives+=[col_half]
		oscillation_periods+=[col_osci]
	return oscillation_periods+half_lives[1:]

def analyze_eigenvalues(path,name,Ne):
	r'''This function can be called after calling ``run_diagonalization`` if the option ``save_eigenvalues`` 
is set to ``True``. It will return the shortest and longest oscillation period, and the shortest and 
longest half life. The results are lists ordered as:
	``[detunings, shortest_oscillations, longest_oscillations, shortest_half_lifes, longest_half_lifes]``.

The shortest oscillation period can be interpreted as a suggested size for the time step needed for a 
fine-grained picture of time evolution (an order of magnitude smaller is likely to be optimal).

The longest half life can be interpreted as a suggested time for the stationary state to be reached
(an order of magnitude larger is likely to be optimal).
'''
	times=characteristic_times(path,name)

	min_osci=min([min(times[1+i]) for i in range(Ne**2-1)])
	max_osci=max([max(times[1+i]) for i in range(Ne**2-1)])
	min_half=min([min(times[1+Ne**2-1+i]) for i in range(Ne**2-1)])
	max_half=max([max(times[1+Ne**2-1+i]) for i in range(Ne**2-1)])
	
	return min_osci,max_osci,min_half,max_half

from matplotlib import pyplot
