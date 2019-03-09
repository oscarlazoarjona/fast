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

#from all import sage_included
sage_included = 'sage' in globals().keys()
if not sage_included:
	from misc import write_equations_code
	from time import time
	import os
else:
	from time import time

def analyze_zeros(row_check,col_check,rhs_check,Ne,N_excluded_mu,states):
	band_zeros=False
	zero_row=[]; zero_col=[]
	for mu in range(Ne**2-1-N_excluded_mu):
		if not row_check[mu]:
			band_zeros=True
			zero_row+=[mu+1]
		if not col_check[mu]:
			band_zeros=True
			zero_col+=[mu+1]

	excluded_mu=[]
	if band_zeros:
		#If there are disconnected mu we must exclude them.
		if zero_row==zero_col:
			excluded_mu=zero_row
			#We calculate the mu that must be excluded
			for mu in zero_row:
				i,j,s=IJ(mu,Ne)
				rhs=rhs_check[mu-1]
				if rhs:
					s='The rhs of equation '+str(mu)+' is non zero, but the lhs is zero.'
					raise ValueError,s
				else: rhs=0
				if i==j:
					print 'The population of state',states[i-1],'has equation 0=0.'
				elif s==1:
					print 'The real part of the coherence between states',states[i-1],'and',states[j-1],'has equation 0 =',rhs
				else:
					print 'The imaginary part of the coherence between states',states[i-1],'and',states[j-1],'has equation 0 =',rhs
			print 'ERROR: the following equations',excluded_mu,' must be excluded from A.'
			print 'Try again giving write_stationary the argument excluded_mu=',excluded_mu,'.'
			print 'This argument must also be given to read_result'
			raise ValueError

		else:
			#So far zero_row has been equal to zero_col, but if this is
			#not the case we should look into it.
			print 'zero_row',zero_row
			print 'zero_col',zero_col
			print 'Were diffent, what does this mean pysically?'
			raise NotImplementedError

	return excluded_mu

def write_stationary(path,name,laser,omega,gamma,r,Lij,
				states=None,excluded_mu=[],verbose=1):
	t0=time()
	Ne=len(omega[0])
	Nl=len(laser)
	N_excluded_mu=len(excluded_mu)

	from config import use_netcdf
	
	code0="""program stationary_rho
    implicit none
    real*8, dimension("""+str(Nl)+""") :: E0,detuning_knob,detuning_knobi
    real*8, allocatable, dimension(:,:) :: rho
    real*8, allocatable, dimension(:) :: delta
    integer :: i,ldelta,ndelta,nerrors,info
    real*8 :: ddelta
    logical :: print_steps,save_systems,specific_deltas,use_netcdf
    real*4 :: start_time, end_time
    
    !We load the parameters\n"""
	# We break the path name into several lines if it is needed.
	long_line="""open(unit=2,file='"""+path+name+"""_params.dat',status='unknown')\n"""
	if len(long_line)>=72:
		long_line=long_line[:72]+"&\n&"+long_line[72:]
	code0+=long_line
	code0+="""
    read(2,*) E0
    read(2,*) detuning_knob
	read(2,*) ldelta
	read(2,*) ndelta
	read(2,*) ddelta
	read(2,*) print_steps
	read(2,*) save_systems
	read(2,*) use_netcdf
	read(2,*) specific_deltas

	!We read the deltas if they are given.
	allocate(delta(ndelta),stat=info)
	if (specific_deltas) then
		read(2,*) delta
	else
		!Or else build them.
		delta(1)=detuning_knob(ldelta)
		do i=2,ndelta
			delta(i)=detuning_knob(ldelta)+(i-1)*ddelta
		end do
	end if

    close(2)

	allocate(rho(ndelta,"""+str(Ne**2-1)+"""),stat=info)

	nerrors=0	

	call cpu_time(start_time)

	!$OMP PARALLEL PRIVATE(detuning_knobi)
	!$OMP DO
	do i=1,ndelta
		
		detuning_knobi=detuning_knob
		detuning_knobi(ldelta)=delta(i)
		
		call solve(E0,detuning_knobi,rho(i,:),save_systems)

		if (print_steps) print*,'delta=',detuning_knobi(ldelta)
		
	end do
	!$OMP END DO
	!$OMP END PARALLEL

	call cpu_time(end_time)
	if (print_steps) print*,'total time:',end_time-start_time

	!We write the result to a file.
	"""
	if use_netcdf:
		code0+="""call save_matrix_and_vector('"""+path+name+""".nc',ndelta,"""+str(Ne**2-1)+""",rho,delta)
"""
	else:
		code0+="""open(unit=1,file='"""+path+name+""".dat',status='unknown')
	do i=1,ndelta
		write(1,*) delta(i), rho(i,:)
	end do
	close(1)
	"""
	code0+="""

	deallocate(rho,stat=info)

end program
"""
	
	if use_netcdf:
		code0+="""
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
"""
	code0+="""
subroutine solve(E0,detuning_knob,B,save_systems)
	implicit none
	
	real*8, dimension("""+str(Nl)+"""), intent(in) :: E0,detuning_knob
	real*8, dimension("""+str(Ne**2-1-N_excluded_mu)+""",1), intent(out) :: B
	logical, intent(in) :: save_systems

"""
	####################################################################

	dummy=write_equations_code(path,name,laser,omega,gamma,r,Lij,
				states=states,excluded_mu=excluded_mu,verbose=verbose)
	code,Nd,row_check,col_check,rhs_check,Ne,N_excluded_mu,states,omega_min,detuningsij,omega_rescaled=dummy
	####################################################################
	code0+="""
	real*8, dimension("""+str(Nd)+""") :: detuning
	
	integer :: INFO,j
	real*8, dimension("""+str(Ne**2-1-N_excluded_mu)+""","""+str(Ne**2-1-N_excluded_mu)+""") :: A
	integer, dimension("""+str(Ne**2-1-N_excluded_mu)+""") :: IPIV

	A=0
	B=0
	"""
	#We make the program print A and B
#	code+="""	do j=1,"""+str(Ne**2-1)+"""
#		print*, A(j,:)
#	end do
#		print*,
#		do j=1,"""+str(Ne**2-1)+"""
#			print*,B(j,:)
#		end do
#		print* \n"""
	####################################################################
	#We make LAPACK solve the damn thing.
	code+="	if (save_systems) then\n"
	code+="		open(file='"+path+name+"_AB.dat',unit=4,status='unknown')\n"
	code+="		do j=1,"+str(Ne**2-1-N_excluded_mu)+'\n'
	code+="			write(4,*) A(j,:),B(j,1)\n"
	code+="		end do\n"
	code+="		close(4)\n"
	code+="	end if\n\n"
	
	code+='	call dgesv('+str(Ne**2-1-N_excluded_mu)+', 1, A, '+str(Ne**2-1-N_excluded_mu)+', IPIV, B, '+str(Ne**2-1-N_excluded_mu)+', INFO)\n'
	#code+="	print*,'INFO',INFO\n"
	code+="""	if (INFO>0) B=-11\n"""
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
		code=code0+code+'end subroutine\n'
		f.write(code)
		f.close()
		
		return time()-t0

def run_stationary(path,name,E0,laser_frequencies, spectrum_of_laser,N_delta,
				frequency_step=None,frequency_end=None,print_steps=False,
				specific_deltas=None, save_systems=False,clone=None,use_netcdf=True):
	t0=time()
	
	params=''.join([str(i)+' ' for i in E0])+'\n'
	params+=''.join([str(i)+' ' for i in laser_frequencies])+'\n'
	params+=str(spectrum_of_laser)+'\n'
	params+=str(N_delta)+'\n'
	
	if frequency_end !=None:
		if frequency_step !=None:
			raise ValueError,'both frequency_end and frequency_step were specified.'
		if N_delta==1:
			frequency_step=0.0
		else:
			frequency_step=(frequency_end-laser_frequencies[spectrum_of_laser-1])/(N_delta-1)	
	
	params+=str(frequency_step)+'\n'
	if print_steps:
		params+='.true.\n'
	else:
		params+='.false.\n'

	if save_systems:
		params+='.true.\n'
	else:
		params+='.false.\n'

	if use_netcdf:
		params+='.true.\n'
	else:
		params+='.false.\n'

	if specific_deltas==None:
		params+='.false.\n'
	else:
		params+='.true.\n'
		params+=''.join([str(delta)+' ' for delta in specific_deltas])
	
	if clone!=None:
		clone='_'+str(clone)
	else:
		clone=''
	
	f=file(path+name+'_params'+clone+'.dat','w')
	f.write(params)
	f.close()
	
	com=path+name+clone
	#print 'running',com
	exit_code=os.system(com)
	if exit_code != 0:
		s='command '+com+' returned exit code '+str(exit_code)
		raise RuntimeError,s
	return time()-t0

