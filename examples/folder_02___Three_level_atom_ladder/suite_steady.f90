program stationary_rho
    implicit none
    real*8, dimension(2) :: E0,detuning_knob,detuning_knobi
    real*8, allocatable, dimension(:,:) :: rho
    real*8, allocatable, dimension(:) :: delta
    integer :: i,ldelta,ndelta,nerrors,info
    real*8 :: ddelta
    logical :: print_steps,save_systems,specific_deltas,use_netcdf
    real*4 :: start_time, end_time

    !We load the parameters
open(unit=2,file='folder_02___Three_level_atom_ladder/suite_steady_param&
&s.dat',status='unknown')

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

    allocate(rho(ndelta,8),stat=info)

    nerrors=0

    call cpu_time(start_time)

    !$OMP PARALLEL PRIVATE(detuning_knobi)
    !$OMP DO
    do i=1, ndelta

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
    open(unit=1,file='folder_02___Three_level_atom_ladder/suite_steady.dat',status='unknown')
    do i=1,ndelta
        write(1,*) delta(i), rho(i,:)
    end do
    close(1)
    

    deallocate(rho,stat=info)

end program

subroutine solve(E0,detuning_knob,B,save_systems)
    implicit none

	real*8, dimension(2), intent(in) :: E0,detuning_knob
	real*8, dimension(8,1), intent(out) :: B
    logical, intent(in) :: save_systems


	real*8, dimension(2) :: detuning

    integer :: INFO,j
	real*8, dimension(8,8) :: A
	integer, dimension(8) :: IPIV

    A=0
    B=0
    	!We calculate the detunings.
	!The list of detunings has the following meaning:
	!detuning(1)= delta^1_,
	!detuning(2)= delta^2_,
	detuning(1)=detuning_knob(1) -(0d0)-(0d0)
	detuning(2)=detuning_knob(2) -(0d0)-(0d0)

	!We calculate the independent vector.
	B(6,1)=B(6,1) +E0(1)*(1.0d0)

	B=B/2.0d0

	!We calculate the equations for populations.
	A(1,6)=A(1,6)+E0(1)*(-2.0d0)
	A(1,8)=A(1,8)+E0(2)*(2.0d0)
	A(2,8)=A(2,8)+E0(2)*(-2.0d0)

	!The code to calculate the equations for coherences.
	A(3,7)=A(3,7)+E0(2)*(1.0d0)
	A(6,4)=A(6,4)+E0(2)*(-1.0d0)
	A(6,1)=A(6,1)+E0(1)*(1.0d0)
	A(6,1)=A(6,1)+E0(1)*(1.0d0)
	A(6,2)=A(6,2)+E0(1)*(1.0d0)
	A(4,8)=A(4,8)+E0(1)*(-1.0d0)
	A(4,6)=A(4,6)+E0(2)*(1.0d0)
	A(7,5)=A(7,5)+E0(1)*(1.0d0)
	A(7,3)=A(7,3)+E0(2)*(-1.0d0)
	A(5,7)=A(5,7)+E0(1)*(-1.0d0)
	A(8,4)=A(8,4)+E0(1)*(1.0d0)
	A(8,2)=A(8,2)+E0(2)*(1.0d0)
	A(8,1)=A(8,1)+E0(2)*(-1.0d0)

	A=A/2.0d0

	!We calculate the terms associated with the phase transformation.
    A(3,6)=A(3,6)-(detuning(1))
    A(6,3)=A(6,3)+(detuning(1))
    A(4,7)=A(4,7)-(detuning(1)+detuning(2))
    A(7,4)=A(7,4)+(detuning(1)+detuning(2))
    A(5,8)=A(5,8)-(detuning(2))
    A(8,5)=A(8,5)+(detuning(2))

	!We calculate the terms associated with spontaneous decay.
    A(1,1)=A(1,1)-(6.0d0)
    A(1,2)=A(1,2)-(-0.6d0)
    A(2,2)=A(2,2)-(0.6d0)
    A(5,5)=A(5,5)-(3.0d0)
    A(8,8)=A(8,8)-(3.0d0)
    A(3,3)=A(3,3)-(3.0d0)
    A(6,6)=A(6,6)-(3.0d0)
    A(4,4)=A(4,4)-(0.3d0)
    A(7,7)=A(7,7)-(0.3d0)
    A(5,5)=A(5,5)-(0.3d0)
    A(8,8)=A(8,8)-(0.3d0)

	if (save_systems) then
		open(file='folder_02___Three_level_atom_ladder/suite_steady_AB.dat',unit=4,status='unknown')
		do j=1,8
			write(4,*) A(j,:),B(j,1)
		end do
		close(4)
	end if

	call dgesv(8, 1, A, 8, IPIV, B, 8, INFO)
	if (INFO>0) B=-11
end subroutine
