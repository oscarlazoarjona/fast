program stationary_rho
    implicit none
    real*8, dimension(1) :: E0,detuning_knob,detuning_knobi
    real*8, allocatable, dimension(:,:) :: rho
    real*8, allocatable, dimension(:) :: delta
    integer :: i,ldelta,ndelta,nerrors,info
    real*8 :: ddelta
    logical :: print_steps,save_systems,specific_deltas,use_netcdf
    real*4 :: start_time, end_time

    !We load the parameters
open(unit=2,file='folder_01___Two_level_atom/suite_steady_params.dat',st&
&atus='unknown')

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

    allocate(rho(ndelta,3),stat=info)

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
    open(unit=1,file='folder_01___Two_level_atom/suite_steady.dat',status='unknown')
    do i=1,ndelta
        write(1,*) delta(i), rho(i,:)
    end do
    close(1)
    

    deallocate(rho,stat=info)

end program

subroutine solve(E0,detuning_knob,B,save_systems)
    implicit none

	real*8, dimension(1), intent(in) :: E0,detuning_knob
	real*8, dimension(3,1), intent(out) :: B
    logical, intent(in) :: save_systems


	real*8, dimension(1) :: detuning

    integer :: INFO,j
	real*8, dimension(3,3) :: A
	integer, dimension(3) :: IPIV

    A=0
    B=0
    	!We calculate the detunings.
	!The list of detunings has the following meaning:
	!detuning(1)= delta^1_,
	detuning(1)=detuning_knob(1) -(0.0d0)-(0.0d0)

	!We calculate the independent vector.
	B(3,1)=B(3,1) +E0(1)*(1.0d0)

	B=B/2.0d0

	!We calculate the equations for populations.
	A(1,3)=A(1,3)+E0(1)*(-2.0d0)

	!The code to calculate the equations for coherences.
	A(3,1)=A(3,1)+E0(1)*(1.0d0)
	A(3,1)=A(3,1)+E0(1)*(1.0d0)

	A=A/2.0d0

	!We calculate the terms associated with the phase transformation.
    A(2,3)=A(2,3)-(detuning(1))
    A(3,2)=A(3,2)+(detuning(1))

	!We calculate the terms associated with spontaneous decay.
    A(1,1)=A(1,1)-(1.0d0)
    A(2,2)=A(2,2)-(0.5d0)
    A(3,3)=A(3,3)-(0.5d0)

	if (save_systems) then
		open(file='folder_01___Two_level_atom/suite_steady_AB.dat',unit=4,status='unknown')
		do j=1,3
			write(4,*) A(j,:),B(j,1)
		end do
		close(4)
	end if

	call dgesv(3, 1, A, 3, IPIV, B, 3, INFO)
	if (INFO>0) B=-11
end subroutine
