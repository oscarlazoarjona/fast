program evolution_diagonalization
	implicit none
    real*8 :: dt,ddelta,delta0
	real*8, allocatable, dimension(:) :: t,delta
	integer :: i,j,mu,n,ldelta,ndelta,detuning_index,n_aprox,n_mod,info

	logical :: print_steps,run_spectrum,save_systems,save_eigenvalues,integrate,use_netcdf
    real*8, dimension(2) :: E0,detuning_knob
    real*8, dimension(2) :: detuning

	complex*16, dimension(8) :: r_amp,rho0
	complex*16, allocatable, dimension(:,:) :: rho,rho_spectrum
	real*8, dimension(8) :: rho_inf
	complex*16, dimension(8,8) :: U
	complex*16, dimension(8) :: lam

    !We load the parameters
    n_aprox=1500
    open(unit=2,file='folder_02___Three_level_atom_ladder/suite_evolutio&
&n_params.dat',status='unknown')
    read(2,*) n
    read(2,*) dt
    read(2,*) print_steps
    read(2,*) rho0
    read(2,*) E0
    read(2,*) detuning_knob
    read(2,*) run_spectrum
    read(2,*) save_systems
    read(2,*) save_eigenvalues
    read(2,*) use_netcdf
    read(2,*) integrate

    if (run_spectrum) then
		read(2,*) ldelta
		read(2,*) ndelta
		read(2,*) ddelta


		allocate(rho(n,8),stat=info)
		allocate(rho_spectrum(ndelta,8),stat=info)
		rho(1,:)=rho0
    else
		ldelta=1; ndelta=1; ddelta=0

		allocate(rho(n,8),stat=info)
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

    if (n_mod==0) n_mod=1

    if (save_eigenvalues) then
open(unit=3,file='folder_02___Three_level_atom_ladder/suite_evolution_ei&
&genvalues.dat',status='unknown')
	end if

	!We start the detuning variation, which might consist of a single detuning.
    do j=1,ndelta
		!We calculate the solution for specific E0, detuning_knob, and initial rho,
		!this returns U, r_amp, rho_inf, lam.
		detuning_knob(ldelta)=delta(j)

		call solve(E0,detuning_knob, rho0,U,r_amp,rho_inf,lam,save_systems)

		if (save_eigenvalues) then
			write(3,*) delta,real(lam)
			write(3,*) delta,imag(lam)
		end if


		!We calculate the time evolution with the solution just computed.
		do i=2,n

			if (print_steps.and. .not. run_spectrum) print*,'t=',t,'delta=',delta

			do mu=1,8
				rho(i,mu)=r_amp(mu)*cdexp(lam(mu)*t(i))
			end do
			rho(i,:)= matmul(U , rho(i,:)) + rho_inf

			if (isnan(real(rho(i,1)))) stop 1
			!if (.not. run_spectrum ) WRITE(1,*) t,real(rho(i,:))
		end do

		if (print_steps) print*, 'delta=',delta,'percentage=',100*(delta-delta0)/(ddelta*ndelta)

		!We save the time evolution.
		if (.not. run_spectrum ) then
			open(unit=1,file='folder_02___Three_level_atom_ladder/suite_evolution&
&.dat',status='unknown')

			do i=1,n
				WRITE(1,*) t(i),real(rho(i,:))
			end do
			close(1)
		else
			rho_spectrum(j,:)=rho(n,:)
		end if
		
	end do

	!We save a spectrum.
	if (run_spectrum) then

		open(unit=1,file="folder_02___Three_level_atom_ladder/suite_evolution.dat",status='unknown')
		do i=1,ndelta
			WRITE(1,*) delta(i),real(rho_spectrum(i,:))
		end do
		close(1)
	end if

end program
subroutine solve(E0,detuning_knob,rho0,U,r_amp,rho_inf,lam,save_systems)
	implicit none
    integer :: i,j,k,mu,nu,alpha
	real*8, dimension(2), intent(in) :: E0,detuning_knob
	complex*16, dimension(8,1), intent(in) :: rho0
	complex*16, dimension(8,8), intent(out) :: U
	complex*16, dimension(8,1), intent(out) :: r_amp
	real*8, dimension(8,1), intent(out) :: rho_inf
	complex*16, dimension(8), intent(out) :: lam
	logical, intent(in) :: save_systems
	real*8, dimension(2) :: detuning
	real*8, dimension(8,8) :: A
	complex*16, dimension(8,1) :: B
	integer, dimension(8) :: IPIV

    integer :: Nrho
    real*8, dimension(8) :: WR,WI
    real*8, dimension(1,8) :: VL
    real*8, dimension(8,8) :: VR,VRT
    real*8, dimension(48) :: WORK
    integer :: INFO

    complex*16 :: II
    complex*16, dimension(8,8) :: Ui,Us,Uis
    logical :: band

	Nrho=8
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
    open(unit=4,file='folder_02___Three_level_atom_ladder/suite_evolutio&
&n_AB.dat',status='unknown')
		do j=1,8
			write(4,*) A(j,:),B(j,1)
		end do
		close(4)
	end if


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!We begin the solution process.
    call DGEEV( 'N', 'V', Nrho, A, Nrho, WR, WI, VL, 1, VR, Nrho, WORK, 6*Nrho, INFO )
    II=(0,1)
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
    end do

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
	Ui=Us
	!We calculate d (which will be made by successive modifications of b)
	b=matmul(Ui,b)
	do mu=1,Nrho
		b(mu,1)=b(mu,1)/lam(mu)
	end do
	!We calculate the vector r
	r_amp=matmul(Ui,rho0) - b
	!We calculate rho_inf
	rho_inf=matmul(U,b)
end subroutine
