!=======================================================================================================================
!> @file        global_initialize.f90
!> @brief       Module for initializing and validating the 3D example data.
!> @details     3D example solution is generated, and the corresponding initial data is calculated. 
!>					 The example solution will be used to validate whether the results of solving the TDM are correct.
!>					 Also, as noted in communication.f90, please be aware of the memory size limitations when generating the initial 3D example solution.
!>
!> @author      
!>              - Seung-chan Kim (macks1029@postech.ac.kr), Flow Physics and Engineering Lab., Pohang University of Science and Technology
!>              - Ji-hoo Kim (hugh577@postech.ac.kr), Flow Physics and Engineering Lab., Pohang University of Science and Technology
!>              - Ho-jun Moon (mhj2013@postech.ac.kr), Flow Physics and Engineering Lab., Pohang University of Science and Technology
!> @date        Aug 2024
!> @version     2.0
!> @par         Copyright
!>              Copyright (c) 2020-Present Seugn-chan Kim, Ji-hoo Kim, Ho-jun Moon and Dong-hyun you, Pohang University of Science and Technology
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE in )
!=======================================================================================================================

!>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!
!!    GLOBAL VARIABLES/SUBROUTINES FOR THE HOST (CPU)
!!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module m_global_initialize

contains
!>=====================================================================
subroutine read_infile
    
	use pdd_host

	implicit none
	integer :: fh 
	
	pi = 4.0_cgreal*atan(1.0_cgreal)
	
	!>-- open files
	fh = 50
	open(fh,file='input.dat',status='old')
	read(fh,*)
	read(fh,*) nproc_x, nproc_y , nproc_z
	read(fh,*) 
	read(fh,*) ndim, dir
	read(fh,*)
	read(fh,*) nx0,ny0,nz0

   !>-- Print input
   if (rank == 0) then
		print *,'' 
		print *,&
		'------------------------------------------------------------------'
		write(*,'(a,i10)') ' # of Processors in x-dir : ', nproc_x
		write(*,'(a,i10)') ' # of Processors in y-dir : ', nproc_y
		write(*,'(a,i10)') ' # of Processors in x-dir : ', nproc_z
		print *,''
		write(*,'(a,i16)') ' # of cells in x-dir :', nx0
		write(*,'(a,i16)') ' # of cells in y-dir :', ny0
		write(*,'(a,i16)') ' # of cells in x-dir :', nz0
	endif

end subroutine read_infile
    
!>=====================================================================
subroutine  init_3Ddata
    
	use pdd_host
   use global_host
   use m_pdd_initialize
  	use m_tdma_main
	use m_communication

	implicit none
   include 'mpif.h'
   real :: i, j, k, tmp1, tmp2, tmp3
	
	!>======================================================================
	!!    Generate random data "x" in system of equation Ax=b 
	!!======================================================================
	if (rank == 0) then
		do k = 1,nz0
			do j = 1,ny0
				do i = 1,nx0
					call random_number(tmp1)
					call random_number(tmp2)
					call random_number(tmp3)
					dat_orig(i,j,k) = tmp1
					glb_mat_a(i,j,k) = 0.1 * tmp2
					glb_mat_b(i,j,k) = 1
					glb_mat_c(i,j,k) = 0.1 * tmp3
				enddo
			enddo
		enddo
	endif

	!>======================================================================
	!!    Calculate data "b" in system of equation Ax=b
	!!    Make data in rank =0 and distribute to all processors.
	!!======================================================================
	if (rank == 0) then
		if (dir == 1) then
			do k = 1,nz0
				do j = 1,ny0
					do i = 2,nx0-1
						glb_dat(i,j,k) =   glb_mat_a(i,j,k)*dat_orig(i-1,j,k) &
											  + glb_mat_b(i,j,k)*dat_orig(i  ,j,k) &
											  + glb_mat_c(i,j,k)*dat_orig(i+1,j,k) 
					enddo
				enddo
			enddo
			do k = 1,nz0
				do j = 1,ny0
					glb_dat(1,j,k) =    glb_mat_b(1,j,k)*dat_orig(1,j,k) &
										   + glb_mat_c(1,j,k)*dat_orig(2,j,k) 
					glb_dat(nx0,j,k) =    glb_mat_a(nx0,j,k)*dat_orig(nx0-1,j,k) &
										     + glb_mat_b(nx0,j,k)*dat_orig(nx0  ,j,k) 
				enddo
			enddo
		endif

		if (dir == 2) then
			do k = 1,nz0
				do j = 2,ny0-1
					do i = 1,nx0
						glb_dat(i,j,k) =   glb_mat_a(i,j,k)*dat_orig(i,j-1,k) &
											  + glb_mat_b(i,j,k)*dat_orig(i,j  ,k) &
											  + glb_mat_c(i,j,k)*dat_orig(i,j+1,k) 
					enddo
				enddo
			enddo
			do k = 1,nz0
				do i = 1,nx0
					glb_dat(i,1,k) =    glb_mat_b(i,1,k)*dat_orig(i,1,k) &
										   + glb_mat_c(i,1,k)*dat_orig(i,2,k) 
					glb_dat(i,ny0,k) =    glb_mat_a(i,ny0,k)*dat_orig(i,ny0-1,k) &
										     + glb_mat_b(i,ny0,k)*dat_orig(i,ny0  ,k) 
				enddo
			enddo
		endif

		if (dir == 3) then
			do k = 2,nz0-1
				do j = 1,ny0
					do i = 1,nx0
						glb_dat(i,j,k) =   glb_mat_a(i,j,k)*dat_orig(i,j,k-1) &
											  + glb_mat_b(i,j,k)*dat_orig(i,j,k  ) &
											  + glb_mat_c(i,j,k)*dat_orig(i,j,k+1) 
					enddo
				enddo
			enddo
			do j = 1,ny0
				do i = 1,nx0
					glb_dat(i,j,1) =    glb_mat_b(i,j,1)*dat_orig(i,j,1) &
										   + glb_mat_c(i,j,1)*dat_orig(i,j,2) 
					glb_dat(i,j,nz0) =    glb_mat_a(i,j,nz0)*dat_orig(i,j,nz0-1) &
										     + glb_mat_b(i,j,nz0)*dat_orig(i,j,nz0  ) 
				enddo
			enddo
		endif
	endif
	
	!>-- Data distribution
	call scatter_main  !<-- src subroutine


	call MPI_Barrier(MPI_COMM_WORLD,ierr)

end subroutine init_3Ddata
    
!>=====================================================================
subroutine  validation(maxerr)
    
   use global_host
   use pdd_host
   use m_pdd_initialize
  	use m_tdma_main
	use m_communication

	implicit none
   include 'mpif.h'
   real, intent(out) :: maxerr 
	!---------------------------
   real :: i, j, k
	

	!>-- Data gathering
	call gather_main  !<-- src subroutine

	!>======================================================================
	!!    Generate random data "x" in system of equation Ax=b 
	!!======================================================================
	if (rank == 0) then
		maxerr = 0
		do k = 1,n3f
			do j = 1,n2f
				do i = 1,n1f
					if ( maxerr < abs(glb_dat(i,j,k) - dat_orig(i,j,k)) ) then
						maxerr = abs( glb_dat(i,j,k) - dat_orig(i,j,k) )
					endif
				enddo
			enddo
		enddo
	endif

	call MPI_Barrier(MPI_COMM_WORLD,ierr)

end subroutine validation
    
!>======================================================================
end module m_global_initialize

