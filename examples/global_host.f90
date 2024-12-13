!=======================================================================================================================
!> @file        global_host.f90
!> @brief       Module that contains global CPU variables used in the example folder.
!> @details     all global variables related to example case are defined in here.
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

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!
!!    GLOBAL VARIABLES/SUBROUTINES FOR THE HOST (CPU)
!!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>======================================================================
module global_host

use pdd_host
use m_precision
implicit none

!>--Save original data
real(kind=cgreal), dimension(:,:,:), allocatable :: dat_orig

contains
!>======================================================================
!!    CPU global memory allocation
!!======================================================================
subroutine malloc_global

	implicit none
   include 'mpif.h'

	!>--Save original data
	if (rank == 0) allocate( dat_orig( nx,ny,nz ) )

	!>--global 3D data
	if(rank == 0) allocate( glb_dat(1:nx,1:ny,1:nz) )
	if(rank == 0) allocate( glb_mat_a(1:nx,1:ny,1:nz) )
	if(rank == 0) allocate( glb_mat_b(1:nx,1:ny,1:nz) )
	if(rank == 0) allocate( glb_mat_c(1:nx,1:ny,1:nz) )

	!>--Communication file temporary parameter
	allocate( var_dump( n1i:n1f,n2i:n2f,n3i:n3f ) )
	if (rank == 0) allocate( glb_dat_dump( nx,ny,nz ) )

end subroutine malloc_global

!>======================================================================
!!    Free CPU global memory
!!======================================================================
subroutine free_global

	implicit none

	!>--Save original data
	if (rank == 0) deallocate( dat_orig )
	
	!>--global 3D data
	if(rank == 0) deallocate( glb_dat )
	if(rank == 0) deallocate( glb_mat_a )
	if(rank == 0) deallocate( glb_mat_b )
	if(rank == 0) deallocate( glb_mat_c )
	
	!>--Communication file temporary parameter
	deallocate( var_dump )
	if (rank == 0) deallocate( glb_dat_dump )

end subroutine free_global

!>======================================================================
end module global_host

