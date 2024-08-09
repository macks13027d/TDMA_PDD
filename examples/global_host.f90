!======================================================================================================================
!> @brief       PDD_TDMA - Parallel and Scalable Library for TriDiagonal Matrix Algorithm
!> @details     PaScaL_TDMA includes a CUDA implementation of PaScaL_TDMA, which accelerates 
!>              to solve many tridiagonal systems in multi-dimensional partial differential equations on GPU.
!>              It adopts the pipeline copy within the shared memory for the forward elemination and 
!>              backward substitution procudures of TDMA to reduce global memory access.
!>              For the main algorithm of PaScaL_TDMA, see also https://github.com/MPMC-Lab/PaScaL_TDMA.
!> 
!> @author      
!>              - Mingyu Yang (yang926@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Ki-Ha Kim (k-kiha@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>              - Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>
!> @date        May 2023
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

end subroutine malloc_global

!>======================================================================
!!    Free CPU global memory
!!======================================================================
subroutine free_global

	implicit none

	!>--Save original data
	if (rank == 0) deallocate( dat_orig )

end subroutine free_global

!>======================================================================
end module global_host

