!>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!
!!    PDD RELATED VARIABLES/SUBROUTINES FOR THE HOST (CPU)
!!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module m_precision

! Real kinds
integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real
integer, parameter :: cgreal = selected_real_kind(15,300)     ! double precision real

! Integer kinds
integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

!Complex kinds
integer, parameter :: kc4 = kr4                            ! single precision complex
integer, parameter :: kc8 = kr8                            ! double precision complex

end module m_precision

!>======================================================================
module pdd_host

use m_precision
implicit none

!>--Constant
real(kind=cgreal) :: pi

!>--Input parameter
integer  :: nproc_x, nproc_y, nproc_z
integer  :: ndim, dir
integer  :: nx0,ny0,nz0

!>--MPI variables
integer                               :: rank, nprocs

!>--MPI coordinate description
integer                            :: comm_cart, ierr
integer, dimension(:), allocatable :: pdim, coords
logical, dimension(:), allocatable :: period
integer  						   	  :: n1i, n1f, n2i, n2f, n3i, n3f

!>--MPI communication parameter 1
integer, dimension(:), allocatable :: is_neig,ie_neig,js_neig,je_neig,ks_neig,ke_neig
integer  :: nxloc_mpi, nyloc_mpi, nzloc_mpi, nloc_mpi, nx, ny, nz
integer, dimension(:), allocatable :: thrDR_neig

!>--MPI communication parameter 2 (Data collection)
integer :: comm_new_x,comm_new_y,comm_new_z
integer :: rank_new_x,rank_new_y,rank_new_z
integer :: newtype1,newtype2,newtype3

!>--3D data
real(kind=cgreal), dimension(:,:,:), allocatable :: dat0, dat1_j, dat1_k, glb_dat
real(kind=cgreal), dimension(:,:,:), allocatable :: mat_a, mat_b, mat_c 
real(kind=cgreal), dimension(:,:,:), allocatable :: glb_mat_a, glb_mat_b, glb_mat_c 

!>--TDMA_3d temporary parameter
real(kind=cgreal), dimension(:),     allocatable :: sol_tmp_i, xx_i, vv_i, ww_i
real(kind=cgreal), dimension(:,:),   allocatable :: xxs_i, vvs_i, xxr_i, vvr_i, xxe_i, wwe_i, yv_i, yw_i, yv2_i
real(kind=cgreal), dimension(:,:,:), allocatable :: vv_arr_i, ww_arr_i

real(kind=cgreal), dimension(:),     allocatable :: sol_tmp_j, xx_j, vv_j, ww_j
real(kind=cgreal), dimension(:,:),   allocatable :: xxs_j, vvs_j, xxr_j, vvr_j, xxe_j, wwe_j, yv_j, yw_j, yv2_j
real(kind=cgreal), dimension(:,:,:), allocatable :: vv_arr_j, ww_arr_j

real(kind=cgreal), dimension(:),     allocatable :: sol_tmp_k, xx_k, vv_k, ww_k
real(kind=cgreal), dimension(:,:),   allocatable :: xxs_k, vvs_k, xxr_k, vvr_k, xxe_k, wwe_k, yv_k, yw_k, yv2_k
real(kind=cgreal), dimension(:,:,:), allocatable :: vv_arr_k, ww_arr_k

!>--PDD algorithm temporary parameter
real(kind=cgreal), dimension(:),     allocatable :: gam

!>--Communication file temporary parameter
real(kind=cgreal), dimension(:,:,:),     allocatable :: var_dump, glb_dat_dump 


contains
!>======================================================================
!!    CPU PDD-related memory allocation
!!======================================================================
subroutine malloc_pdd

	implicit none
   include 'mpif.h'

	!>--MPI coordinate description
   allocate( pdim(ndim), period(ndim), coords(ndim) )

	!>--MPI communication parameter 1
	allocate( is_neig(0:nprocs-1), ie_neig(0:nprocs-1) )
	allocate( js_neig(0:nprocs-1), je_neig(0:nprocs-1) )
	allocate( ks_neig(0:nprocs-1), ke_neig(0:nprocs-1) )
	if(rank.eq.0) allocate(thrdr_neig(1:nprocs-1))

	!>===================================================================
	!! 						MPI coordinate description
	!!===================================================================
	!>-- create coordinate variable
   if (ndim == 3) pdim   = (/ nproc_x, nproc_y, nproc_z /)
   if (ndim == 3) period = (/ .false., .false., .false. /)
   
   call MPI_CART_CREATE( MPI_COMM_WORLD, ndim, pdim, period, .true., comm_cart, ierr )
   call MPI_CART_GET( comm_cart, ndim, pdim, period, coords, ierr )
       
	!>-- MPI x coordinate description
!	if(coords(1).lt.mod((nx),nproc_x))then
!	  nxloc_mpi   = int((nx)/nproc_x)+1
!	  n1i      = coords(1)*nxloc_mpi+1
!	  n1f      = n1i + nxloc_mpi-1
!	else
!	  nxloc_mpi   = int((nx)/nproc_x)
!	  n1i      =   mod((nx),nproc_x)*(nxloc_mpi+1)&
!	                +(coords(1)-mod((nx),nproc_x))*nxloc_mpi+1
!	  n1f      = n1i + nxloc_mpi-1
!	endif
   

   nxloc_mpi   =  ceiling(real(nx0)/nproc_x)
   n1i      =  coords(1)*nxloc_mpi+1
   n1f      =  n1i + nxloc_mpi-1
	

!>-- MPI y coordinate description
!	if(coords(2).lt.mod((ny),nproc_y))then
!	  nyloc_mpi   = int((ny)/nproc_y)+1
!	  n2i      = coords(2)*nyloc_mpi+1
!	  n2f      = n2i + nyloc_mpi-1
!	else
!	  nyloc_mpi   = int((ny)/nproc_y)
!	  n2i      =   mod((ny),nproc_y)*(nyloc_mpi+1)&
!	                +(coords(2)-mod((ny),nproc_y))*nyloc_mpi+1
!	  n2f      = n2i + nyloc_mpi-1
!	endif

   nyloc_mpi   =  ceiling(real(ny0)/nproc_y)
   n2i      =  coords(2)*nyloc_mpi+1
   n2f      =  n2i + nyloc_mpi-1
	
	!>-- MPI z coordinate description
!	if(coords(3).lt.mod((nz),nproc_z))then
!	  nzloc_mpi   = int((nz)/nproc_z)+1
!	  n3i      = coords(3)*nzloc_mpi+1
!	  n3f      = n3i + nzloc_mpi-1
!	else
!	  nzloc_mpi   = int((nz)/nproc_z)
!	  n3i      =   mod((nz),nproc_z)*(nzloc_mpi+1)&
!	                +(coords(3)-mod((nz),nproc_z))*nzloc_mpi+1
!	  n3f      = n3i + nzloc_mpi-1
!	endif

   nzloc_mpi   =  ceiling(real(nz0)/nproc_z)
   n3i      =  coords(3)*nzloc_mpi+1
   n3f      =  n3i + nzloc_mpi-1

	nx = nxloc_mpi * nproc_x
	ny = nyloc_mpi * nproc_y
	nz = nzloc_mpi * nproc_z
	
	!>--3D data
	allocate( mat_a(n1i:n1f,n2i:n2f,n3i:n3f) )
	allocate( mat_b(n1i:n1f,n2i:n2f,n3i:n3f) )
	allocate( mat_c(n1i:n1f,n2i:n2f,n3i:n3f) )
	allocate( dat0(n1i:n1f,n2i:n2f,n3i:n3f) )

	!>--TDMA_3d temporary parameter
	!> dir = 1
	allocate( sol_tmp_i( n1i:n1f ), xx_i( n1i:n1f ), vv_i( n1i:n1f ), ww_i( n1i:n1f ) )
	allocate( xxs_i( n2i:n2f,n3i:n3f ), vvs_i( n2i:n2f,n3i:n3f ), xxr_i( n2i:n2f,n3i:n3f ) )
	allocate( vvr_i( n2i:n2f,n3i:n3f ), xxe_i( n2i:n2f,n3i:n3f ), wwe_i( n2i:n2f,n3i:n3f ) )
	allocate(  yv_i( n2i:n2f,n3i:n3f ),  yw_i( n2i:n2f,n3i:n3f ), yv2_i( n2i:n2f,n3i:n3f ) )
	allocate( vv_arr_i( n1i:n1f,n2i:n2f,n3i:n3f ), ww_arr_i( n1i:n1f,n2i:n2f,n3i:n3f ) )

	!> dir = 2
	allocate( sol_tmp_j( n2i:n2f ), xx_j( n2i:n2f ), vv_j( n2i:n2f ), ww_j( n2i:n2f ) )
	allocate( xxs_j( n1i:n1f,n3i:n3f ), vvs_j( n1i:n1f,n3i:n3f ), xxr_j( n1i:n1f,n3i:n3f ) )
	allocate( vvr_j( n1i:n1f,n3i:n3f ), xxe_j( n1i:n1f,n3i:n3f ), wwe_j( n1i:n1f,n3i:n3f ) )
	allocate(  yv_j( n1i:n1f,n3i:n3f ),  yw_j( n1i:n1f,n3i:n3f ), yv2_j( n1i:n1f,n3i:n3f ) )
	allocate( vv_arr_j( n2i:n2f,n1i:n1f,n3i:n3f ), ww_arr_j( n2i:n2f,n1i:n1f,n3i:n3f ) )

	allocate(dat1_j(n2i:n2f,n1i:n1f,n3i:n3f))

	!> dir = 3
	allocate( sol_tmp_k( n3i:n3f ), xx_k( n3i:n3f ), vv_k( n3i:n3f ), ww_k( n3i:n3f ) )
	allocate( xxs_k( n1i:n1f,n2i:n2f ), vvs_k( n1i:n1f,n2i:n2f ), xxr_k( n1i:n1f,n2i:n2f ) )
	allocate( vvr_k( n1i:n1f,n2i:n2f ), xxe_k( n1i:n1f,n2i:n2f ), wwe_k( n1i:n1f,n2i:n2f ) )
	allocate(  yv_k( n1i:n1f,n2i:n2f ),  yw_k( n1i:n1f,n2i:n2f ), yv2_k( n1i:n1f,n2i:n2f ) )
	allocate( vv_arr_k( n3i:n3f,n1i:n1f,n2i:n2f ), ww_arr_k( n3i:n3f,n1i:n1f,n2i:n2f ) )

	allocate(dat1_k(n3i:n3f,n1i:n1f,n2i:n2f))


	!>--PDD algorithm temporary parameter
	allocate( gam( -1:max(nx,ny,nz)+2 ) )
	
	!>--Communication file temporary parameter
	allocate( var_dump( n1i:n1f,n2i:n2f,n3i:n3f ) )
	if (rank == 0) allocate( glb_dat_dump( nx,ny,nz ) )
	
	!>--global 3D data
	if(rank == 0) allocate( glb_dat(1:nx,1:ny,1:nz) )
	if(rank == 0) allocate( glb_mat_a(1:nx,1:ny,1:nz) )
	if(rank == 0) allocate( glb_mat_b(1:nx,1:ny,1:nz) )
	if(rank == 0) allocate( glb_mat_c(1:nx,1:ny,1:nz) )

end subroutine malloc_pdd

!>======================================================================
!!    Free CPU PDD-related host memory
!!======================================================================
subroutine free_pdd

	implicit none

	!>--MPI coordinate description
   deallocate( pdim, period, coords )

	!>--MPI communication parameter 1
	deallocate( is_neig, ie_neig )
	deallocate( js_neig, je_neig )
	deallocate( ks_neig, ke_neig )

	!>--3D data
	deallocate( mat_a, mat_b, mat_c )
	deallocate( dat0 )

	!>--TDMA_3d temporary parameter
	deallocate( sol_tmp_i, xx_i, vv_i, ww_i )
	deallocate( xxs_i, vvs_i, xxr_i, vvr_i, xxe_i, wwe_i, yv_i, yw_i, yv2_i ) 
	deallocate( vv_arr_i, ww_arr_i ) 

	deallocate( sol_tmp_j, xx_j, vv_j, ww_j )
	deallocate( xxs_j, vvs_j, xxr_j, vvr_j, xxe_j, wwe_j, yv_j, yw_j, yv2_j ) 
	deallocate( vv_arr_j, ww_arr_j ) 

	deallocate( sol_tmp_k, xx_k, vv_k, ww_k )
	deallocate( xxs_k, vvs_k, xxr_k, vvr_k, xxe_k, wwe_k, yv_k, yw_k, yv2_k ) 
	deallocate( vv_arr_k, ww_arr_k ) 

	deallocate( dat1_j )
	deallocate( dat1_k )

	!>--PDD algorithm temporary parameter
	deallocate( gam )
	
	!>--Communication file temporary parameter
	deallocate( var_dump )
	if (rank == 0) deallocate( glb_dat_dump )
	
	!>--global 3D data
	if(rank == 0) deallocate( glb_dat )
	if(rank == 0) deallocate( glb_mat_a )
	if(rank == 0) deallocate( glb_mat_b )
	if(rank == 0) deallocate( glb_mat_c )

end subroutine free_pdd

!>======================================================================
end module pdd_host






   MODULE grid_arrays

    USE pdd_host
use m_precision
    
    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: x,y,z,xc,yc,zc
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dx,dy,dz
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxinv,dyinv,dzinv
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxc,dyc,dzc
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxcinv,dycinv,dzcinv
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: fx,fy,fz
    INTEGER :: order
    INTEGER,DIMENSION(26) :: dir_flag,add_flag
   
   END MODULE grid_arrays
!------------------------------------------------------

 

