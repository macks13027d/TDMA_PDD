!--------------------------------------------------------
!  SOLVE THE TRIDIAGONAL MATRIX 
!    [b1  c1           ]
!    [a2  b2  c2       ]
!    [...              ]
!    [...              ]
!    [           aN  bN]
! Reference: Numerical Recipes Fortran 77
!--------------------------------------------------------
module m_tdma_pdd

use pdd_host
contains
!>======================================================================
subroutine tdma(a,b,c,r,u,n1,n2)

   implicit none

   integer, intent(in)  :: n1, n2
   real(kind=cgreal), dimension(n1:n2), intent(in)  :: a, b, c
   real(kind=cgreal), dimension(n1:n2), intent(in)  :: r
   real(kind=cgreal), dimension(n1:n2), intent(out) :: u
	!-------------------------------------
   real(kind=cgreal)                       :: binv
   integer                                 :: j

   binv  = 1.0_cgreal/b(n1)
   u(n1) = r(n1)*binv

   do j = n1+1,n2
   	gam(j)  = c(j-1)*binv
   	binv    = 1.0_cgreal/(b(j)-a(j)*gam(j))
   	u(j)    = (r(j)-a(j)*u(j-1))*binv
   enddo

   do j = n2-1,n1,-1
   	u(j) = u(j) - gam(j+1)*u(j+1)
   enddo

   return
end subroutine tdma

!>======================================================================
subroutine tdma_pdd_x(a,b,c,sxx,svv,sww)

   implicit none
   include 'mpif.h'

   real(kind=cgreal), dimension(n1i:n1f), intent(in) :: a, b, c
   real(kind=cgreal), dimension(n1i:n1f), intent(inout) :: sxx
   real(kind=cgreal), dimension(n1i:n1f), intent(out) :: svv, sww
	!-------------------------------------
   real(kind=cgreal)                               :: binv
   integer :: i

	!>-------------------------------------------------------------------
	!!   Forward elimination
	!!-------------------------------------------------------------------
   binv  = 1.d0/b(n1i)
   sxx(n1i) = sxx(n1i)*binv
   svv(n1i) = a(n1i)*binv

   do i = n1i+1,n1f
   	gam(i) = c(i-1)*binv
   	binv   = 1.d0/(b(i)-a(i)*gam(i))
   	sxx(i) = (sxx(i)-a(i)*sxx(i-1))*binv
   	svv(i) = -a(i)*svv(i-1)*binv
   enddo

   sww(n1f) = c(n1f)*binv

	!>-------------------------------------------------------------------
	!!  Backward elimination
	!!-------------------------------------------------------------------
   do i = n1f-1,n1i,-1
   	sxx(i) = sxx(i)-gam(i+1)*sxx(i+1)
   	svv(i) = svv(i)-gam(i+1)*svv(i+1)
   	sww(i) = -gam(i+1)*sww(i+1)
   enddo

end subroutine tdma_pdd_x

!>======================================================================
subroutine pdd_comm1_x(sxxs,svvs,sxxr,svvr)

   implicit none
   include 'mpif.h'

   real(kind=cgreal), dimension(n2i:n2f,n3i:n3f), intent(in)  :: sxxs, svvs
   real(kind=cgreal), dimension(n2i:n2f,n3i:n3f), intent(out) :: sxxr, svvr
	!-------------------------------------
   integer :: status(mpi_status_size)
   integer :: reqs1, reqs2, reqr1, reqr2

   !>--------------------------------------
   !! Communication from i to i-1 rank
   !!--------------------------------------
   if (n1i == 1) then
   	call MPI_Irecv(sxxr,1,newtype1,rank_new_x+1,11,COMM_NEW_X,reqr1,ierr)
   	call MPI_Irecv(svvr,1,newtype1,rank_new_x+1,12,COMM_NEW_X,reqr2,ierr)
   	call MPI_Wait(reqr1,status,ierr)
   	call MPI_Wait(reqr2,status,ierr)

   elseif (n1f == nx) then
   	call MPI_ISEND(sxxs,1,newtype1,rank_new_x-1,11,COMM_NEW_X,reqs1,ierr)
   	call MPI_ISEND(svvs,1,newtype1,rank_new_x-1,12,COMM_NEW_X,reqs2,ierr)
   	call MPI_WAIT(reqs1,status,ierr)
   	call MPI_WAIT(reqs2,status,ierr)

   else
   	call MPI_ISEND(sxxs,1,newtype1,rank_new_x-1,11,COMM_NEW_X,reqs1,ierr)
   	call MPI_ISEND(svvs,1,newtype1,rank_new_x-1,12,COMM_NEW_X,reqs2,ierr)
   	call MPI_IRECV(sxxr,1,newtype1,rank_new_x+1,11,COMM_NEW_X,reqr1,ierr)
   	call MPI_IRECV(svvr,1,newtype1,rank_new_x+1,12,COMM_NEW_X,reqr2,ierr)
   	call MPI_WAIT(reqs1,status,ierr)
   	call MPI_WAIT(reqs2,status,ierr)
   	call MPI_WAIT(reqr1,status,ierr)
   	call MPI_WAIT(reqr2,status,ierr)

   endif

end subroutine pdd_comm1_x

!>======================================================================
subroutine pdd_ycal_x(sxxr,svvr,sww,sxx,syv,syw)

	implicit none
	include 'mpif.h'
	
	real(kind=cgreal), dimension(n2i:n2f,n3i:n3f), intent(in)  :: sxxr, svvr, sww,sxx
	real(kind=cgreal), dimension(n2i:n2f,n3i:n3f), intent(out) :: syv, syw
	!-------------------------------------
	integer :: jj, kk
	integer :: status(mpi_status_size)

	!---------------------------------------
	!  Calculation of y  
	!---------------------------------------
	if (n1f /= nx) then
		do kk=n3i,n3f;do jj=n2i,n2f
			syw(jj,kk) = (sxxr(jj,kk)-svvr(jj,kk)*sxx(jj,kk))/(1.d0-sww(jj,kk)*svvr(jj,kk))
			syv(jj,kk) = sxx(jj,kk)-sww(jj,kk)*syw(jj,kk)
		enddo;enddo
	endif

end subroutine pdd_ycal_x

!>======================================================================
subroutine pdd_comm2_x(syv,syv2)

   implicit none
   include 'mpif.h'

   integer :: status(MPI_STATUS_SIZE)

   real(kind=cgreal),dimension(n2i:n2f,n3i:n3f),intent(in)  :: syv
   real(kind=cgreal),dimension(n2i:n2f,n3i:n3f),intent(out) :: syv2
	!-------------------------------------
   integer :: reqs1,reqr1

   !---------------------------------------
   !  Communication from i to i+1 rank
   !---------------------------------------
   if (n1i == 1) then
   	call MPI_ISEND(syv,1,newtype1,rank_new_x+1,13,COMM_NEW_X,reqs1,ierr)
   	call MPI_WAIT(reqs1,status,ierr)

   elseif (n1f == nx) then
   	call MPI_IRECV(syv2,1,newtype1,rank_new_x-1,13,COMM_NEW_X,reqr1,ierr)
   	call MPI_WAIT(reqr1,status,ierr)

   else
   	call MPI_ISEND( syv,1,newtype1,rank_new_x+1,13,COMM_NEW_X,reqs1,ierr)
   	call MPI_IRECV(syv2,1,newtype1,rank_new_x-1,13,COMM_NEW_X,reqr1,ierr)
   	call MPI_WAIT(reqs1,status,ierr)
   	call MPI_WAIT(reqr1,status,ierr)

   endif

end subroutine pdd_comm2_x

!>======================================================================
subroutine pdd_subs_x(sxx_arr,svv_arr,sww_arr,syv2,syw)


   implicit none
   include 'mpif.h'

   integer :: status(mpi_status_size)

   !inputs
   real(kind=cgreal),dimension(n1i:n1f,n2i:n2f,n3i:n3f),intent(inout) :: sxx_arr
   real(kind=cgreal),dimension(n1i:n1f,n2i:n2f,n3i:n3f),intent(in)    :: svv_arr,sww_arr
   real(kind=cgreal),dimension(n2i:n2f,n3i:n3f),intent(in) :: syv2,syw
	!-------------------------------------
   real(kind=cgreal) :: ryv2,ryw
   integer :: ii,jj,kk 

   do kk=n3i,n3f
   	do jj=n2i,n2f

    		ryv2 = syv2(jj,kk);ryw = syw(jj,kk)
    		do ii=n1i,n1f
      		sxx_arr(ii,jj,kk) = sxx_arr(ii,jj,kk) - svv_arr(ii,jj,kk)*ryv2 - sww_arr(ii,jj,kk)*ryw
    		enddo!ii

   	enddo!jj
   enddo!kk

end subroutine pdd_subs_x

!>======================================================================
subroutine tdma_pdd_y(a,b,c,sxx,svv,sww)


    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: status(MPI_STATUS_SIZE)

    !Inputs
    REAL(KIND=CGREAL),DIMENSION(n2i:n2f),INTENT(IN) :: a,b,c
    REAL(KIND=CGREAL),DIMENSION(n2i:n2f),INTENT(INOUT) :: sxx
    REAL(KIND=CGREAL),DIMENSION(n2i:n2f),INTENT(OUT)   :: svv,sww
    !local
    REAL(KIND=CGREAL)                              :: binv
    INTEGER :: j

!--------------------------------------------------------------------
!   Forward elimination
!--------------------------------------------------------------------
    binv  = 1.d0/b(n2i)
    sxx(n2i) = sxx(n2i)*binv
    svv(n2i) = a(n2i)*binv
    do j=n2i+1,n2f
      gam(j) = c(j-1)*binv
      binv   = 1.d0/(b(j)-a(j)*gam(j))
      sxx(j) = (sxx(j)-a(j)*sxx(j-1))*binv
      svv(j) = -a(j)*svv(j-1)*binv
    enddo
    sww(n2f) = c(n2f)*binv

!--------------------------------------------------------------------
!   Backward elimination
!--------------------------------------------------------------------
    do j=n2f-1,n2i,-1
      sxx(j) = sxx(j)-gam(j+1)*sxx(j+1)
      svv(j) = svv(j)-gam(j+1)*svv(j+1)
      sww(j) = -gam(j+1)*sww(j+1)
    enddo

end subroutine tdma_pdd_y

!>======================================================================
subroutine pdd_comm1_y(sxxs,svvs,sxxr,svvr)


    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: status(MPI_STATUS_SIZE)

    !Inputs
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n3i:n3f),INTENT(IN)  :: sxxs,svvs
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n3i:n3f),INTENT(OUT) :: sxxr,svvr
    INTEGER :: reqs1,reqs2,reqr1,reqr2
    !---------------------------------------
    !  Communication from i to i-1 rank
    !---------------------------------------
    if (n2i == 1) then
 
     call MPI_IRECV(sxxr,1,newtype2,rank_new_y+1,11,COMM_NEW_Y,reqr1,ierr)
     call MPI_IRECV(svvr,1,newtype2,rank_new_y+1,12,COMM_NEW_Y,reqr2,ierr)
     call MPI_WAIT(reqr1,status,ierr)
     call MPI_WAIT(reqr2,status,ierr)

    else if (n2f == ny) then

     call MPI_ISEND(sxxs,1,newtype2,rank_new_y-1,11,COMM_NEW_Y,reqs1,ierr)
     call MPI_ISEND(svvs,1,newtype2,rank_new_y-1,12,COMM_NEW_Y,reqs2,ierr)
     call MPI_WAIT(reqs1,status,ierr)
     call MPI_WAIT(reqs2,status,ierr)

    else

     call MPI_ISEND(sxxs,1,newtype2,rank_new_y-1,11,COMM_NEW_Y,reqs1,ierr)
     call MPI_ISEND(svvs,1,newtype2,rank_new_y-1,12,COMM_NEW_Y,reqs2,ierr)
     call MPI_IRECV(sxxr,1,newtype2,rank_new_y+1,11,COMM_NEW_Y,reqr1,ierr)
     call MPI_IRECV(svvr,1,newtype2,rank_new_y+1,12,COMM_NEW_Y,reqr2,ierr)
     call MPI_WAIT(reqs1,status,ierr)
     call MPI_WAIT(reqs2,status,ierr)
     call MPI_WAIT(reqr1,status,ierr)
     call MPI_WAIT(reqr2,status,ierr)

    end if

end subroutine pdd_comm1_y

!>======================================================================
subroutine pdd_ycal_y(sxxr,svvr,sww,sxx,syv,syw)


    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: status(MPI_STATUS_SIZE)

    !Inputs
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n3i:n3f),INTENT(IN)  :: sxxr,svvr,sww,sxx
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n3i:n3f),INTENT(OUT) :: syv,syw
    INTEGER :: ii,kk
    !---------------------------------------
    !  Calculation of y  
    !---------------------------------------
    IF(n2f /= ny)then
     do kk=n3i,n3f;do ii=n1i,n1f
      syw(ii,kk) = (sxxr(ii,kk)-svvr(ii,kk)*sxx(ii,kk))/(1.d0-sww(ii,kk)*svvr(ii,kk))
      syv(ii,kk) = sxx(ii,kk)-sww(ii,kk)*syw(ii,kk)
     enddo;enddo
    endif

end subroutine pdd_ycal_y

!>======================================================================
subroutine pdd_comm2_y(syv,syv2)


    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: status(MPI_STATUS_SIZE)

    !Inputs
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n3i:n3f),INTENT(IN)  :: syv
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n3i:n3f),INTENT(OUT) :: syv2
    INTEGER :: reqs1,reqr1
    !---------------------------------------
    !  Communication from i to i+1 rank
    !---------------------------------------
    IF(n2i == 1)then

     call MPI_ISEND(syv,1,newtype2,rank_new_y+1,13,COMM_NEW_Y,reqs1,ierr)
     call MPI_WAIT(reqs1,status,ierr)

    elseif(n2f == ny)then

     call MPI_IRECV(syv2,1,newtype2,rank_new_y-1,13,COMM_NEW_Y,reqr1,ierr)
     call MPI_WAIT(reqr1,status,ierr)

    else

     call MPI_ISEND(syv,1,newtype2,rank_new_y+1,13,COMM_NEW_Y,reqs1,ierr)
     call MPI_IRECV(syv2,1,newtype2,rank_new_y-1,13,COMM_NEW_Y,reqr1,ierr)
     call MPI_WAIT(reqs1,status,ierr)
     call MPI_WAIT(reqr1,status,ierr)

    endif

end subroutine pdd_comm2_y

!>======================================================================
subroutine pdd_subs_y(sxx_arr,svv_arr,sww_arr,syv2,syw)


    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: status(MPI_STATUS_SIZE)

    !Inputs
    REAL(KIND=CGREAL),DIMENSION(n2i:n2f,n1i:n1f,n3i:n3f),INTENT(INOUT) :: sxx_arr
    REAL(KIND=CGREAL),DIMENSION(n2i:n2f,n1i:n1f,n3i:n3f),INTENT(IN)    :: svv_arr,sww_arr
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n3i:n3f),INTENT(IN) :: syv2,syw
    REAL(KIND=CGREAL) :: ryv2,ryw
    INTEGER :: ii,jj,kk 


    do kk=n3i,n3f
    do ii=n1i,n1f

     ryv2 = syv2(ii,kk);ryw = syw(ii,kk)
     do jj=n2i,n2f
       sxx_arr(jj,ii,kk) = sxx_arr(jj,ii,kk) - svv_arr(jj,ii,kk)*ryv2 - sww_arr(jj,ii,kk)*ryw
     enddo!jj

    enddo!ii
    enddo!kk

end subroutine pdd_subs_y

!>======================================================================
subroutine tdma_pdd_z(a,b,c,sxx,svv,sww)


    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: status(MPI_STATUS_SIZE)

    !Inputs
    REAL(KIND=CGREAL),DIMENSION(n3i:n3f),INTENT(IN)    :: a,b,c
    REAL(KIND=CGREAL),DIMENSION(n3i:n3f),INTENT(INOUT) :: sxx
    REAL(KIND=CGREAL),DIMENSION(n3i:n3f),INTENT(OUT)   :: svv,sww
    !local
    REAL(KIND=CGREAL)                              :: binv
    INTEGER :: k

!--------------------------------------------------------------------
!   Forward elimination
!--------------------------------------------------------------------
    binv = 1.d0/b(n3i)
    sxx(n3i) = sxx(n3i)*binv
    svv(n3i) = a(n3i)*binv
    do k=n3i+1,n3f
      gam(k) = c(k-1)*binv
      binv   = 1.d0/(b(k)-a(k)*gam(k))
      sxx(k) = (sxx(k)-a(k)*sxx(k-1))*binv
      svv(k) = -a(k)*svv(k-1)*binv
    enddo
    sww(n3f) = c(n3f)*binv

!--------------------------------------------------------------------
!   Backward elimination
!--------------------------------------------------------------------
    do k=n3f-1,n3i,-1
      sxx(k) = sxx(k)-gam(k+1)*sxx(k+1)
      svv(k) = svv(k)-gam(k+1)*svv(k+1)
      sww(k) = -gam(k+1)*sww(k+1)
    enddo

end subroutine tdma_pdd_z

!>======================================================================
subroutine pdd_comm1_z(sxxs,svvs,sxxr,svvr)


    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: status(MPI_STATUS_SIZE)

    !Inputs
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n2i:n2f),INTENT(IN)  :: sxxs,svvs
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n2i:n2f),INTENT(OUT) :: sxxr,svvr
    INTEGER :: reqs1,reqs2,reqr1,reqr2
    !---------------------------------------
    !  Communication from i to i-1 rank
    !---------------------------------------
    if (n3i == 1)then
 
     call MPI_IRECV(sxxr,1,newtype3,rank_new_z+1,11,COMM_NEW_Z,reqr1,ierr)
     call MPI_IRECV(svvr,1,newtype3,rank_new_z+1,12,COMM_NEW_Z,reqr2,ierr)
     call MPI_WAIT(reqr1,status,ierr)
     call MPI_WAIT(reqr2,status,ierr)

    elseif(n3f == nz)then

     call MPI_ISEND(sxxs,1,newtype3,rank_new_z-1,11,COMM_NEW_Z,reqs1,ierr)
     call MPI_ISEND(svvs,1,newtype3,rank_new_z-1,12,COMM_NEW_Z,reqs2,ierr)
     call MPI_WAIT(reqs1,status,ierr)
     call MPI_WAIT(reqs2,status,ierr)

    else

     call MPI_ISEND(sxxs,1,newtype3,rank_new_z-1,11,COMM_NEW_Z,reqs1,ierr)
     call MPI_ISEND(svvs,1,newtype3,rank_new_z-1,12,COMM_NEW_Z,reqs2,ierr)
     call MPI_IRECV(sxxr,1,newtype3,rank_new_z+1,11,COMM_NEW_Z,reqr1,ierr)
     call MPI_IRECV(svvr,1,newtype3,rank_new_z+1,12,COMM_NEW_Z,reqr2,ierr)
     call MPI_WAIT(reqs1,status,ierr)
     call MPI_WAIT(reqs2,status,ierr)
     call MPI_WAIT(reqr1,status,ierr)
     call MPI_WAIT(reqr2,status,ierr)

    end if

end subroutine pdd_comm1_z

!>======================================================================
subroutine pdd_ycal_z(sxxr,svvr,sww,sxx,syv,syw)


    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: status(MPI_STATUS_SIZE)

    !Inputs
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n2i:n2f),INTENT(IN)  :: sxxr,svvr,sww,sxx
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n2i:n2f),INTENT(OUT) :: syv,syw
    INTEGER :: ii,jj
    !---------------------------------------
    !  Calculation of y   !n3f
    !---------------------------------------
    IF(n3f /= nz)then
     do jj=n2i,n2f;do ii=n1i,n1f
      syw(ii,jj) = (sxxr(ii,jj)-svvr(ii,jj)*sxx(ii,jj))/(1.d0-sww(ii,jj)*svvr(ii,jj))
      syv(ii,jj) = sxx(ii,jj)-sww(ii,jj)*syw(ii,jj)
     enddo;enddo
    endif

end subroutine pdd_ycal_z

!>======================================================================
subroutine pdd_comm2_z(syv,syv2)


	IMPLICIT NONE
	INCLUDE 'mpif.h'
	
	INTEGER :: status(MPI_STATUS_SIZE)
	
	!Inputs
	REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n2i:n2f),INTENT(IN)  :: syv
	REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n2i:n2f),INTENT(OUT) :: syv2
	INTEGER :: reqs1,reqr1
	!---------------------------------------
	!  Communication from i to i+1 rank
	!---------------------------------------
	if(n3i == 1) then
		call MPI_ISEND(syv,1,newtype3,rank_new_z+1,13,COMM_NEW_Z,reqs1,ierr)
		call MPI_WAIT(reqs1,status,ierr)
	
	elseif(n3f == nz) then
		call MPI_IRECV(syv2,1,newtype3,rank_new_z-1,13,COMM_NEW_Z,reqr1,ierr)
		call MPI_WAIT(reqr1,status,ierr)
	
	else
		call MPI_ISEND(syv,1,newtype3,rank_new_z+1,13,COMM_NEW_Z,reqs1,ierr)
		call MPI_IRECV(syv2,1,newtype3,rank_new_z-1,13,COMM_NEW_Z,reqr1,ierr)
		call MPI_WAIT(reqs1,status,ierr)
		call MPI_WAIT(reqr1,status,ierr)
	
	endif

end subroutine pdd_comm2_z

!>======================================================================
subroutine pdd_subs_z(sxx_arr,svv_arr,sww_arr,syv2,syw)


    IMPLICIT NONE
    INCLUDE 'mpif.h'


    INTEGER :: status(MPI_STATUS_SIZE)

    !Inputs
    REAL(KIND=CGREAL),DIMENSION(n3i:n3f,n1i:n1f,n2i:n2f),INTENT(INOUT) :: sxx_arr
    REAL(KIND=CGREAL),DIMENSION(n3i:n3f,n1i:n1f,n2i:n2f),INTENT(IN)    :: svv_arr,sww_arr
    REAL(KIND=CGREAL),DIMENSION(n1i:n1f,n2i:n2f),INTENT(IN) :: syv2,syw
    REAL(KIND=CGREAL) :: ryv2,ryw
    INTEGER :: ii,jj,kk,counter,vcounter 


    do jj=n2i,n2f
    do ii=n1i,n1f

     ryv2 = syv2(ii,jj);ryw = syw(ii,jj)
     do kk=n3i,n3f
       sxx_arr(kk,ii,jj) = sxx_arr(kk,ii,jj) - svv_arr(kk,ii,jj)*ryv2 - sww_arr(kk,ii,jj)*ryw
     enddo!kk

    enddo!ii
    enddo!jj

end subroutine pdd_subs_z

!>======================================================================
end module m_tdma_pdd






