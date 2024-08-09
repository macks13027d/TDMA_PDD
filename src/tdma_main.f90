module m_tdma_main

use m_tdma_pdd
contains
!==========================================================================================
subroutine tdma_3d_main

	use pdd_host

	implicit none
	include 'mpif.h'

	if (dir == 1) then
		call tdma_3d_x

	elseif (dir == 2) then
		call transpose_ij
		call tdma_3d_y
		call transpose_ji

	elseif (dir == 3) then
		call transpose_ik
		call tdma_3d_z
		call transpose_ki

	endif

	call MPI_Barrier(MPI_COMM_WORLD,ierr)

end subroutine

!==========================================================================================
subroutine tdma_3d_x

	use pdd_host
	
	implicit none
	include 'mpif.h'
	
	integer :: i, j, k 
	real :: t1, t2, t_tridag
	

	call MPI_Barrier(MPI_COMM_WORLD,ierr)
	t1 = MPI_Wtime()
   !>------------------------------------------
   !! If # of x-dir processors = 1
   !!------------------------------------------
	if (nproc_x == 1) then
		do k = n3i,n3f
			do j = n2i,n2f

				!>-- Solve Eq
				xx_i(n1i:n1f) = dat0(n1i:n1f,j,k)
				call tdma(mat_a(:,j,k),mat_b(:,j,k),mat_c(:,j,k),xx_i,sol_tmp_i,n1i,n1f)
if(k==1 .and. j== 1)print*,xx_i-sol_tmp_i,'vvv'

				!>-- Substitute solution
				dat0(n1i:n1f,j,k) = sol_tmp_i(n1i:n1f)

			enddo !<--j
		enddo !<--k

   !>------------------------------------------
   !! If # of x-dir processors > 1
   !!------------------------------------------
   else 
		do k=n3i,n3f
			do j=n2i,n2f

			!>-- Solve Eq
			xx_i(n1i:n1f) = dat0(n1i:n1f,j,k)
			call tdma_pdd_x(mat_a(:,j,k),mat_b(:,j,k),mat_c(:,j,k),xx_i,vv_i,ww_i)

			!>-- Substitute solution
			dat0(n1i:n1f,j,k) = xx_i(n1i:n1f)
			vv_arr_i(n1i:n1f,j,k) = vv_i(n1i:n1f)
			ww_arr_i(n1i:n1f,j,k) = ww_i(n1i:n1f)
			xxs_i(j,k) = xx_i(n1i)
			xxe_i(j,k) = xx_i(n1f)
			vvs_i(j,k) = vv_i(n1i)
			wwe_i(j,k) = ww_i(n1f)

			enddo !<--j
		enddo !<--k
	end if

   t2 = MPI_WTIME()
   t_tridag = t2 - t1


   if(nproc_x /= 1)then      
		!>------------------------------------------
		!! Communication between processors
		!!------------------------------------------
   	t1 = MPI_Wtime()

		call pdd_comm1_x(xxs_i,vvs_i,xxr_i,vvr_i)
		call pdd_ycal_x(xxr_i,vvr_i,wwe_i,xxe_i,yv_i,yw_i)
		call pdd_comm2_x(yv_i,yv2_i)

		t2 = MPI_WTIME()
  		if (rank == 0) print*,'    TDMA (communication) : ', t2 - t1

		!>------------------------------------------
		!! Inversion of tridiagonal matrix
		!!------------------------------------------
		t1 = MPI_Wtime()

		call pdd_subs_x(dat0,vv_arr_i,ww_arr_i,yv2_i,yw_i)

		t2 = MPI_WTIME()
   	t_tridag = t_tridag + t2 - t1

   endif

	!>-- Print 
	if (rank == 0) print*,'    TDMA (tridiagonal matrix inversion) : ', t_tridag

   return
end subroutine tdma_3d_x 

!==========================================================================================
subroutine tdma_3d_y

   use pdd_host

   implicit none
   include 'mpif.h'
   
   integer :: i, j, k
	real :: t1, t2, t_tridag


   !>------------------------------------------
   !! If # of y-dir processors = 1
   !!------------------------------------------
	t1 = MPI_Wtime()

   if (nproc_y == 1) then
   	do k = n3i,n3f
   		do i = n1i,n1f

				!>-- Solve Eq
        		xx_j(n2i:n2f) = dat1_j(n2i:n2f,i,k)
         	call tdma(mat_a(i,:,k),mat_b(i,:,k),mat_c(i,:,k),xx_j,sol_tmp_j,n2i,n2f)

				!>-- Substitute solution
				dat1_j(n2i:n2f,i,k) = sol_tmp_j(n2i:n2f)

			enddo !<--i
		enddo !<--k

   !>------------------------------------------
   !! If # of y-dir processors > 1
   !!------------------------------------------
   else 
   	do k = n3i,n3f
   		do i = n1i,n1f

			!>-- Solve Eq
        	xx_j(n2i:n2f) = dat1_j(n2i:n2f,i,k)
        	call tdma_pdd_y(mat_a(i,:,k),mat_b(i,:,k),mat_c(i,:,k),xx_j,vv_j,ww_j)


			!>-- Substitute solution
        	dat1_j(n2i:n2f,i,k) = xx_j(n2i:n2f)
        	vv_arr_j(n2i:n2f,i,k) = vv_j(n2i:n2f)
        	ww_arr_j(n2i:n2f,i,k) = ww_j(n2i:n2f)

        	xxs_j(i,k) = xx_j(n2i);xxe_j(i,k) = xx_j(n2f)
        	vvs_j(i,k) = vv_j(n2i);wwe_j(i,k) = ww_j(n2f)

			enddo !<--i
		enddo !<--k
	endif

   t2 = MPI_WTIME()
   t_tridag = t2 - t1

   if(nproc_y /= 1) then
		!>------------------------------------------
		!! Communication between processors
		!!------------------------------------------
   	t1 = MPI_Wtime()

      call pdd_comm1_y(xxs_j,vvs_j,xxr_j,vvr_j)
      call pdd_ycal_y(xxr_j,vvr_j,wwe_j,xxe_j,yv_j,yw_j)
      call pdd_comm2_y(yv_j,yv2_j)

		t2 = MPI_WTIME()
  		if (rank == 0) print*,'    TDMA (communication) : ', t2 - t1

		!>------------------------------------------
		!! Inversion of tridiagonal matrix
		!!------------------------------------------
		t1 = MPI_Wtime()

		call pdd_subs_y(dat1_j,vv_arr_j,ww_arr_j,yv2_j,yw_j)

		t2 = MPI_WTIME()
   	t_tridag = t_tridag + t2 - t1

   endif
      
	!>-- Print 
	if (rank == 0) print*,'    TDMA (tridiagonal matrix inversion) : ', t_tridag

   return
end subroutine tdma_3d_y

!==========================================================================================
subroutine tdma_3d_z

   use pdd_host

   implicit none
   include 'mpif.h'
   
   integer :: i, j, k
	real :: t1, t2, t_tridag

   !>------------------------------------------
   !! If # of z-dir processors = 1
   !!------------------------------------------
	t1 = MPI_Wtime()

   if (nproc_z == 1) then
   	do j = n2i,n2f
   		do i = n1i,n1f

				!>-- Solve Eq
        		xx_k(n3i:n3f) = dat1_k(n3i:n3f,i,j)
           	call tdma(mat_a(i,j,:),mat_b(i,j,:),mat_c(i,j,:),xx_k,sol_tmp_k,n3i,n3f)

				!>-- Substitute solution
        		dat1_k(n3i:n3f,i,j) = sol_tmp_k(n3i:n3f)

			enddo !<--i
		enddo !<--j

   !>------------------------------------------
   !! If # of z-dir processors > 1
   !!------------------------------------------
   else 
   	do j = n2i,n2f
   		do i = n1i,n1f

			!>-- Solve Eq
        	xx_k(n3i:n3f) = dat1_k(n3i:n3f,i,j)
         call tdma_pdd_z(mat_a(i,j,:),mat_b(i,j,:),mat_c(i,j,:),xx_k,vv_k,ww_k)

			!>-- Substitute solution
        	dat1_k(n3i:n3f,i,j) = xx_k(n3i:n3f)
        	vv_arr_k(n3i:n3f,i,j) = vv_k(n3i:n3f)
        	ww_arr_k(n3i:n3f,i,j) = ww_k(n3i:n3f)
        	xxs_k(i,j) = xx_k(n3i);xxe_k(i,j) = xx_k(n3f)
        	vvs_k(i,j) = vv_k(n3i);wwe_k(i,j) = ww_k(n3f)

			enddo !<--i
		enddo !<--j
	endif

   t2 = MPI_WTIME()
   t_tridag = t2 - t1


   IF(nproc_z /= 1)THEN
		!>------------------------------------------
		!! Communication between processors
		!!------------------------------------------
   	t1 = MPI_Wtime()

   	call pdd_comm1_z(xxs_k,vvs_k,xxr_k,vvr_k)
   	call pdd_ycal_z(xxr_k,vvr_k,wwe_k,xxe_k,yv_k,yw_k)
   	call pdd_comm2_z(yv_k,yv2_k)

		t2 = MPI_WTIME()
  		if (rank == 0) print*,'    TDMA (communication) : ', t2 - t1

		!>------------------------------------------
		!! Inversion of tridiagonal matrix
		!!------------------------------------------
		t1 = MPI_Wtime()

      call pdd_subs_z(dat1_k,vv_arr_k,ww_arr_k,yv2_k,yw_k)


		t2 = MPI_WTIME()
   	t_tridag = t_tridag + t2 - t1

   endif

	!>-- Print 
	if (rank == 0) print*,'    TDMA (tridiagonal matrix inversion) : ', t_tridag

   return
end subroutine tdma_3d_z

!==========================================================================================
subroutine transpose_ij

   use pdd_host

   implicit none
   include 'mpif.h'
   
   integer :: i, j, k
	real :: t1, t2, t_trans


   !>------------------------------------------
   !! Matrix transpose
   !!------------------------------------------
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
	t1 = MPI_Wtime()

   do k = n3i,n3f
		do i = n1i,n1f
     		dat1_j(n2i:n2f,i,k) = dat0(i,n2i:n2f,k)
   	enddo
	enddo

	t2 = MPI_WTIME()
   t_trans = t2 - t1

	!>-- Print 
	if (rank == 0) print*,'    TDMA (transpose matrix (i,j,k)->(j,i,k)) : ', t_trans

   return
end subroutine transpose_ij

!==========================================================================================
subroutine transpose_ji

   use pdd_host

   implicit none
   include 'mpif.h'
   
   integer :: i, j, k
	real :: t1, t2, t_trans


   !>------------------------------------------
   !! Matrix transpose
   !!------------------------------------------
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
	t1 = MPI_Wtime()

	do k = n3i,n3f
		do j = n2i,n2f
			dat0(n1i:n1f,j,k) = dat1_j(j,n1i:n1f,k)
		enddo
	enddo

	t2 = MPI_WTIME()
   t_trans = t2 - t1

	!>-- Print 
	if (rank == 0) print*,'    TDMA (transpose matrix (j,i,k)->(i,j,k)) : ', t_trans

   return
end subroutine transpose_ji

!==========================================================================================
subroutine transpose_ik

   use pdd_host

   implicit none
   include 'mpif.h'
   
   integer :: i, j, k
	real :: t1, t2, t_trans


   !>------------------------------------------
   !! Matrix transpose
   !!------------------------------------------
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
	t1 = MPI_Wtime()

   do j = n2i,n2f
		do i = n1i,n1f
      	dat1_k(n3i:n3f,i,j) = dat0(i,j,n3i:n3f)
   	enddo
	enddo

	t2 = MPI_WTIME()
   t_trans = t2 - t1

	!>-- Print 
	if (rank == 0) print*,'    TDMA (transpose matrix (i,j,k)->(k,i,j)) : ', t_trans

   return
end subroutine transpose_ik

!==========================================================================================
subroutine transpose_ki

   use pdd_host

   implicit none
   include 'mpif.h'
   
   integer :: i, j, k
	real :: t1, t2, t_trans


   !>------------------------------------------
   !! Matrix transpose
   !!------------------------------------------
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
	t1 = MPI_Wtime()

	do k = n3i,n3f
		do j = n2i,n2f
      	dat0(n1i:n1f,j,k) = dat1_k(k,n1i:n1f,j)
   	enddo
	enddo

	t2 = MPI_WTIME()
   t_trans = t2 - t1

	!>-- Print 
	if (rank == 0) print*,'    TDMA (transpose matrix (k,i,j)->(i,j,k)) : ', t_trans

   return
end subroutine transpose_ki

!==========================================================================================
subroutine transpose_jk

   use pdd_host

   implicit none
   include 'mpif.h'
   
   integer :: i, j, k
	real :: t1, t2, t_trans


   !>------------------------------------------
   !! Matrix transpose
   !!------------------------------------------
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
	t1 = MPI_Wtime()

   do j = n2i,n2f
		do i = n1i,n1f
      	dat1_k(n3i:n3f,i,j) = dat1_j(j,i,n3i:n3f)
   	enddo
	enddo

	t2 = MPI_WTIME()
   t_trans = t2 - t1

	!>-- Print 
	if (rank == 0) print*,'    TDMA (transpose matrix (j,i,k)->(k,i,j)) : ', t_trans

   return
end subroutine transpose_jk

!==========================================================================================
subroutine transpose_kj

   use pdd_host

   implicit none
   include 'mpif.h'
   
   integer :: i, j, k
	real :: t1, t2, t_trans


   !>------------------------------------------
   !! Matrix transpose
   !!------------------------------------------
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
	t1 = MPI_Wtime()

	do k = n3i,n3f
		do i = n1i,n1f
      	dat1_j(n2i:n2f,i,k) = dat1_k(k,i,n2i:n2f)
   	enddo
	enddo

	t2 = MPI_WTIME()
   t_trans = t2 - t1

	!>-- Print 
	if (rank == 0) print*,'    TDMA (transpose matrix (k,i,j)->(j,i,k)) : ', t_trans

   return
end subroutine transpose_kj

!==========================================================================================
end module m_tdma_main
