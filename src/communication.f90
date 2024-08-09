module m_communication

contains
!>=================================================================================
subroutine gather_main

	use pdd_host
	
	implicit none
	include 'mpif.h'
	
	if(rank == 0)then
		!>-- recieve data
		call gather_recv
	elseif(rank /= 0)then
		!>-- send data
		call gather_send
	endif
    
end subroutine gather_main

!>=================================================================================
subroutine gather_recv

	use pdd_host
	
	implicit none
	include 'mpif.h'
	INTEGER :: status(MPI_STATUS_SIZE)
	integer :: i
	
	glb_dat(n1i:n1f,n2i:n2f,n3i:n3f) = dat0

	do i = 1, nprocs-1
		call MPI_Recv(glb_dat_dump(is_neig(i),js_neig(i),ks_neig(i)),1 &
									 ,thrDR_neig(i),i,301,MPI_COMM_WORLD,status,ierr)
		glb_dat(is_neig(i):ie_neig(i),js_neig(i):je_neig(i),ks_neig(i):ke_neig(i)) = &
		  glb_dat_dump(is_neig(i):ie_neig(i),js_neig(i):je_neig(i),ks_neig(i):ke_neig(i))
	end do
    
end subroutine gather_recv

!>=================================================================================
subroutine gather_send

	use pdd_host
	
	implicit none
	include 'mpif.h'
	integer :: i
	
	call MPI_Send(dat0,(nxloc_mpi)*(nyloc_mpi)*(nzloc_mpi) & 
                     ,MPI_DOUBLE_PRECISION,0,301, MPI_COMM_WORLD,ierr)
    
end subroutine gather_send

!>=================================================================================
subroutine scatter_main

	use pdd_host
	
	implicit none
	include 'mpif.h'
	integer :: i, stype
	

	do i=1,nxloc_mpi
		do stype=1,4
			!>-- send data
			if (rank == 0)	call scatter_send(i,stype)
			!>-- recieve data
			call scatter_recv(i,stype)

			call MPI_Barrier(MPI_COMM_WORLD,ierr)
		enddo
	enddo
    
end subroutine scatter_main

!>=================================================================================
subroutine scatter_recv( cnt, stype )

	use pdd_host
	
	implicit none
	include 'mpif.h'
	integer :: status(MPI_STATUS_SIZE)
	integer, intent(in) :: cnt, stype
	
	if (rank == 0) then
		if(stype == 1) dat0(n1i+cnt-1,:,:) = glb_dat(n1i+cnt-1,n2i:n2f,n3i:n3f)
		if(stype == 2) mat_a(n1i+cnt-1,:,:) = glb_mat_a(n1i+cnt-1,n2i:n2f,n3i:n3f)
		if(stype == 3) mat_b(n1i+cnt-1,:,:) = glb_mat_b(n1i+cnt-1,n2i:n2f,n3i:n3f)
		if(stype == 4) mat_c(n1i+cnt-1,:,:) = glb_mat_c(n1i+cnt-1,n2i:n2f,n3i:n3f)

	elseif (rank /= 0 ) then
		call MPI_Recv(var_dump(n1i+cnt-1,:,:),(nyloc_mpi)*(nzloc_mpi) &
									 ,MPI_DOUBLE_PRECISION,0,rank,MPI_COMM_WORLD,status,ierr)

		if(stype == 1) dat0(n1i+cnt-1,:,:) = var_dump(n1i+cnt-1,:,:)
		if(stype == 2) mat_a(n1i+cnt-1,:,:) = var_dump(n1i+cnt-1,:,:)
		if(stype == 3) mat_b(n1i+cnt-1,:,:) = var_dump(n1i+cnt-1,:,:)
		if(stype == 4) mat_c(n1i+cnt-1,:,:) = var_dump(n1i+cnt-1,:,:)
	endif
    
end subroutine scatter_recv

!>=================================================================================
subroutine scatter_send( cnt, stype )

	use pdd_host
	
	implicit none
	include 'mpif.h'
	integer, intent(in) :: cnt
	integer :: i, stype
	
	do i = 1,nprocs-1
		if (stype == 1) then
			call MPI_Send(glb_dat(is_neig(i)+cnt-1, js_neig(i):je_neig(i), &
								ks_neig(i):ke_neig(i)), (nyloc_mpi)*(nzloc_mpi) & 
								,MPI_DOUBLE_PRECISION,i,i, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
		elseif (stype == 2) then
			call MPI_Send(glb_mat_a(is_neig(i)+cnt-1, js_neig(i):je_neig(i), &
								ks_neig(i):ke_neig(i)), (nyloc_mpi)*(nzloc_mpi) & 
								,MPI_DOUBLE_PRECISION,i,i, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
		elseif (stype == 3) then
			call MPI_Send(glb_mat_b(is_neig(i)+cnt-1, js_neig(i):je_neig(i), &
								ks_neig(i):ke_neig(i)), (nyloc_mpi)*(nzloc_mpi) & 
								,MPI_DOUBLE_PRECISION,i,i, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
		elseif (stype == 4) then
			call MPI_Send(glb_mat_c(is_neig(i)+cnt-1, js_neig(i):je_neig(i), &
								ks_neig(i):ke_neig(i)), (nyloc_mpi)*(nzloc_mpi) & 
								,MPI_DOUBLE_PRECISION,i,i, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
		endif
	enddo
    
end subroutine scatter_send
!>=================================================================================
end module m_communication
