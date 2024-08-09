!  subroutine read_infile()
!  subroutine mpi_topology()

module m_pdd_initialize 
  
contains
!>=====================================================================
subroutine mpi_topology

   use pdd_host

   implicit none
   include 'mpif.h'

   integer :: i, oneslice, twoslice, sizeofreal, sizeofint
   integer :: nxloc_dump, nyloc_dump, nzloc_dump
   integer :: is_recv, ie_recv, js_recv, je_recv, ks_recv, ke_recv
   !>-- Sub topology
   logical,dimension(0:2) :: belongs
   integer :: newtype_dum1,newtype_dum2,newtype_dum3

	!>-- Check whether # of processora is correct
   if ( nproc_x*nproc_y*nproc_z /= nprocs ) then
     if ( rank == 0 ) then
       print*,""
       print*,"error:  p1*p2 =/= # of cpu's!"
       print*,"please correct your input file and run again"
       print*,""
     endif
     call mpi_finalize( ierr )
     stop
   endif

	!>===================================================================
	!! 						MPI communication parameter 1
	!!					collect data in rank 0 after calculation	
	!!===================================================================
	!>--Share domain index information among processors
	do i = 0,nprocs-1
	   if (i .eq. rank) then  !<-- just save when proceeding rank is equal to 'i' 
	   	is_neig(i) = n1i
	   	js_neig(i) = n2i
	   	ks_neig(i) = n3i
	   	ie_neig(i) = n1f
	   	je_neig(i) = n2f
	   	ke_neig(i) = n3f   
	   	
	   	is_recv = n1i
	   	js_recv = n2i
	   	ks_recv = n3i
	   	ie_recv = n1f
	   	je_recv = n2f
	   	ke_recv = n3f
	   	
	   	call MPI_BCAST(is_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(js_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(ks_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(ie_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(je_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(ke_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	
	   else  !<-- otherwise broadcast and save
	   	call MPI_BCAST(is_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(js_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(ks_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(ie_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(je_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	call MPI_BCAST(ke_recv,1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
	   	
	   	is_neig(i) = is_recv
	   	js_neig(i) = js_recv
	   	ks_neig(i) = ks_recv
	   	ie_neig(i) = ie_recv
	   	je_neig(i) = je_recv
	   	ke_neig(i) = ke_recv  
	     
	   endif
	enddo

	!>-- Define matrix type for communicating sub-domain information
	if (rank.eq.0) then
		
		do i = 1,nprocs-1
			nxloc_dump = ie_neig(i)-is_neig(i)+1
			nyloc_dump = je_neig(i)-js_neig(i)+1
			nzloc_dump = ke_neig(i)-ks_neig(i)+1
			
			call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,sizeofreal,ierr)
			call MPI_TYPE_VECTOR(nxloc_dump,1,1,MPI_DOUBLE_PRECISION,oneslice,ierr)
			call MPI_TYPE_HVECTOR(nyloc_dump,1,(nx)*sizeofreal,oneslice,twoslice,ierr)
			call MPI_TYPE_HVECTOR(nzloc_dump,1,(nx)*(ny)*sizeofreal,twoslice,thrDR_neig(i),ierr)
			call MPI_TYPE_COMMIT(thrDR_neig(i),ierr)
		enddo

	endif
    
	!>===================================================================
	!! 						MPI communication parameter 2
	!!					    recv and send during PDD algorithm	
	!!===================================================================
	!>--Kind 1
   if (nproc_x > 1) then
   	belongs(0) = .true.
   	belongs(1) = .false.
   	belongs(2) = .false.
   	call MPI_CART_SUB(comm_cart,belongs,comm_new_x,ierr)
   	call MPI_COMM_RANK(comm_new_x,rank_new_x,ierr)
   end if

   if (nproc_y > 1) then
   	belongs(1) = .true.
   	belongs(0) = .false.
   	belongs(2) = .false.
   	call MPI_CART_SUB(comm_cart,belongs,comm_new_Y,ierr)
   	call MPI_COMM_RANK(comm_new_Y,rank_new_y,ierr)
   end if

   if (nproc_z > 1) then
   	belongs(2) = .true.
   	belongs(0) = .false.
   	belongs(1) = .false.
   	call MPI_CART_SUB(comm_cart,belongs,comm_new_Z,ierr)
   	call MPI_COMM_RANK(comm_new_Z,rank_new_z,ierr)
   end if

	!>--Kind 2
   call MPI_TYPE_CONTIGUOUS((n2f-n2i+1),MPI_DOUBLE_PRECISION,newtype_dum1,ierr)
   call MPI_TYPE_CONTIGUOUS((n3f-n3i+1),newtype_dum1,newtype1,ierr) !x
   call MPI_TYPE_COMMIT(newtype1,ierr)

   call MPI_TYPE_CONTIGUOUS((n3f-n3i+1),MPI_DOUBLE_PRECISION,newtype_dum2,ierr)
   call MPI_TYPE_CONTIGUOUS((n1f-n1i+1),newtype_dum2,newtype2,ierr) !y
   call MPI_TYPE_COMMIT(newtype2,ierr)

   call MPI_TYPE_CONTIGUOUS((n1f-n1i+1),MPI_DOUBLE_PRECISION,newtype_dum3,ierr)
   call MPI_TYPE_CONTIGUOUS((n2f-n2i+1),newtype_dum3,newtype3,ierr) !z
   call MPI_TYPE_COMMIT(newtype3,ierr)

end subroutine mpi_topology

!>=============================================================================
subroutine matrix_extension

   use pdd_host

   implicit none
   include 'mpif.h'

   if (nproc_x-1.eq.coords(1).and.mod((nx),nproc_x).gt.0) then
	  mat_a(nx+1:n1f,n2i:n2f,n3i:n3f)=0
     mat_b(nx+1:n1f,n2i:n2f,n3i:n3f)=1
     mat_c(nx+1:n1f,n2i:n2f,n3i:n3f)=0
     dat0(nx+1:n1f,n2i:n2f,n3i:n3f)=0   
   endif

   if (nproc_y-1.eq.coords(2).and.mod((ny),nproc_y).gt.0) then
     mat_a(n1i:n1f,ny+1:n2f,n3i:n3f)=0
     mat_b(n1i:n1f,ny+1:n2f,n3i:n3f)=1
     mat_c(n1i:n1f,ny+1:n2f,n3i:n3f)=0
     dat0(n1i:n1f,ny+1:n2f,n3i:n3f)=0
   endif

   if (nproc_z-1.eq.coords(3).and.mod((nz),nproc_z).gt.0) then
     mat_a(n1i:n1f,n2i:n2f,nz+1:n3f)=0
     mat_b(n1i:n1f,n2i:n2f,nz+1:n3f)=1
     mat_c(n1i:n1f,n2i:n2f,nz+1:n3f)=0
     dat0(n1i:n1f,n2i:n2f,nz+1:n3f)=0
   endif

end subroutine matrix_extension

!>=============================================================================
end module m_pdd_initialize
