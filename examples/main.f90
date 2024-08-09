!>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!
!! 						Tri-Diagonal matrice Algorithm.
!! 						   Last updated: Apr 16, 2024
!!
!!
!!                                           By Seungchan Kim, Jihoo Kim
!!                                  Department of Mechanical Engineering
!!                                                               POSTECH
!!                                         Email: mack1029@postech.ac.kr 
!!        
!!        
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main

   use pdd_host
   use global_host
   use m_global_initialize
   use m_pdd_initialize
  	use m_tdma_main
	use m_communication

   implicit none
   include 'mpif.h'
   real :: t1, t2, maxerr


	!>-- MPI Initialize
   call MPI_INIT(ierr)
   call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
   call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

	!>-- Initialize variable
   call read_infile   !<-- src subroutine
	call malloc_pdd    !<-- src subroutine
	call malloc_global 
   call mpi_topology  !<-- src subroutine
	if (rank == 0) print*,''
	if (rank == 0) write(*,319) 

	!>-- Initialize 3D data (generate and distribute data from rank 0)
	t1 = MPI_Wtime()
	call init_3Ddata   !<-- contain src subroutine
	t2 = MPI_Wtime()
	if (rank == 0) write(*,320) t2-t1

	!>-- TDMA algorithm
	t1 = MPI_Wtime()
	call tdma_3d_main  !<-- src subroutine
	t2 = MPI_Wtime()
	if (rank == 0) write(*,321) t2-t1

	!>-- Validation (gather all data in rank 0)
	t1 = MPI_Wtime()
	call validation(maxerr)  !<-- contain src subroutine
	t2 = MPI_Wtime()
	if (rank == 0) then
		write(*,322) t2-t1
		print*,''
		write(*,323) maxerr
		print*,''
	endif

	!>-- Free memory and finalize MPI
   call free_pdd      !<-- src subroutine
   call free_global
   call MPI_Finalize( ierr )


 319 format('Time calculation')
 320 format(' Time for scattering data =',f10.6)
 321 format(' Time for TDMA solving =',f10.6)
 322 format(' Time for validating result =',f10.6)
 323 format('Maximum error =',e21.12)
    
end program main
