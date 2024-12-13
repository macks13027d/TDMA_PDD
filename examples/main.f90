!=======================================================================================================================
!> @file        main.f90
!> @brief       Main module for running example codes
!> @details     This module calls the subroutines for running the example codes
!>					 First, it generates the 3D example solution and calculates the corresponding initial data. 
!>					 Then, the data is distributed to all processors, and the TDM is solved in the selected direction. 
!> 				 Finally, the results are gathered in a single processor and validated using the example solution generated at the initial step.
!>					 Also, as noted in communication.f90, please be aware of the memory size limitations when generating the 3D example solution.
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
   real :: time(3), t1, t2, maxerr


	!>-- MPI Initialize
   call MPI_INIT(ierr)
   call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
   call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

	!>-- Initialize variable
   call read_infile   !<-- src subroutine
	call malloc_pdd    !<-- src subroutine
	call malloc_global 
   call mpi_topology  !<-- src subroutine

	!>-- Initialize 3D data (generate and distribute data from rank 0)
	t1 = MPI_Wtime()
	call init_3Ddata   !<-- contain src subroutine
	t2 = MPI_Wtime()
	time(1) = t2-t1

	!>-- TDMA algorithm
	t1 = MPI_Wtime()
	call tdma_3d_main  !<-- src subroutine
	t2 = MPI_Wtime()
	time(2) = t2-t1

	!>-- Validation (gather all data in rank 0)
	t1 = MPI_Wtime()
	call validation(maxerr)  !<-- contain src subroutine
	t2 = MPI_Wtime()
	time(3) = t2-t1


	!>-- Print results
	if (rank == 0) then
		print *,'' 
		print *,&
		'------------------------------------------------------------------'
		write(*,'(a,e10.3,a)') ' Scattering data needs   ', time(1), ' s'
		write(*,'(a,e10.3,a)') ' Solving TDMA needs      ', time(2), ' s'
		write(*,'(a,e10.3,a)') ' Validationg result needs', time(3), ' s'
		print *,''
		write(*,'(a,e13.5)') ' Maximum error:', maxerr
		print *,''
		print *,& 
		'------------------------------------------------------------------'
		print *,''
	endif

	!>-- Free memory and finalize MPI
   call free_pdd      !<-- src subroutine
   call free_global
   call MPI_Finalize( ierr )

end program main
