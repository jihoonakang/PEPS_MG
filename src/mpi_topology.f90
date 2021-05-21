!======================================================================================================================
!> @file        mpi_topology.f90
!> @brief       This file contains a module of communication topology for the example problem of PEPS_MG.
!> @author      
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>
!> @date        May 2021
!> @version     1.0
!> @par         Copyright
!>              Copyright (c) 2021 Ji-Hoon Kang, Korea Institute of Science and Technology Information,
!>              All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE in )
!======================================================================================================================

!>
!> @brief       Module for creating the cartesian topology of the MPI processes and subcommunicators.
!> @details     This module has three subcommunicators in each-direction and related subroutines.
!>
module mpi_topology

    implicit none

    private

    integer(kind=4), public :: mpi_world_cart   !< Communicator for cartesian topology
    integer(kind=4), public :: nprocs, myrank   !< Number of processes and rank ID in MPI_COMM_WORLD
    integer(kind=4), public :: np_dim(0:2)      !< Number of MPI processes in 3D topology
    logical, public         :: period(0:2)      !< Periodicity in each direction

    !> @brief   Type variable for the information of 1D communicator
    type, public :: cart_comm_1d
        integer(kind=4), :: myrank              !< Rank ID in current communicator
        integer(kind=4), :: nprocs              !< Number of processes in current communicator
        integer(kind=4), :: west_rank           !< Previous rank ID in current communicator
        integer(kind=4), :: east_rank           !< Next rank ID in current communicator
        integer(kind=4), :: mpi_comm            !< Current communicator
    end type cart_comm_1d

    ! Target attribute for the change of aggregation dimension in CGPSA
    type(cart_comm_1d), public, target :: comm_1d_x     !< Subcommunicator information in x-direction
    type(cart_comm_1d), public, target :: comm_1d_y     !< Subcommunicator information in y-direction
    type(cart_comm_1d), public, target :: comm_1d_z     !< Subcommunicator information in z-direction

    public  :: mpi_topology_create
    public  :: mpi_topology_destroy

    contains

    !>
    !> @brief       Create the cartesian topology for the MPI processes and subcommunicators.
    !>
    subroutine mpi_topology_create()
        
        use mpi

        implicit none

        logical         :: remain(0:2)
        integer(kind=4) :: ierr

        ! Create the cartesian topology.
        call MPI_Cart_create( MPI_COMM_WORLD    &!  input  | integer(kind=4)      | Input communicator (handle).
                            , 3                 &!  input  | integer(kind=4)      | Number of dimensions of Cartesian grid (integer(kind=4)).
                            , np_dim            &!  input  | integer(kind=4)(1:3) | integer(kind=4) array of size ndims specifying the number of processes in each dimension.
                            , period            &!  input  | logical(1:3) | Logical array of size ndims specifying whether the grid is periodic (true=1) or not (false=0) in each dimension.
                            , .false.           &!  input  | logical      | Ranking may be reordered (true=1) or not (false=0) (logical).
                            , mpi_world_cart    &! *output | integer(kind=4)      | Communicator with new Cartesian topology (handle).
                            , ierr              &!  output | integer(kind=4)      | Fortran only: Error status
                            )

        ! Create subcommunicators and assign two neighboring processes in the x-direction.
        remain(0) = .true.
        remain(1) = .false.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_x%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x%mpi_comm, comm_1d_x%myrank, ierr)
        call MPI_Comm_size(comm_1d_x%mpi_comm, comm_1d_x%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_x%mpi_comm, 0, 1, comm_1d_x%west_rank, comm_1d_x%east_rank, ierr)
    
        ! Create subcommunicators and assign two neighboring processes in the y-direction.
        remain(0) = .false.
        remain(1) = .true.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_y%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_y%mpi_comm, comm_1d_y%myrank, ierr)
        call MPI_Comm_size(comm_1d_y%mpi_comm, comm_1d_y%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_y%mpi_comm, 0, 1, comm_1d_y%west_rank, comm_1d_y%east_rank, ierr)

        ! Create subcommunicators and assign two neighboring processes in the z-direction.
        remain(0) = .false.
        remain(1) = .false.
        remain(2) = .true.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_z%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_z%mpi_comm, comm_1d_z%myrank, ierr)
        call MPI_Comm_size(comm_1d_z%mpi_comm, comm_1d_z%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_z%mpi_comm, 0, 1, comm_1d_z%west_rank, comm_1d_z%east_rank, ierr)

    end subroutine mpi_topology_create

    !>
    !> @brief       Destroy the communicator for cartesian topology.
    !>
    subroutine mpi_topology_destroy()

        implicit none
        integer(kind=4) :: ierr

        call MPI_Comm_free(mpi_world_cart, ierr)

    end subroutine mpi_topology_destroy

end module mpi_topology
