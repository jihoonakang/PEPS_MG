module mpi_topology

    implicit none

    private

    integer(kind=4), public :: nprocs, myrank
    integer(kind=4), public :: mpi_world_cart
    integer(kind=4), public :: np_dim(0:2)  !!!!x1:0, x2:1, x3:2
    logical, public         :: period(0:2) !!!!x1:0, x2:1, x3:2

    type, public :: cart_comm_1d
        integer(kind=4) :: myrank, nprocs
        integer(kind=4) :: west_rank, east_rank
        integer(kind=4) :: mpi_comm
    end type cart_comm_1d

    type(cart_comm_1d), public, target      :: comm_1d_x, comm_1d_y, comm_1d_z

    public  :: mpi_topology_create
    public  :: mpi_topology_destroy

    contains

    subroutine mpi_topology_create()
        
        use mpi

        implicit none

        logical         :: remain(0:2)
        integer(kind=4) :: ierr

        call MPI_Cart_create( MPI_COMM_WORLD    &!  input  | integer(kind=4)      | Input communicator (handle).
                            , 3                 &!  input  | integer(kind=4)      | Number of dimensions of Cartesian grid (integer(kind=4)).
                            , np_dim            &!  input  | integer(kind=4)(1:3) | integer(kind=4) array of size ndims specifying the number of processes in each dimension.
                            , period            &!  input  | logical(1:3) | Logical array of size ndims specifying whether the grid is periodic (true=1) or not (false=0) in each dimension.
                            , .false.           &!  input  | logical      | Ranking may be reordered (true=1) or not (false=0) (logical).
                            , mpi_world_cart    &! *output | integer(kind=4)      | Communicator with new Cartesian topology (handle).
                            , ierr              &!  output | integer(kind=4)      | Fortran only: Error status
                            )

        remain(0) = .true.
        remain(1) = .false.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_x%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x%mpi_comm, comm_1d_x%myrank, ierr)
        call MPI_Comm_size(comm_1d_x%mpi_comm, comm_1d_x%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_x%mpi_comm, 0, 1, comm_1d_x%west_rank, comm_1d_x%east_rank, ierr)
    
        remain(0) = .false.
        remain(1) = .true.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_y%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_y%mpi_comm, comm_1d_y%myrank, ierr)
        call MPI_Comm_size(comm_1d_y%mpi_comm, comm_1d_y%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_y%mpi_comm, 0, 1, comm_1d_y%west_rank, comm_1d_y%east_rank, ierr)

        remain(0) = .false.
        remain(1) = .false.
        remain(2) = .true.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_z%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_z%mpi_comm, comm_1d_z%myrank, ierr)
        call MPI_Comm_size(comm_1d_z%mpi_comm, comm_1d_z%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_z%mpi_comm, 0, 1, comm_1d_z%west_rank, comm_1d_z%east_rank, ierr)

    end subroutine mpi_topology_create

    subroutine mpi_topology_destroy()

        implicit none
        integer(kind=4) :: ierr

        call MPI_Comm_free(mpi_world_cart, ierr)

    end subroutine mpi_topology_destroy

end module mpi_topology
