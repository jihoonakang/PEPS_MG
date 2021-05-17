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

    type, public :: boundary_comm
        integer(kind=4) :: myrank, nprocs
        integer(kind=4) :: mpi_comm
    end type boundary_comm

    type(boundary_comm), public     :: comm_boundary

    public  :: mpi_topology_create
    public  :: mpi_topology_destroy
    public  :: mpi_boundary_create

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

    subroutine mpi_boundary_create()

        use mpi

        implicit none

        integer(kind=4) :: i, j, k
        integer(kind=4) :: ierr
        integer(kind=4) :: coords(0:2)
        integer(kind=4) :: group_boundary, group_cart
        integer(kind=4), allocatable    :: rank(:)
        integer(kind=4) :: curr_rank, n_rank

        call MPI_Comm_group(mpi_world_cart, group_cart, ierr)

        allocate(rank(product(np_dim)))

        rank = -1
        n_rank = 0
        do i = 0, np_dim(0)-1
            do j = 0, np_dim(1)-1
                do k = 0, np_dim(2)-1
                    coords(0) = i
                    coords(1) = j
                    coords(2) = k

                    if((i.eq.0 .or. i.eq.np_dim(0)-1).and.(period(0).eq..false.)) then
                        call MPI_Cart_rank(mpi_world_cart, coords, curr_rank, ierr)
                        n_rank = n_rank + 1
                        rank(n_rank) = curr_rank
                        ! print *, myrank, coords(0), coords(1), coords(2), n_rank, rank(n_rank)
                    else if((j.eq.0 .or. j.eq.np_dim(1)-1).and.(period(1).eq..false.)) then
                        call MPI_Cart_rank(mpi_world_cart, coords, curr_rank, ierr)
                        n_rank = n_rank + 1
                        rank(n_rank) = curr_rank
                    else if((k.eq.0 .or. k.eq.np_dim(2)-1).and.(period(2).eq..false.)) then
                        call MPI_Cart_rank(mpi_world_cart, coords, curr_rank, ierr)
                        n_rank = n_rank + 1
                        rank(n_rank) = curr_rank
                    endif
                enddo
            enddo
        enddo
        if(n_rank.ne.0) then
            call MPI_Group_incl(group_cart, n_rank, rank, group_boundary, ierr)
        endif
        
        call MPI_Comm_create(mpi_world_cart, group_boundary, comm_boundary%mpi_comm, ierr)

        comm_boundary%nprocs = -1
        comm_boundary%myrank = -1

        if(comm_boundary%mpi_comm.ne.MPI_COMM_NULL) then
            call MPI_Comm_size(comm_boundary%mpi_comm, comm_boundary%nprocs, ierr)
            call MPI_Comm_rank(comm_boundary%mpi_comm, comm_boundary%myrank, ierr)
        endif

        deallocate(rank)
        call MPI_Group_free(group_cart,ierr)
        call MPI_Group_free(group_boundary,ierr)

    end subroutine mpi_boundary_create

end module mpi_topology
