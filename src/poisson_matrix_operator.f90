module poisson_matrix_operator

    implicit none

    contains

    subroutine  vv_dot_3d_matrix(result, x, y, nx, ny, nz, is_serial)

        use mpi
        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z

        implicit none

        real(kind=8),   intent(out)     :: result
        real(kind=8),   intent(in)      :: x(0:,0:,0:), y(0:,0:,0:)
        integer(kind=4),intent(in)      :: nx, ny, nz
        logical, intent(in)             :: is_serial(0:2)

        real(kind=8)                    :: result_local, result_x, result_xy
        integer(kind=4)                 :: i, j, k, ierr

        result_local = 0.0d0

!$omp parallel do shared(x,y) reduction(+:result_local)
        do k=1, nz
            do j=1, ny
                do i=1, nx
                    result_local = result_local + x(i,j,k) * y(i,j,k)
                enddo
            enddo
        enddo

        if(all(is_serial)) then
            result = result_local
        else if(.not.any(is_serial)) then
            call MPI_Allreduce(result_local, result, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        else
            if(is_serial(0) .eq. .false.) then
                call MPI_Allreduce(result_local, result_x, 1, MPI_REAL8, MPI_SUM, comm_1d_x%mpi_comm, ierr)
            else
                result_x = result_local
            endif
            if(is_serial(1) .eq. .false.) then
                call MPI_Allreduce(result_x, result_xy, 1, MPI_REAL8, MPI_SUM, comm_1d_y%mpi_comm, ierr)
            else
                result_xy = result_x
            endif
            if(is_serial(2) .eq. .false.) then
                call MPI_Allreduce(result_xy, result, 1, MPI_REAL8, MPI_SUM, comm_1d_z%mpi_comm, ierr)
            else
                result = result_xy
            endif
        endif

    end subroutine  vv_dot_3d_matrix

    subroutine mv_mul_poisson_matrix(y, a_poisson, x, dm, is_serial)

        use mpi
        use geometry, only      : subdomain, geometry_halocell_update_selectively
        use matrix, only        : matrix_heptadiagonal

        implicit none

        real(kind=8), intent(out)               :: y(0:,0:,0:)
        type(matrix_heptadiagonal), intent(in)  :: a_poisson
        real(kind=8), intent(inout)             :: x(0:,0:,0:)
        type(subdomain), intent(in)             :: dm
        logical, intent(in)                     :: is_serial(0:2)

        integer(kind=4)     :: i, j, k, nx, ny, nz

        nx = dm%nx
        ny = dm%ny
        nz = dm%nz

        if(.not.all(is_serial)) then
            call geometry_halocell_update_selectively(x, dm, is_serial)
        endif

!$omp parallel do shared(y)
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1       ! initialize y
                    y(i,j,k) = 0.0d0
                enddo
            enddo
        enddo

!$omp parallel do shared(y, a_poisson, x)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    y(i,j,k) = a_poisson%coeff(0,i,j,k) * x(i,j,k) &
                            + a_poisson%coeff(1,i,j,k) * x(i-1,j,k) &
                            + a_poisson%coeff(2,i,j,k) * x(i+1,j,k) &
                            + a_poisson%coeff(3,i,j,k) * x(i,j-1,k) &
                            + a_poisson%coeff(4,i,j,k) * x(i,j+1,k) &
                            + a_poisson%coeff(5,i,j,k) * x(i,j,k-1) &
                            + a_poisson%coeff(6,i,j,k) * x(i,j,k+1)
                enddo
            enddo
        enddo

    end subroutine mv_mul_poisson_matrix

end module poisson_matrix_operator
