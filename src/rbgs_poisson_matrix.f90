module rbgs_poisson_matrix

    use mpi
    use timer
    use poisson_matrix_operator

    implicit none

    private

    public  ::  rbgs_solver_poisson_matrix
    public  ::  rbgs_iterator_poisson_matrix

    contains

    subroutine rbgs_solver_poisson_matrix(sol, a_poisson, rhs, dm, maxiteration, tolerance, omega, is_aggregated)

        use matrix, only        : matrix_poisson
        use geometry, only      : subdomain, geometry_halocell_update_selectively
        use mpi_topology, only  : myrank, comm_1d_x, comm_1d_y, comm_1d_z

        implicit none

        real(kind=8),   intent(inout)       :: sol(0:,0:,0:)
        type(matrix_poisson), intent(in)    :: a_poisson
        real(kind=8),   intent(in)          :: rhs(0:,0:,0:)
        type(subdomain), intent(in)         :: dm
        integer(kind=4), intent(in)         :: maxiteration
        real(kind=8), intent(in)            :: tolerance
        real(kind=8), intent(in)            :: omega
        logical, intent(in)                 :: is_aggregated(0:2)

        integer(kind=4)     :: i, j, k, iter, ista, offset
        real(kind=8)        :: rsd0tol=0.0d0, temp, rsd_norm = 0.0d0
        real(kind=8), allocatable :: rsd(:,:,:)
        real(kind=8)        :: t0

        t0 = MPI_Wtime()

        allocate( rsd(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )

        rsd(:,:,:)  = 0.0d0

        call mv_mul_poisson_matrix0(rsd, a_poisson, sol, dm, is_aggregated)
!$omp parallel do shared(rsd,rhs)
        do k = 1, dm%nz
            do j = 1, dm%ny
                do i = 1, dm%nx
                    rsd(i,j,k) = rhs(i,j,k) - rsd(i,j,k)
                enddo
            enddo
        enddo
        call vv_dot_3d_matrix0(rsd0tol, rsd, rsd, dm%nx, dm%ny, dm%nz, is_aggregated)
        rsd_norm = rsd0tol

        offset = mod(mod(dm%nx,2)*mod(comm_1d_x%myrank,2) &
                    +mod(dm%ny,2)*mod(comm_1d_y%myrank,2) &
                    +mod(dm%nz,2)*mod(comm_1d_z%myrank,2),2 )
        do iter=0, maxiteration-1

            if ((mod(iter,20).EQ.0).AND.(myrank.eq.0)) then
                print 101, sqrt(rsd_norm), sqrt(rsd_norm/rsd0tol), sqrt(rsd0tol), iter
                101   format('   [RBGS solver] mse: ',e14.6,x,', r_mse: ',e14.6,x,', rsd0: ',e14.6,' at ',i5,' iterations.')
            endif
            call timer_comm_stamp0
            call geometry_halocell_update_selectively(sol, dm, is_aggregated)
            call timer_comm_stamp(34)

!$omp parallel do shared(a_poisson, sol)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    ista = 1 + mod(j+k+offset,2)
                    do i = ista, dm%nx, 2
                        temp= a_poisson%coeff(1,i,j,k) * sol(i-1,j,k) &
                            + a_poisson%coeff(2,i,j,k) * sol(i+1,j,k) &
                            + a_poisson%coeff(3,i,j,k) * sol(i,j-1,k) &
                            + a_poisson%coeff(4,i,j,k) * sol(i,j+1,k) &
                            + a_poisson%coeff(5,i,j,k) * sol(i,j,k-1) &
                            + a_poisson%coeff(6,i,j,k) * sol(i,j,k+1)
                        sol(i,j,k) = omega * ( rhs(i,j,k) - temp ) / a_poisson%coeff(0,i,j,k) + (1.0d0-omega)*sol(i,j,k)
                    enddo
                enddo
            enddo
            call timer_comm_stamp0
            call geometry_halocell_update_selectively(sol, dm, is_aggregated)
            call timer_comm_stamp(34)

!$omp parallel do shared(a_poisson, sol)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    ista = 1 + mod(j+k+1+offset,2)
                    do i = ista, dm%nx, 2
                        temp= a_poisson%coeff(1,i,j,k) * sol(i-1,j,k) &
                            + a_poisson%coeff(2,i,j,k) * sol(i+1,j,k) &
                            + a_poisson%coeff(3,i,j,k) * sol(i,j-1,k) &
                            + a_poisson%coeff(4,i,j,k) * sol(i,j+1,k) &
                            + a_poisson%coeff(5,i,j,k) * sol(i,j,k-1) &
                            + a_poisson%coeff(6,i,j,k) * sol(i,j,k+1)
                        sol(i,j,k) = omega * ( rhs(i,j,k) - temp ) / a_poisson%coeff(0,i,j,k) + (1.0d0-omega)*sol(i,j,k)
                    enddo
                enddo
            enddo

            call mv_mul_poisson_matrix0(rsd, a_poisson, sol, dm, is_aggregated)
!$omp parallel do shared(rsd,rhs)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    do i = 1, dm%nx
                        rsd(i,j,k) = rhs(i,j,k) - rsd(i,j,k)
                    enddo
                enddo
            enddo

            call vv_dot_3d_matrix0(rsd_norm, rsd, rsd, dm%nx, dm%ny, dm%nz, is_aggregated)
            if(sqrt(rsd_norm/rsd0tol).le.tolerance) then
                exit
            endif
           
        enddo

        if(myrank.eq.0) print '(a,i5,a,e15.7,a,e15.7,a,e15.7)','   [RBGS solver] Solution obtained. Iter. = ',iter,', mse = ',sqrt(rsd_norm),', r_mse = ',sqrt(rsd_norm/rsd0tol), ', time = ', MPI_Wtime()-t0

        deallocate(rsd)

    end subroutine rbgs_solver_poisson_matrix

    subroutine rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, dm, maxiteration, omega, is_aggregated)

        use matrix, only        : matrix_poisson
        use geometry, only      : subdomain, geometry_halocell_update_selectively
        use mpi_topology, only  : myrank, comm_1d_x, comm_1d_y, comm_1d_z

        implicit none

        real(kind=8),   intent(inout)       :: sol(0:,0:,0:)
        type(matrix_poisson), intent(in)    :: a_poisson
        real(kind=8),   intent(in)          :: rhs(0:,0:,0:)
        type(subdomain), intent(in)         :: dm
        integer(kind=4),intent(in)          :: maxiteration
        real(kind=8), intent(in)            :: omega
        logical, intent(in)                 :: is_aggregated(0:2)

        integer(kind=4)     :: i, j, k, iter, ista, offset, is_same(0:2)
        real(kind=8)        :: temp

        if(is_aggregated(0) .eq. .true.) then
            is_same(0) = 0
        else
            is_same(0) = 1
        endif

        if(is_aggregated(1) .eq. .true.) then
            is_same(1) = 0
        else
            is_same(1) = 1
        endif

        if(is_aggregated(2) .eq. .true.) then
            is_same(2) = 0
        else
            is_same(2) = 1
        endif

        offset = mod(mod(dm%nx,2)*mod(comm_1d_x%myrank*is_same(0),2) &
                    +mod(dm%ny,2)*mod(comm_1d_y%myrank*is_same(1),2) &
                    +mod(dm%nz,2)*mod(comm_1d_z%myrank*is_same(2),2),2 )

        do iter=0, maxiteration-1

            call timer_comm_stamp0
            call geometry_halocell_update_selectively(sol, dm, is_aggregated)
            call timer_comm_stamp(33)

!$omp parallel do shared(a_poisson, sol)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    ista = 1 + mod(j+k+offset,2)
                    do i = ista, dm%nx, 2
                        temp= a_poisson%coeff(1,i,j,k) * sol(i-1,j,k) &
                            + a_poisson%coeff(2,i,j,k) * sol(i+1,j,k) &
                            + a_poisson%coeff(3,i,j,k) * sol(i,j-1,k) &
                            + a_poisson%coeff(4,i,j,k) * sol(i,j+1,k) &
                            + a_poisson%coeff(5,i,j,k) * sol(i,j,k-1) &
                            + a_poisson%coeff(6,i,j,k) * sol(i,j,k+1)
                        sol(i,j,k) = omega * ( rhs(i,j,k) - temp ) / a_poisson%coeff(0,i,j,k) + (1.0d0-omega)*sol(i,j,k)
                    enddo
                enddo
            enddo

            call timer_comm_stamp0
            call geometry_halocell_update_selectively(sol, dm, is_aggregated)
            call timer_comm_stamp(33)

!$omp parallel do shared(a_poisson, sol)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    ista = 1 + mod(j+k+1+offset,2)
                    do i = ista, dm%nx, 2
                        temp= a_poisson%coeff(1,i,j,k) * sol(i-1,j,k) &
                            + a_poisson%coeff(2,i,j,k) * sol(i+1,j,k) &
                            + a_poisson%coeff(3,i,j,k) * sol(i,j-1,k) &
                            + a_poisson%coeff(4,i,j,k) * sol(i,j+1,k) &
                            + a_poisson%coeff(5,i,j,k) * sol(i,j,k-1) &
                            + a_poisson%coeff(6,i,j,k) * sol(i,j,k+1)
                        sol(i,j,k) = omega * ( rhs(i,j,k) - temp ) / a_poisson%coeff(0,i,j,k) + (1.0d0-omega)*sol(i,j,k)
                    enddo
                enddo
            enddo

        enddo

        if(myrank.eq.0) print '(a,i4)','   [RBGS iterator] Iteration completed. Iter. = ',iter
        
    end subroutine rbgs_iterator_poisson_matrix

end module rbgs_poisson_matrix
