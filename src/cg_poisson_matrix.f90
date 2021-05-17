module cg_poisson_matrix

    use poisson_matrix_operator

    implicit none

    private

    integer(kind=4)     :: row_size

    public  :: cg_solver_poisson_matrix
    public  :: cg_iterator_poisson_matrix

    contains

    subroutine cg_solver_poisson_matrix(sol, a_poisson, rhs, dm, maxiteration, tolerance, is_aggregated)

        use matrix, only        : matrix_poisson
        use geometry, only      : subdomain
        use mpi_topology, only  : myrank

        implicit none

        real(kind=8),   intent(inout)       :: sol(0:,0:,0:)
        type(matrix_poisson),intent(in)     :: a_poisson
        real(kind=8),   intent(in)          :: rhs(0:,0:,0:)
        type(subdomain), intent(in)         :: dm
        integer(kind=4),intent(in)          :: maxiteration
        real(kind=8),   intent(in)          :: tolerance
        logical, intent(in)                 :: is_aggregated(0:2)

        integer(kind=4)     :: i, j, k, iter
        real(kind=8)        :: alpha=0.0d0, beta=0.0d0, temp1=0.0d0, temp2=0.0d0, rsd0tol=0.0d0
        real(kind=8), allocatable :: p(:,:,:), Ax(:,:,:), Ap(:,:,:), rsd(:,:,:)

        row_size = a_poisson%dof

        allocate(   p(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )
        allocate(  Ax(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )
        allocate(  Ap(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )
        allocate( rsd(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )

        p(:,:,:)    = 0.0d0
        ax(:,:,:)   = 0.0d0
        ap(:,:,:)   = 0.0d0
        rsd(:,:,:)  = 0.0d0
        call mv_mul_poisson_matrix(ax, a_poisson, sol, dm, is_aggregated)

    ! #ifdef USE_MKL
    !     call dcopy(row_size, rhs, 1, rsd, 1)
    !     call daxpy(row_size, -1.0d0, ax, 1, rsd, 1)
    !     call dcopy(row_size, rsd, 1, p, 1)
    ! #else
    !$omp parallel do shared(rsd,rhs,ax,p)
        do k = 1, dm%nz
            do j = 1, dm%ny
                do i = 1, dm%nx
                    rsd(i,j,k) = rhs(i,j,k) - ax(i,j,k)
                    p(i,j,k)   = rsd(i,j,k)
                enddo
            enddo
        enddo
    ! #endif

        call vv_dot_3d_matrix(rsd0tol, rsd, rsd, dm%nx, dm%ny, dm%nz, is_aggregated)
        if(myrank.eq.0) print '(a)', '[CG] Conjugate gradient is started.'
        if(myrank.eq.0) print '(a,e15.7)', '[CG] Initial residual : ', rsd0tol


        do iter=0, maxiteration-1

            if ((mod(iter,10).EQ.0).AND.(myrank.eq.0)) then
                print 101, sqrt(temp2), sqrt(temp2/rsd0tol), sqrt(rsd0tol), tolerance, iter
                101   format('[CG] mse: ',e14.6,x,', r_mse: ',e14.6,x,', rsd0: ',e14.6,' with a tolerance criteria of ',e12.4,' at ',i5,' iterations.')
            endif

            call vv_dot_3d_matrix(temp1, rsd, rsd, dm%nx, dm%ny, dm%nz, is_aggregated)

            call mv_mul_poisson_matrix(ap, a_poisson, p, dm, is_aggregated)
            call vv_dot_3d_matrix(temp2, Ap, p, dm%nx, dm%ny, dm%nz, is_aggregated)

            alpha=temp1/temp2

    ! #ifdef USE_MKL
    !         call daxpy(row_size, alpha, p, 1, sol, 1)
    !         call daxpy(row_size, -alpha, Ap, 1, rsd, 1)
    ! #else
    !$omp parallel do shared(sol,p,rsd,Ap) firstprivate(alpha)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    do i = 1, dm%nx
                        sol(i,j,k) = sol(i,j,k) + alpha * p(i,j,k)
                        rsd(i,j,k) = rsd(i,j,k) - alpha * Ap(i,j,k)
                    enddo
                enddo
            enddo
    ! #endif

            call vv_dot_3d_matrix(temp2, rsd, rsd, dm%nx, dm%ny, dm%nz, is_aggregated)

            if (sqrt(temp2/rsd0tol) .LT. tolerance) exit 

            beta = temp2/temp1

    ! #ifdef USE_MKL
    !         call dscal(row_size, beta, p, 1)
    !         call daxpy(row_size, 1.0d0, rsd, 1, p, 1)
    ! #else
    !$omp parallel do shared(p,rsd) firstprivate(beta)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    do i = 1, dm%nx
                        p(i,j,k)= rsd(i,j,k) + beta * p(i,j,k)
                    enddo
                enddo
            enddo
    ! #endif

        enddo
        if(myrank.eq.0) print '(a,i5,a,e15.7,a,e15.7)','[CG] Finished with total iteration = ',iter,', mse = ',sqrt(temp2),'r_mse = ',sqrt(temp2/rsd0tol)

        deallocate(p)
        deallocate(Ax)
        deallocate(Ap)
        deallocate(rsd)

    end subroutine cg_solver_poisson_matrix

    subroutine cg_iterator_poisson_matrix(sol, a_poisson, rhs, dm, maxiteration, is_aggregated)

        use matrix, only        : matrix_poisson
        use geometry, only      : subdomain
        use mpi_topology, only  : myrank

        implicit none

        real(kind=8),   intent(inout)       :: sol(0:,0:,0:)
        type(matrix_poisson), intent(in)    :: a_poisson
        real(kind=8),   intent(in)          :: rhs(0:,0:,0:)
        type(subdomain), intent(in)         :: dm
        integer(kind=4),intent(in)          :: maxiteration
        logical, intent(in)                 :: is_aggregated(0:2)

        integer(kind=4)     :: i, j, k, iter
        real(kind=8)        :: alpha=0.0d0, beta=0.0d0, temp1=0.0d0, temp2=0.0d0
        real(kind=8), allocatable :: p(:,:,:), Ax(:,:,:), Ap(:,:,:), rsd(:,:,:)

        row_size = a_poisson%dof

        allocate(   p(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )
        allocate(  Ax(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )
        allocate(  Ap(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )
        allocate( rsd(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )

        p(:,:,:)    = 0.0d0
        ax(:,:,:)   = 0.0d0
        ap(:,:,:)   = 0.0d0
        rsd(:,:,:)  = 0.0d0
        call mv_mul_poisson_matrix(ax, a_poisson, sol, dm, is_aggregated)

    ! #ifdef USE_MKL
    !     call dcopy(row_size, rhs, 1, rsd, 1)
    !     call daxpy(row_size, -1.0d0, ax, 1, rsd, 1)
    !     call dcopy(row_size, rsd, 1, p, 1)
    ! #else
    !$omp parallel do shared(rsd,rhs,ax,p)
        do k = 1, dm%nz
            do j = 1, dm%ny
                do i = 1, dm%nx
                    rsd(i,j,k) = rhs(i,j,k) - ax(i,j,k)
                    p(i,j,k)   = rsd(i,j,k)
                enddo
            enddo
        enddo
    ! #endif

        do iter=0, maxiteration-1

            call vv_dot_3d_matrix(temp1, rsd, rsd, dm%nx, dm%ny, dm%nz, is_aggregated)

            call mv_mul_poisson_matrix(ap, a_poisson, p, dm, is_aggregated)
            call vv_dot_3d_matrix(temp2, Ap, p, dm%nx, dm%ny, dm%nz, is_aggregated)

            alpha=temp1/temp2

    ! #ifdef USE_MKL
    !         call daxpy(row_size, alpha, p, 1, sol, 1)
    !         call daxpy(row_size, -alpha, Ap, 1, rsd, 1)
    ! #else
    !$omp parallel do shared(sol,p,rsd,Ap) firstprivate(alpha)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    do i = 1, dm%nx
                        sol(i,j,k) = sol(i,j,k) + alpha * p(i,j,k)
                        rsd(i,j,k) = rsd(i,j,k) - alpha * Ap(i,j,k)
                    enddo
                enddo
            enddo
    ! #endif

            call vv_dot_3d_matrix(temp2, rsd, rsd, dm%nx, dm%ny, dm%nz, is_aggregated)

            beta = temp2/temp1

    ! #ifdef USE_MKL
    !         call dscal(row_size, beta, p, 1)
    !         call daxpy(row_size, 1.0d0, rsd, 1, p, 1)
    ! #else
    !$omp parallel do shared(p,rsd) firstprivate(beta)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    do i = 1, dm%nx
                        p(i,j,k)= rsd(i,j,k) + beta * p(i,j,k)
                    enddo
                enddo
            enddo
    ! #endif

        enddo
        if(myrank.eq.0) print '(a,i5,a,e15.7)','   [CG] CG iterator completed with total iteration = ',iter,', mse = ',sqrt(temp2)

        deallocate(p)
        deallocate(Ax)
        deallocate(Ap)
        deallocate(rsd)

    end subroutine cg_iterator_poisson_matrix

end module cg_poisson_matrix
