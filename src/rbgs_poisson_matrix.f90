!======================================================================================================================
!> @file        rbgs_poisson_matrix.f90
!> @brief       This file contains a module for red-black Gauss-Seidel method with heptadiagonal poisson matrix
!> @details     This file contains a module for red-black Gauss-Seidel (RBGS) solver with convergence criteria and
!>              RBGS iterator with iteration number
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
!> @brief       Module for for red-black Gauss-Seidel method with heptadiagonal poisson matrix
!> @details     Red-black Gauss-Seidel (RBGS) solver with convergence criteria and RBGS iterator with iteration number
!>
module rbgs_poisson_matrix

    use mpi
    use poisson_matrix_operator

    implicit none

    private

    public  ::  rbgs_solver_poisson_matrix
    public  ::  rbgs_iterator_poisson_matrix

    contains

    !>
    !> @brief       Red-black Gauss-Seidel solver with convergence criteria
    !> @param       sol             Result solution having a shape of 3D matrix
    !> @param       a_poisson       Heptadiagonal poisson matrix
    !> @param       rhs             RHS vector having a shape of 3D matrix
    !> @param       dm              Subdomain
    !> @param       maxiteration    Maximum number of iterations
    !> @param       tolerance       Convergence criteria
    !> @param       omega           Relexation factor
    !> @param       is_aggregated   Boolean whether domain is aggregated (.true.) or not (.false.)
    !>
    subroutine rbgs_solver_poisson_matrix(sol, a_poisson, rhs, dm, maxiteration, tolerance, omega, is_aggregated)

        use matrix, only        : matrix_heptadiagonal
        use geometry, only      : subdomain, geometry_halocell_update_selectively
        use mpi_topology, only  : myrank, comm_1d_x, comm_1d_y, comm_1d_z

        implicit none

        real(kind=8),   intent(inout)           :: sol(0:,0:,0:)
        type(matrix_heptadiagonal), intent(in)  :: a_poisson
        real(kind=8),   intent(in)              :: rhs(0:,0:,0:)
        type(subdomain), intent(in)             :: dm
        integer(kind=4), intent(in)             :: maxiteration
        real(kind=8), intent(in)                :: tolerance
        real(kind=8), intent(in)                :: omega
        logical, intent(in)                     :: is_aggregated(0:2)

        ! Local temporary variables
        integer(kind=4)     :: i, j, k, iter, ista, offset, is_same(0:2)
        real(kind=8)        :: rsd0tol=0.0d0, temp, rsd_norm = 0.0d0
        real(kind=8), allocatable :: rsd(:,:,:)

        allocate( rsd(0:dm%nx+1, 0:dm%ny+1, 0:dm%nz+1) )

        rsd(:,:,:)  = 0.0d0

        call mv_mul_poisson_matrix(rsd, a_poisson, sol, dm, is_aggregated)
!$omp parallel do shared(rsd,rhs)
        do k = 1, dm%nz
            do j = 1, dm%ny
                do i = 1, dm%nx
                    rsd(i,j,k) = rhs(i,j,k) - rsd(i,j,k)
                enddo
            enddo
        enddo
        call vv_dot_3d_matrix(rsd0tol, rsd, rsd, dm%nx, dm%ny, dm%nz, is_aggregated)
        rsd_norm = rsd0tol

        ! Starting offset calculation : Offset is zero if domain is aggregated
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

        ! Main solver
        do iter=0, maxiteration-1

#ifdef MESSAGE_DETAIL
            if ((mod(iter,20).EQ.0).AND.(myrank.eq.0)) then
                print 101, sqrt(rsd_norm), sqrt(rsd_norm/rsd0tol), sqrt(rsd0tol), iter
                101   format('   [RBGS solver] mse: ',e14.6,x,', r_mse: ',e14.6,x,', rsd0: ',e14.6,' at ',i5,' iterations.')
            endif
#endif
            ! Update ghostcell in the direction without aggregation
            call geometry_halocell_update_selectively(sol, dm, is_aggregated)

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
            ! Update ghostcell in the direction without aggregation
            call geometry_halocell_update_selectively(sol, dm, is_aggregated)

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

            call mv_mul_poisson_matrix(rsd, a_poisson, sol, dm, is_aggregated)

            ! Residual calculation and check whether convergence criteria is satisfied.
!$omp parallel do shared(rsd,rhs)
            do k = 1, dm%nz
                do j = 1, dm%ny
                    do i = 1, dm%nx
                        rsd(i,j,k) = rhs(i,j,k) - rsd(i,j,k)
                    enddo
                enddo
            enddo

            call vv_dot_3d_matrix(rsd_norm, rsd, rsd, dm%nx, dm%ny, dm%nz, is_aggregated)
            if(sqrt(rsd_norm/rsd0tol).le.tolerance) then
                exit
            endif
           
        enddo

        if(myrank.eq.0) print '(a,i5,a,e15.7,a,e15.7,a,e15.7)','   [RBGS solver] Solution obtained. Iter. = ',iter,', mse = ',sqrt(rsd_norm),', r_mse = ',sqrt(rsd_norm/rsd0tol), ', time = ', MPI_Wtime()-t0

        deallocate(rsd)

    end subroutine rbgs_solver_poisson_matrix

    !>
    !> @brief       Red-black Gauss-Seidel solver with convergence criteria
    !> @param       sol             Result solution having a shape of 3D matrix
    !> @param       a_poisson       Heptadiagonal poisson matrix
    !> @param       rhs             RHS vector having a shape of 3D matrix
    !> @param       dm              Subdomain
    !> @param       maxiteration    Maximum number of iterations
    !> @param       omega           Relexation factor
    !> @param       is_aggregated   Boolean whether domain is aggregated (.true.) or not (.false.)
    !>
    subroutine rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, dm, maxiteration, omega, is_aggregated)

        use matrix, only        : matrix_heptadiagonal
        use geometry, only      : subdomain, geometry_halocell_update_selectively
        use mpi_topology, only  : myrank, comm_1d_x, comm_1d_y, comm_1d_z

        implicit none

        real(kind=8),   intent(inout)           :: sol(0:,0:,0:)
        type(matrix_heptadiagonal), intent(in)  :: a_poisson
        real(kind=8),   intent(in)              :: rhs(0:,0:,0:)
        type(subdomain), intent(in)             :: dm
        integer(kind=4),intent(in)              :: maxiteration
        real(kind=8), intent(in)                :: omega
        logical, intent(in)                     :: is_aggregated(0:2)

        ! Local temporary variables
        integer(kind=4)     :: i, j, k, iter, ista, offset, is_same(0:2)
        real(kind=8)        :: temp

        ! Starting offset calculation : Offset is zero if domain is aggregated
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

        ! Main solver : Residual calculation and convergence check are not necessary
        do iter=0, maxiteration-1

            call geometry_halocell_update_selectively(sol, dm, is_aggregated)

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

            call geometry_halocell_update_selectively(sol, dm, is_aggregated)

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

#ifdef MESSAGE_DETAIL
        if(myrank.eq.0) print '(a,i4)','   [RBGS iterator] Iteration completed. Iter. = ',iter
#endif
        
    end subroutine rbgs_iterator_poisson_matrix

end module rbgs_poisson_matrix
