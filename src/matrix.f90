!======================================================================================================================
!> @file        matrix.f90
!> @brief       This file contains a module that defines a matrix type.
!> @details     Heptadiagonal matrix supported for finite difference method using 7-stencil points
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
!> @brief       Module for matrix type for finite difference method with 7-stencil points.
!> @details     The module containes the matrix type, its creator and destroyer.
!>
module matrix

    implicit none

    !> @brief   Heptadiagonam matrix for finite difference method with 7-stencil points.
    !> @details It contains the degree of freedom (row size) and its seven coefficients for each coorinate.
    type, public        :: matrix_heptadiagonal
        integer(kind=4)                 :: dof              !< Degree of freedom, nx*ny*nz
        real(kind=8), allocatable       :: coeff(:,:,:,:)   !< Coefficient matrix. Seven elements for each coordinate
    end type matrix_heptadiagonal

    contains

    !>
    !> @brief       Create the heptadiagonal matrix for finite difference method with 7-stencil points
    !> @param       a_poisson   Matrix of matrix_heptadiagonal type
    !> @param       sdm         Subdomain
    !>
    subroutine matrix_heptadiagonal_create(a_poisson, sdm)

        use geometry, only : subdomain
        implicit none

        type(matrix_heptadiagonal), intent(inout)   :: a_poisson
        type(subdomain), intent(in)                 :: sdm

        integer(kind=4) :: i, j, k
        real(kind=8)    :: dxm2i, dym2i, dzm2i
        real(kind=8)    :: dxmp2i, dxmn2i, dymp2i, dymn2i, dzmp2i, dzmn2i

        ! Degree of freedom (row size) is equal to the number of grids in subdomain
        a_poisson%dof = sdm%nx * sdm%ny * sdm%nz

        ! Heptadiagonal matrix allocation
        ! First index 0 : (i, j, k)
        ! First index 1 : (i-1, j, k)
        ! First index 2 : (i+1, j, k)
        ! First index 3 : (i, j-1, k)
        ! First index 4 : (i, j+1, k)
        ! First index 5 : (i, j, k-1)
        ! First index 6 : (i, j, k+1)
        allocate( a_poisson%coeff(0:6, sdm%nx, sdm%ny, sdm%nz) )

        a_poisson%coeff = 0.0d0

        ! Build coefficient matrix
        do k=1, sdm%nz
            dzm2i = 1.0d0 / (sdm%dzm(k) * sdm%dzm(k))
            dzmp2i = 1.0d0 / (sdm%dzm(k) * sdm%dzg(k))
            dzmn2i = 1.0d0 / (sdm%dzm(k) * sdm%dzg(k+1))

            do j=1, sdm%ny
                dym2i = 1.0d0 / (sdm%dym(j) * sdm%dym(j))
                dymp2i = 1.0d0 / (sdm%dym(j) * sdm%dyg(j))
                dymn2i = 1.0d0 / (sdm%dym(j) * sdm%dyg(j+1))

                do i=1, sdm%nx
                    dxm2i = 1.0d0 / (sdm%dxm(i) * sdm%dxm(i))
                    dxmp2i = 1.0d0 / (sdm%dxm(i) * sdm%dxg(i))
                    dxmn2i = 1.0d0 / (sdm%dxm(i) * sdm%dxg(i+1))

                    a_poisson%coeff(0,i,j,k) = -(dxmp2i + dxmn2i + dymp2i + dymn2i + dzmp2i + dzmn2i)
                    a_poisson%coeff(1,i,j,k) = 1.0d0 * dxmp2i
                    a_poisson%coeff(2,i,j,k) = 1.0d0 * dxmn2i
                    a_poisson%coeff(3,i,j,k) = 1.0d0 * dymp2i
                    a_poisson%coeff(4,i,j,k) = 1.0d0 * dymn2i
                    a_poisson%coeff(5,i,j,k) = 1.0d0 * dzmp2i
                    a_poisson%coeff(6,i,j,k) = 1.0d0 * dzmn2i

                enddo
            enddo
        enddo

    end subroutine matrix_heptadiagonal_create

    !>
    !> @brief       Destroy the heptadiagonal matrix for finite difference method with 7-stencil points
    !> @param       a_poisson   Matrix of matrix_heptadiagonal type
    !>
    subroutine matrix_heptadiagonal_destroy(a_poisson)

        implicit none

        type(matrix_heptadiagonal), intent(inout)     :: a_poisson

        a_poisson%dof = 0
        deallocate( a_poisson%coeff )
        
    end subroutine matrix_heptadiagonal_destroy

end module matrix