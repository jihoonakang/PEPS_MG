module matrix

    implicit none

    type, public        :: matrix_poisson
        integer(kind=4)                 :: dof
        real(kind=8), allocatable       :: coeff(:,:,:,:)
    end type matrix_poisson

    contains

    subroutine matrix_poisson_create(a_poisson, sdm)

        use geometry, only : subdomain
        implicit none

        type(matrix_poisson), intent(inout)     :: a_poisson
        type(subdomain), intent(in)         :: sdm

        integer(kind=4) :: i, j, k
        real(kind=8)    :: dxm2i, dym2i, dzm2i
        real(kind=8)    :: dxmp2i, dxmn2i, dymp2i, dymn2i, dzmp2i, dzmn2i

        a_poisson%dof = sdm%nx * sdm%ny * sdm%nz
        allocate( a_poisson%coeff(0:6, sdm%nx, sdm%ny, sdm%nz) )

        a_poisson%coeff = 0.0d0

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

    end subroutine matrix_poisson_create

    subroutine matrix_poisson_destroy(a_poisson)

        implicit none

        type(matrix_poisson), intent(inout)     :: a_poisson

        a_poisson%dof = 0
        deallocate( a_poisson%coeff )
        
    end subroutine matrix_poisson_destroy

end module matrix