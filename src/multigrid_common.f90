module multigrid_common

    implicit none

    contains

    subroutine multigrid_common_print_subdomain_grid_info(sdm, mg_sdm, n_levels)

        use mpi_topology, only  : myrank
        use geometry, only      : subdomain

        implicit none

        type(subdomain), intent(in)     :: sdm
        type(subdomain), intent(in)     :: mg_sdm(:)
        integer(kind=4), intent(in)     :: n_levels
        integer(kind=4)                 :: l
        character(256)                  :: myfilename

        write(myfilename,'(a,i0.3)') './result/mg.',myrank
        open(myrank, file=myfilename, form='formatted')
        write(myrank,*) '========== 0 ========='
        write(myrank,'(3(a,i5))'), 'Mesh: Nx =', sdm%nx, ', Ny =',  sdm%ny,  ', Nz =', sdm%nz
        write(myrank,101) 'dxm:',myrank, sdm%dxm
        write(myrank,101) 'dxg:',myrank, sdm%dxg
        write(myrank,101) 'crx:',myrank, sdm%xg
        write(myrank,101) 'dym:',myrank, sdm%dym
        write(myrank,101) 'dyg:',myrank, sdm%dyg
        write(myrank,101) 'cry:',myrank, sdm%yg
        write(myrank,101) 'dzm:',myrank, sdm%dzm
        write(myrank,101) 'dzg:',myrank, sdm%dzg
        write(myrank,101) 'crz:',myrank, sdm%zg
    101 format (a4,i3,130f10.5)

        do l = 1, n_levels

            write(myrank,*) '==========',l,'========='
            write(myrank,'(3(a,i5))'), 'Mesh: Nx =', mg_sdm(l)%nx, ', Ny =',  mg_sdm(l)%ny,  ', Nz =', mg_sdm(l)%nz
            write(myrank,102) 'dxm:',myrank, mg_sdm(l)%dxm
            write(myrank,102) 'dxg:',myrank, mg_sdm(l)%dxg
            write(myrank,102) 'crx:',myrank, mg_sdm(l)%xg
            write(myrank,102) 'dym:',myrank, mg_sdm(l)%dym
            write(myrank,102) 'dyg:',myrank, mg_sdm(l)%dyg
            write(myrank,102) 'cry:',myrank, mg_sdm(l)%yg
            write(myrank,102) 'dzm:',myrank, mg_sdm(l)%dzm
            write(myrank,102) 'dzg:',myrank, mg_sdm(l)%dzg
            write(myrank,102) 'crz:',myrank, mg_sdm(l)%zg
        102 format (a4,i3,130f10.5)

        enddo

        close(myrank)

    end subroutine multigrid_common_print_subdomain_grid_info

    subroutine multigrid_common_print_poisson_matrix_info(mg_sdm, mg_a_poisson, lv)

        use geometry, only      : subdomain
        use matrix, only        : matrix_poisson
        use mpi_topology, only  : myrank

        implicit none

        type(subdomain), intent(in)         :: mg_sdm
        type(matrix_poisson), intent(in)    :: mg_a_poisson
        integer(kind=4), intent(in)         :: lv

        integer(kind=4)                     :: i, j, k
        character(256)                      :: myfilename

        write(myfilename,'(a,i0.3,a,i0.3)') './result/csr.', lv, '.',myrank
        open(myrank, file=myfilename, form='formatted')
            
        do k = 1, mg_sdm%nz
            do j = 1, mg_sdm%ny
                do i = 1, mg_sdm%nx
                    write(myrank,'(3(i5,x),7(e14.7,x))') i,j,k,mg_a_poisson%coeff(0,i,j,k),mg_a_poisson%coeff(1,i,j,k),mg_a_poisson%coeff(2,i,j,k) &
                    ,mg_a_poisson%coeff(3,i,j,k), mg_a_poisson%coeff(4,i,j,k), mg_a_poisson%coeff(5,i,j,k), mg_a_poisson%coeff(6,i,j,k)
                enddo
            enddo
        enddo
        close(myrank)

    end subroutine multigrid_common_print_poisson_matrix_info

    subroutine multigrid_common_print_restriction_info(val_c, val_f, vol_c, vol_f, level, nx_c, ny_c, nz_c, i_gl_c, j_gl_c, k_gl_c, &
                                                      i_stride_f, j_stride_f, k_stride_f, i_offset_f, j_offset_f, k_offset_f)

        use mpi_topology, only : myrank

        implicit none

        real(kind=8), intent(in)        :: val_c(0:,0:,0:)
        real(kind=8), intent(in)        :: val_f(0:,0:,0:)
        real(kind=8), intent(in)        :: vol_f(2,2,2), vol_c
        integer(kind=4), intent(in)     :: level
        integer(kind=4), intent(in)     :: i_gl_c, j_gl_c, k_gl_c
        integer(kind=4), intent(in)     :: i_stride_f, i_offset_f
        integer(kind=4), intent(in)     :: j_stride_f, j_offset_f
        integer(kind=4), intent(in)     :: k_stride_f, k_offset_f
        integer(kind=4), intent(in)     :: nx_c, ny_c, nz_c

        integer(kind=4)                 :: i, j, k
        integer(kind=4)                 :: i_c, j_c, k_c
        integer(kind=4)                 :: iz_f, jz_f, kz_f
        integer(kind=4)                 :: ip_f, jp_f, kp_f

        character(256)  :: myfilename

        write(myfilename,'(a,i0.3,a,i0.3)') './result/restriction_lv_',level,'.',myrank
        open(myrank, file=myfilename, form='formatted')
            
        write(myrank,'(4(a18,x))') 'val_f', 'vol_f', 'val_c', 'vol_c'
        do k = 1, nz_c
            k_c  = k + k_gl_c
            kp_f = k * k_stride_f
            kz_f = kp_f - k_offset_f
            do j = 1, ny_c
                j_c  = j + j_gl_c
                jp_f = j * j_stride_f
                jz_f = jp_f - j_offset_f
                do i = 1, nx_c
                    i_c  = i + i_gl_c
                    ip_f = i * i_stride_f
                    iz_f = ip_f - i_offset_f
                    write(myrank,'(4(e15.7,x))') val_f(iz_f, jz_f, kz_f), vol_f(2,2,2), val_c(i_c, j_c, k_c), vol_c
                    write(myrank,'(2(e15.7,x))') val_f(ip_f, jz_f, kz_f), vol_f(1,2,2) 
                    write(myrank,'(2(e15.7,x))') val_f(iz_f, jp_f, kz_f), vol_f(2,1,2) 
                    write(myrank,'(2(e15.7,x))') val_f(ip_f, jp_f, kz_f), vol_f(1,1,2) 
                    write(myrank,'(2(e15.7,x))') val_f(iz_f, jz_f, kp_f), vol_f(2,2,1) 
                    write(myrank,'(2(e15.7,x))') val_f(ip_f, jz_f, kp_f), vol_f(1,2,1) 
                    write(myrank,'(2(e15.7,x))') val_f(iz_f, jp_f, kp_f), vol_f(2,1,1)
                    write(myrank,'(2(e15.7,x))') val_f(ip_f, jp_f, kp_f), vol_f(1,1,1)
                enddo
            enddo
        enddo

        close(myrank)

        write(myfilename,'(a,i0.3,a,i0.3)') './result/restriction_coarse_mesh_lv_',level,'.',myrank
        open(myrank, file=myfilename, form='formatted')
            
        write(myrank,'(a18)') 'val_c'
        do k = 1, nz_c
            do j = 1, ny_c
                do i = 1,nx_c
                    write(myrank,'(3(i5,x),e15.7,x)') i, j, k, val_c(i,j,k)
                enddo
            enddo
        enddo

        close(myrank)

    end subroutine multigrid_common_print_restriction_info

    subroutine multigrid_common_print_prolongation_info(val_c, val_f, level, nx_f, ny_f, nz_f, i_gl_c, j_gl_c, k_gl_c, &
                                                        i_stride_f, j_stride_f, k_stride_f, i_offset_f, j_offset_f, k_offset_f)

        use mpi_topology, only : myrank

        implicit none

        real(kind=8), intent(in)        :: val_c(0:,0:,0:)
        real(kind=8), intent(in)        :: val_f(0:,0:,0:)
        integer(kind=4), intent(in)     :: level
        integer(kind=4), intent(in)     :: nx_f, ny_f, nz_f
        integer(kind=4), intent(in)     :: i_gl_c, j_gl_c, k_gl_c
        integer(kind=4), intent(in)     :: i_stride_f, i_offset_f
        integer(kind=4), intent(in)     :: j_stride_f, j_offset_f
        integer(kind=4), intent(in)     :: k_stride_f, k_offset_f

        integer(kind=4)                 :: i, j, k
        integer(kind=4)                 :: iz_c, jz_c, kz_c
        integer(kind=4)                 :: ip_c, jp_c, kp_c

        character(256)  :: myfilename

        write(myfilename,'(a,i0.3,a,i0.3)') './result/prolongation_lv_',level,'.',myrank
        open(myrank, file=myfilename, form='formatted')
            
        write(myrank,'(2(a18,x))') 'val_f', 'val_c'
        do k = 1, nz_f
            kz_c = int(k/k_stride_f) + k_gl_c
            kp_c = kz_c + k_offset_f
            do j = 1, ny_f
                jz_c = int(j/j_stride_f) + j_gl_c
                jp_c = jz_c + j_offset_f
            do i = 1, nx_f
                    iz_c = int(i/i_stride_f) + i_gl_c
                    ip_c = iz_c + i_offset_f
                    ! i = 1, iz_c = 0,1
                    ! i = 2, iz_c = 1,2
                    ! i = 3, iz_c = 1,2
                    ! i = 4, iz_c = 2,3
                    write(myrank,'(9(e15.7,x))') val_f(i,j,k), &
                                                val_c(iz_c,jz_c,kz_c), &
                                                val_c(ip_c,jz_c,kz_c), &
                                                val_c(iz_c,jp_c,kz_c), &
                                                val_c(ip_c,jp_c,kz_c), &
                                                val_c(iz_c,jz_c,kp_c), &
                                                val_c(ip_c,jz_c,kp_c), &
                                                val_c(iz_c,jp_c,kp_c), &
                                                val_c(ip_c,jp_c,kp_c)
                                                
                enddo
            enddo
        enddo

        close(myrank)

    end subroutine multigrid_common_print_prolongation_info

    subroutine multigrid_common_print_coarsest_level_solution(cyc, mg_sdm)

        use mpi_topology, only  : myrank
        use geometry, only      : subdomain

        implicit none

        integer(kind=4), intent(in)             :: cyc
        type(subdomain), intent(in)             :: mg_sdm
        integer(kind=4) :: i, j, k
        character(256)  :: myfilename

        write(myfilename,'(a,i0.3,a,i0.3)') './result/solution_lv_last_cycle_',cyc,'.',myrank
        open(myrank, file=myfilename, form='formatted')

        do k=1, mg_sdm%nz
            do j=1, mg_sdm%ny
                do i=1, mg_sdm%nx
                    write(myrank,101) mg_sdm%xg(i), mg_sdm%yg(j), mg_sdm%zg(k), &
                                      mg_sdm%r(i,j,k), mg_sdm%x(i,j,k), mg_sdm%b(i,j,k)
                enddo
            enddo
        enddo
        101   format(3(f9.5,x),3(e16.8))

        close(myrank)

    end subroutine multigrid_common_print_coarsest_level_solution

end module multigrid_common
