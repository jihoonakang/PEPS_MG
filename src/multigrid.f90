!======================================================================================================================
!> @file        multigrid.f90
!> @brief       This file contains a module for V-cycle multigrid method with CGS, CGA, and CGPSA
!> @detail      CGS : Coarse grid solver without grid aggregation.
!>              CGA : Coarse grid aggregation in a target level. All grid variables are aggregated into a singel core.
!>              CGPSA : Coarse grid partial semi-aggregation in multiple levels.
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
!> @brief       Module for solving Poisson equation with V-cycle multigrid(MG) method.
!> @detail      CGS, CGA and CGPSA are supported in this module.
!>              CGS : Coarse grid solver without grid aggregation.
!>              CGA : Coarse grid aggregation in a target level. All grid variables are aggregated into a singel core.
!>              CGPSA : Coarse grid partial semi-aggregation in multiple levels.!>

module multigrid

    use mpi
    use geometry
    use matrix
    use cg_poisson_matrix
    use rbgs_poisson_matrix
    use poisson_matrix_operator
    use multigrid_debug

    implicit none

    private

    integer(kind=4)         :: lv_gdm_coarsest_x    !< Possible coarsest level in x-direction (odd grid number)
    integer(kind=4)         :: lv_gdm_coarsest_y    !< Possible coarsest level in y-direction (odd grid number
    integer(kind=4)         :: lv_gdm_coarsest_z    !< Possible coarsest level in z-direction (odd grid number
    integer(kind=4)         :: lv_gdm_coarsest_max  !< Max(lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z)
    integer(kind=4)         :: n_levels             !< Number of MG level
    integer(kind=4)         :: n_vcycles            !< Number of V-cycle
    integer(kind=4)         :: n_smooting           !< Number of iteration in smoothing operation
    integer(kind=4)         :: lv_aggregation       !< Aggregation level in CGA
    integer(kind=4), target :: lv_aggregation_x     !< Aggregation level in x-direction in CGPSA 
    integer(kind=4), target :: lv_aggregation_y     !< Aggregation level in y-direction in CGPSA
    integer(kind=4), target :: lv_aggregation_z     !< Aggregation level in x-direction in CGPSA
    integer(kind=4)         :: lv_aggregation_max   !< Maximum number of aggregation levels in CGPSA, for error check
    integer(kind=4)         :: aggregation_type     !< 0: CGS, 1: CGA, 2: CGPSA

    type(subdomain), allocatable, target            :: mg_sdm(:)        !< MG subdomain in each MG level
    type(matrix_heptadiagonal), allocatable         :: mg_a_poisson(:)  !< Poisson matrix in each MG level

    public  :: multigrid_create
    public  :: multigrid_solve_vcycle
    public  :: multigrid_destroy

    contains

    !>
    !> @brief       Create MG subdomains and Poisson matrix
    !> @detail      MG subdomain depends on aggregation type. Different subroutines are called according to 
    !>              aggregation type.
    !> @param       sdm                 Subdomain
    !> @param       nlevel              Number of MG level
    !> @param       ncycle              Number of V-cycle
    !> @param       aggr_method         Aggregation method, 0: CGS, 1: CGA, 2: CGPSA
    !> @param       single_aggr_level   Aggregation level for CGA
    !> @param       adaptive_aggr_level Aggregation levels in x-, y-, and z-directions for CGPSA
    !>

    subroutine multigrid_create(sdm, nlevel, ncycle, aggr_method, single_aggr_level, adaptive_aggr_level)

        use mpi
        use mpi_topology, only : myrank

        implicit none

        type(subdomain), intent(in)     :: sdm
        integer(kind=4), intent(in)     :: nlevel
        integer(kind=4), intent(in)     :: ncycle
        integer(kind=4), intent(in)     :: aggr_method
        integer(kind=4), intent(in)     :: single_aggr_level
        integer(kind=4), intent(in)     :: adaptive_aggr_level(3)

        integer(kind=4)                 :: l
        integer(kind=4)                 :: ierr

        n_levels = nlevel
        n_vcycles = ncycle
        aggregation_type = aggr_method
        lv_aggregation = single_aggr_level
        lv_aggregation_x = adaptive_aggr_level(1)
        lv_aggregation_y = adaptive_aggr_level(2)
        lv_aggregation_z = adaptive_aggr_level(3)
        lv_aggregation_max = maxval(adaptive_aggr_level)

        if(n_levels.le.1) then
            if(myrank.eq.0) then
                print '(a,i2)',    '[Error] The number of level should be larger than 1. Current number: ', n_levels
            endif
            call MPI_Finalize(ierr)
            stop
        endif

        ! Alloocate as many MG subdomains as n_levels
        allocate(mg_sdm(n_levels))

        ! Create MG subdomains according to aggregation_type
        select case(aggregation_type)
        case(0)
            call multigrid_subdomain_create_CGS(sdm)
        case(1)
            call multigrid_subdomain_create_CGA(sdm)
        case(2)
            call multigrid_subdomain_create_CGPSA(sdm)
        case default
            if(myrank.eq.0) print '(a,i2)', '[Error] Aggregation method should be 0, 1, or 2.', aggregation_type
            call MPI_Finalize(ierr)
            stop
        end select

        ! Allocate grid variables, x, b, and r
        call multigrid_allocate_subdomain_variables()

        if(myrank.eq.0) print '(a,3l2)', '[MG] Aggregation info. of subdomain: ',sdm%is_aggregated

        do l = 1, n_levels
            if(myrank.eq.0) print '(a,i3,a,3l2)', '[MG] Aggregation info. of level ',l,': ',mg_sdm(l)%is_aggregated
        enddo

        if(myrank.eq.0) print '(a)','[MG] Multigrid geometry constructed.'

        ! Alloocate as many Poisson matrices as n_levels
        allocate(mg_a_poisson(n_levels))

        ! Create as many MG Poisson matrices as n_levels
        do l = 1, n_levels
            call matrix_heptadiagonal_create(mg_a_poisson(l), mg_sdm(l))
        enddo
        
        if(myrank.eq.0) print '(a)','[MG] Poisson matrix in multigrid constructed.'

        ! Print debuging info
#ifdef DEBUG_GRID
        call multigrid_debug_print_subdomain_grid_info(sdm, mg_sdm, n_levels)
#endif

        ! Print debuging info
#ifdef DEBUG_MATRIX
        do l = 1, n_levels
            call multigrid_debug_print_poisson_matrix_info(mg_sdm(l), mg_a_poisson(l), l)
        enddo
#endif

    end subroutine multigrid_create

    !>
    !> @brief       Create MG subdomains for CGS
    !> @detail      Define MG domain and assign grid dimensions such as mesh size, grid spacing 
    !>              and grid coordinates in each MG levels.
    !> @param       sdm                 Subdomain
    !>
    subroutine multigrid_subdomain_create_CGS(sdm)

        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z
        use mpi_topology, only : myrank

        implicit none

        type(subdomain), intent(in), target     :: sdm

        integer(kind=4)             :: l, ierr
        integer(kind=4)             :: nx, ny, nz
        !> @{ Pointer for mesh size in current MG level
        real(kind=8), pointer       :: dxm_f(:), dym_f(:), dzm_f(:)
        !> @}
        !> @{ Pointer for grid spacing in current MG level
        real(kind=8), pointer       :: dxg_f(:), dyg_f(:), dzg_f(:)
        !> @}
        !> @{ Pointer for grid coordinates in current MG level
        real(kind=8), pointer       :: xg_f(:),  yg_f(:),  zg_f(:)
        !> @}

        if(myrank.eq.0) print '(a)', '[MG] Grid coarsening without aggregation.'

        ! Error check
        if(lv_aggregation.ne.0) then
            if(myrank.eq.0) print '(a)', '[Error] Aggretation level should be 0 for CGS method.'
            call MPI_Finalize(ierr)
            stop
        endif

        nx = sdm%nx
        ny = sdm%ny
        nz = sdm%nz

        ! is_aggregated is .true. if the number of process is equal to 1 (same as serial)
        if(comm_1d_x%nprocs.eq.1) then
            mg_sdm(1:n_levels)%is_aggregated(0) = .true.
        else
            mg_sdm(1:n_levels)%is_aggregated(0) = .false.
        endif

        if(comm_1d_y%nprocs.eq.1) then
            mg_sdm(1:n_levels)%is_aggregated(1) = .true.
        else
            mg_sdm(1:n_levels)%is_aggregated(1) = .false.
        endif

        if(comm_1d_z%nprocs.eq.1) then
            mg_sdm(1:n_levels)%is_aggregated(2) = .true.
        else
            mg_sdm(1:n_levels)%is_aggregated(2) = .false.
        endif

        ! Define the number of grids in MG domains and coarsest level in each coordinate dimensions.
        do l = 1, n_levels
            ! Divid the domain size by two as the MG level increases.
            ! Current level is assigned as coarsest level when the number of grid is odd.
            if(myrank.eq.0) print '(a,i2)', '[MG] Grid coarsening at level ', l
            ! Define coarsest level in x-dimension. 
            if(mod(nx,2).EQ.0) then
                nx = nx/2   ! Divide the domain size by two
                lv_gdm_coarsest_x = l
                if(myrank.eq.0) print '(3(a,i6),a)', '[MG] X-grid coarsening at level ', l, ': ',nx,' grids reduced to ',nx/2,' grids.'
            else
                if(mg_sdm(l)%is_aggregated(0).eq..true.) then
                    ! Number of grid can be odd beyond coarsest level only if domain is aggregated. (serial)
                    if(myrank.eq.0) then
                        print '(a,i2,a)',  '[MG] No more X-grid coarsening from level ', l,'. Keeping X-grid number.'
                        print '(2(a,i6))', '[MG] The number of grid is ', nx, ' and number processes is ', comm_1d_x%nprocs
                    endif
                else
                    ! Error if number of grid is beyond coarsest level when domain is not aggregated. (not serial)
                    if(myrank.eq.0) then
                        print '(a,i2)',    '[Error] X-grid coarsening impossible for multiple processes from level ', l
                        print '(2(a,i6))', '[Error] The number of grid per process is ', nx, ' and number processes is ', comm_1d_x%nprocs
                    endif
                    call MPI_Finalize(ierr)
                    stop
                endif
            endif

            ! Define coarsest level in y-dimension. 
            if(mod(ny,2).EQ.0) then
                ny = ny/2   ! Divide the domain size by two
                lv_gdm_coarsest_y = l
                if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Y-grid coarsening at level ', l, ': ',ny,' grids reduced to ',ny/2,' grids.'
            else
                if(mg_sdm(l)%is_aggregated(1).eq..true.) then
                    ! Number of grid can be odd beyond coarsest level only if domain is aggregated. (serial)
                    if(myrank.eq.0) then
                        print '(a,i2,a)',  '[MG] No more Y-grid coarsening from level ', l,'. Keeping Y-grid number.'
                        print '(2(a,i6))', '[MG] The number of grid is ', ny, ' and number processes is ', comm_1d_y%nprocs
                    endif
                else
                    ! Error if number of grid is beyond coarsest level when domain is not aggregated. (not serial)
                    if(myrank.eq.0) then
                        print '(a,i2)',    '[Error] Y-grid coarsening impossible for multiple processes from level ', l
                        print '(2(a,i6))', '[Error] The number of grid per process is ', ny, ' and number processes is ', comm_1d_y%nprocs
                    endif
                    call MPI_Finalize(ierr)
                    stop
                endif
            endif

            ! Define coarsest level in z-dimension. 
            if(mod(nz,2).EQ.0) then
                nz = nz/2   ! Divide the domain size by two
                lv_gdm_coarsest_z = l
                if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Z-grid coarsening at level ', l, ': ',nz,' grids reduced to ',nz/2,' grids.'
            else
                if(mg_sdm(l)%is_aggregated(2).eq..true.) then
                    ! Number of grid can be odd beyond coarsest level only if domain is aggregated. (serial)
                    if(myrank.eq.0) then
                        print '(a,i2,a)',  '[MG] No more Z-grid coarsening from level ', l,'. Keeping Z-grid number.'
                        print '(2(a,i6))', '[MG] The number of grid is ', nz, ' and number processes is ', comm_1d_z%nprocs
                    endif
                else
                    ! Error if number of grid is beyond coarsest level when domain is not aggregated. (not serial)
                    if(myrank.eq.0) then
                        print '(a,i2)',    '[Error] Z-grid coarsening impossible for multiple processes from level ', l
                        print '(2(a,i6))', '[Error] The number of grid per process is ', nz, ' and number processes is ', comm_1d_z%nprocs
                    endif
                    call MPI_Finalize(ierr)
                    stop
                endif
            endif

            ! Assign the number of grids in the current MG domain.
            mg_sdm(l)%nx = nx
            mg_sdm(l)%ny = ny
            mg_sdm(l)%nz = nz

        enddo

        ! Assign the max coarsest level and error check
        lv_gdm_coarsest_max = max(lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z)
        if(myrank.eq.0) then 
            print '(a,4i3)', '[MG] Number of levels : ',n_levels
            print '(a,4i3)', '[MG] Final coarsest levels max, x, y, z directions : ',lv_gdm_coarsest_max, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z
        endif
        if(n_levels.gt.lv_gdm_coarsest_max) then
            if(myrank.eq.0) then 
                print '(a)', '[MG] Number of levels is greater than the coarsest level in global domain. '
                print '(a)', '[MG] It is not allowed and reduce the number of levels.'
            endif
            call MPI_Finalize(ierr)
            stop
        endif

        if(myrank.eq.0) print '(a)', '[MG] Generating grid dimension '

        ! Setup the mesh size, grid spacing and grid coordinates in x-direction.
        ! Begin with the partitioned subdomain as parent (finer) MG subdomain
        dxm_f => sdm%dxm
        dxg_f => sdm%dxg
        xg_f  => sdm%xg

        do l = 1, n_levels

            nx = mg_sdm(l)%nx
            allocate(mg_sdm(l)%dxm(0:nx+1))
            allocate(mg_sdm(l)%dxg(0:nx+1))
            allocate(mg_sdm(l)%xg(0:nx+1))

            ! Setup the current MG domain from parent (finer) MG domain.
            call multigrid_subdomain_CGS_make_grid(mg_sdm(l)%dxm, mg_sdm(l)%dxg, mg_sdm(l)%xg, nx, &
                                                                    dxm_f, dxg_f, xg_f, sdm%ox, &
                                                                    l, lv_gdm_coarsest_x, comm_1d_x,'x')
            nullify(dxm_f, dxg_f, xg_f)

            ! Change the current MG subdomain to the next parent (finer) MG subdomain
            dxm_f => mg_sdm(l)%dxm
            dxg_f => mg_sdm(l)%dxg
            xg_f  => mg_sdm(l)%xg

        enddo
        nullify(dxm_f, dxg_f, xg_f)

        ! Setup the mesh size, grid spacing and grid coordinates in y-direction.
        ! Begin with the partitioned subdomain as parent (finer) MG subdomain
        dym_f => sdm%dym
        dyg_f => sdm%dyg
        yg_f  => sdm%yg

        do l = 1, n_levels

            ny = mg_sdm(l)%ny
            allocate(mg_sdm(l)%dym(0:ny+1))
            allocate(mg_sdm(l)%dyg(0:ny+1))
            allocate(mg_sdm(l)%yg(0:ny+1))

            ! Setup the current MG domain from parent (finer) MG domain.
            call multigrid_subdomain_CGS_make_grid(mg_sdm(l)%dym, mg_sdm(l)%dyg, mg_sdm(l)%yg, ny, &
                                                                    dym_f, dyg_f, yg_f, sdm%oy, &
                                                                    l, lv_gdm_coarsest_y, comm_1d_y,'y')


            nullify(dym_f, dyg_f, yg_f)

            ! Change the current MG subdomain to the next parent (finer) MG subdomain
            dym_f => mg_sdm(l)%dym
            dyg_f => mg_sdm(l)%dyg
            yg_f  => mg_sdm(l)%yg

        enddo
        nullify(dym_f, dyg_f, yg_f)

        ! Setup the mesh size, grid spacing and grid coordinates in z-direction.
        ! Begin with the partitioned subdomain as parent (finer) MG subdomain
        dzm_f => sdm%dzm
        dzg_f => sdm%dzg
        zg_f  => sdm%zg

        do l = 1, n_levels
            
            nz = mg_sdm(l)%nz
            allocate(mg_sdm(l)%dzm(0:nz+1))
            allocate(mg_sdm(l)%dzg(0:nz+1))
            allocate(mg_sdm(l)%zg(0:nz+1))

            ! Setup the current MG domain from parent (finer) MG domain.
            call multigrid_subdomain_CGS_make_grid(mg_sdm(l)%dzm, mg_sdm(l)%dzg, mg_sdm(l)%zg, nz, &
                                                                    dzm_f, dzg_f, zg_f, sdm%oz, &
                                                                    l, lv_gdm_coarsest_z, comm_1d_z,'z')

            
            nullify(dzm_f, dzg_f, zg_f)
            
            ! Change the current MG subdomain to the next parent (finer) MG subdomain
            dzm_f => mg_sdm(l)%dzm
            dzg_f => mg_sdm(l)%dzg
            zg_f  => mg_sdm(l)%zg
        enddo
        nullify(dzm_f, dzg_f, zg_f)

    end subroutine multigrid_subdomain_create_CGS

    !>
    !> @brief       Make meshes and grids of the current MG subdomain form the parent (finer) MG subdomain
    !>              in the defined direction
    !> @detail      Mesh and grid creation proceeds until lv_coarsest using communication information and direction
    !> @param       dxm             Mesh size in current MG subdomain
    !> @param       dxg             Grid spacing in current MG subdomain
    !> @param       xg              Grid coordinate in current MG subdomain
    !> @param       nx              Number of grid in current MG subdomain
    !> @param       dxm_f           Mesh size in finer MG subdomain
    !> @param       dxg_f           Grid spacing in finer MG subdomain
    !> @param       xg_f            Grid coordinate in finer MG subdomain
    !> @param       ox              Starting point of subdomain
    !> @param       lv_cur          Level of current MG subdomain
    !> @param       lv_coarsest     Coarsest level in the target direction, dir
    !> @param       comm_1d         1D Subcommunicator in the target direction, dir
    !> @param       dir             Target direction
    !>
    subroutine multigrid_subdomain_CGS_make_grid(dxm, dxg, xg, nx, dxm_f, dxg_f, xg_f, &
                                                                ox, lv_cur, lv_coarsest, comm_1d, dir)

        use mpi_topology, only : cart_comm_1d, myrank

        implicit none

        real(kind=8), intent(inout)         :: dxm(0:)
        real(kind=8), intent(inout)         :: dxg(0:)
        real(kind=8), intent(inout)         :: xg(0:)
        integer(kind=4), intent(in)         :: nx
        real(kind=8), pointer, intent(in)   :: dxm_f(:)
        real(kind=8), pointer, intent(in)   :: dxg_f(:)
        real(kind=8), pointer, intent(in)   :: xg_f(:)
        real(kind=8), intent(in)            :: ox
        integer(kind=4), intent(in)         :: lv_cur, lv_coarsest
        type(cart_comm_1d), intent(in)      :: comm_1d
        character, intent(in)               :: dir

        integer(kind = 4)           :: i, ierr
        integer(kind=4)             :: request1(2), request2(2)

        if(myrank.eq.0) print *, '[MG] level : ',lv_cur,', size : ', size(dxm_f), size(dxm), ', dir :', dir

        dxm = 0.0d0
        dxg = 0.0d0
        xg  = 0.0d0

        if(lv_cur.le.lv_coarsest) then
            ! Grid coarsening if lv_cur <= lv_coarsest
            ! Mese size of current level by summing two neighboring mesh sizes in finer level
            do i = 1, nx
                dxm(i) = dxm_f(2*i-1) + dxm_f(2*i)
            enddo
            ! Update the mesh size of ghostcell
            if(comm_1d%west_rank.ne.MPI_PROC_NULL) then
                ! If not boundary at x,y,z=0
                call MPI_Isend(dxm(1), 1, MPI_REAL8, comm_1d%west_rank, 1, comm_1d%mpi_comm, request1(1), ierr)
                call MPI_Irecv(dxm(0), 1, MPI_REAL8, comm_1d%west_rank, 2, comm_1d%mpi_comm, request1(2), ierr)
            else
                ! If boundary at x,y,z=0
                dxm(0) = dxm_f(0)
            endif
            if(comm_1d%east_rank.ne.MPI_PROC_NULL) then
                ! If not boundary at x,y,z=L
                call MPI_Isend(dxm(nx+0), 1, MPI_REAL8, comm_1d%east_rank, 2, comm_1d%mpi_comm, request2(1), ierr)
                call MPI_Irecv(dxm(nx+1), 1, MPI_REAL8, comm_1d%east_rank, 1, comm_1d%mpi_comm, request2(2), ierr)
            else
                ! If boundary at x,y,z=L
                dxm(nx+1) = dxm_f(2*nx+1)
            endif

            if(comm_1d%west_rank.ne.MPI_PROC_NULL) then
                call MPI_Waitall(2, request1, MPI_STATUSES_IGNORE, ierr)
            endif
            if(comm_1d%east_rank.ne.MPI_PROC_NULL) then
                call MPI_Waitall(2, request2, MPI_STATUSES_IGNORE, ierr)
            endif

            ! Update grid spacing and coordinate using mesh size information
            xg(0) = ox - 0.5d0 * dxm(0)
            do i = 1, nx
                dxg(i) = 0.5d0 * (dxm(i) + dxm(i-1))
                xg(i) = xg(i-1) + dxg(i)
            enddo
            dxg(nx+1) = 0.5d0 * (dxm(nx+1) + dxm(nx))
            xg(nx+1) = xg(nx) + dxg(nx+1)
        else
            ! No grid coarsening if lv_cur > lv_coarsest
            dxm(0:nx+1) = dxm_f(0:nx+1)
            dxg(0:nx+1) = dxg_f(0:nx+1)
            xg(0:nx+1)  = xg_f(0:nx+1)
        endif

    end subroutine multigrid_subdomain_CGS_make_grid

    !>
    !> @brief       Create MG subdomains for CGA
    !> @detail      Define MG domain and assign grid dimensions such as mesh size, grid spacing 
    !>              and grid coordinates in each MG levels.
    !> @param       sdm                 Subdomain
    !>
    subroutine multigrid_subdomain_create_CGA(sdm)

        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z
        use mpi_topology, only : myrank

        implicit none

        type(subdomain), intent(in), target     :: sdm

        integer(kind=4)             :: l, ierr
        integer(kind=4)             :: nx, ny, nz
        !> @{ Number of grid in aggregation levels, nx * number of MPI processes
        integer(kind=4)             :: nx_lv_aggregation, ny_lv_aggregation, nz_lv_aggregation
        !> @}
        !> @{ Pointer for mesh size in current MG level
        real(kind=8), pointer       :: dxm_f(:), dym_f(:), dzm_f(:)
        !> @}
        !> @{ Pointer for grid spacing in current MG level
        real(kind=8), pointer       :: dxg_f(:), dyg_f(:), dzg_f(:)
        !> @}
        !> @{ Pointer for grid coordinates in current MG level
        real(kind=8), pointer       :: xg_f(:),  yg_f(:),  zg_f(:)
        !> @}

        if(myrank.eq.0) print '(a,i2)', '[MG] Grid coarsening with aggretation level = ', lv_aggregation

        ! Error check
        if(lv_aggregation.eq.0) then
            if(myrank.eq.0) print '(a,i2)', '[Error] Aggretation level should be larger than 0'
            call MPI_Finalize(ierr)
            stop
        endif

        nx = sdm%nx
        ny = sdm%ny
        nz = sdm%nz

        ! Define the number of grids in MG domains and coarsest level in each coordinate dimensions.
        do l = 1, n_levels
            if(myrank.eq.0) print '(a,i2)', '[MG] Grid coarsening at level ', l
            if(l.lt.lv_aggregation) then
                ! Before aggregation level : Divide the number of grids in parent level by two (# is even) or return error (# is odd)
                if(mod(nx,2).eq.1) then
                    ! Error if the number of grids is odd before aggregation level
                    if(myrank.eq.0) then
                        print '(2(a,i6))', '[Error] X-grid coarsening impossible at level = ',l,', less than aggregation level = ', lv_aggregation
                        print '(2(a,i6))', '[Error] The number of grid per process is ', nx, ' and number processes is ', comm_1d_x%nprocs
                    endif
                    call MPI_Finalize(ierr)
                    stop
                else
                    ! Divide the parent grid size by two before aggregation level
                    mg_sdm(l)%nx = nx/2
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] X-grid coarsening at level ', l, ': ',nx,' grids reduced to ',nx/2,' grids.'
                endif

                if(mod(ny,2).eq.1) then
                    ! Error if the number of grids is odd before aggregation level
                    if(myrank.eq.0) then
                        print '(2(a,i2))', '[Error] Y-grid coarsening impossible at level = ',l,', less than aggregation level = ', lv_aggregation
                        print '(2(a,i6))', '[Error] The number of grid per process is ', ny, ' and number processes is ', comm_1d_y%nprocs
                    endif
                    call MPI_Finalize(ierr)
                    stop
                else
                    ! Divide the parent grid size by two before aggregation level
                    mg_sdm(l)%ny = ny/2
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Y-grid coarsening at level ', l, ': ',ny,' grids reduced to ',ny/2,' grids.'
                endif

                if(mod(nz,2).eq.1) then
                    ! Error if the number of grids is odd before aggregation level
                    if(myrank.eq.0) then
                        print '(2(a,i2))', '[Error] z-grid coarsening impossible at level = ',l,', less than aggregation level = ', lv_aggregation
                        print '(2(a,i6))', '[Error] The number of grid per process is ', nz, ' and number processes is ', comm_1d_z%nprocs
                    endif
                    call MPI_Finalize(ierr)
                    stop
                else
                    ! Divide the parent grid size by two before aggregation level
                    mg_sdm(l)%nz = nz/2
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Z-grid coarsening at level ', l, ': ',nz,' grids reduced to ',nz/2,' grids.'
                endif
            else if(l.eq.lv_aggregation) then
                ! At aggregation level : Enlarge the grid size, multiplied by the number of MPI processes
                if(mod(nx,2).eq.1) then
                    ! Error if the parent grid size is equal to one: No more processing
                    if(myrank.eq.0) then
                        print '(2(a,i2))', '[Error] X-grid coarsening impossible at level = ',l,', equal to aggregation level = ', lv_aggregation
                        print '(2(a,i6))', '[Error] The number of grid per process is ', nx, ' and number processes is ', comm_1d_x%nprocs
                    endif
                    call MPI_Finalize(ierr)
                    stop
                else
                    ! Enlarge the grid size, multiplied by the number of MPI processes
                    nx_lv_aggregation = int(nx/2)
                    mg_sdm(l)%nx = nx/2 * comm_1d_x%nprocs
                    lv_gdm_coarsest_x = l
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] X-grid coarsening at level ', l, ': ',nx,' grids reduced to ',nx/2,' grids.'
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] X-grid aggregation at level ', l, ': ',nx/2,' grids aggregated to ',mg_sdm(l)%nx ,' grids.'
                endif

                if(mod(ny,2).eq.1) then
                    ! Error if the parent grid size is equal to one: No more processing
                    if(myrank.eq.0) then
                        print '(2(a,i2))', '[Error] Y-grid coarsening impossible at level = ',l,', equal to aggregation level = ', lv_aggregation
                        print '(2(a,i6))', '[Error] The number of grid per process is ', ny, ' and number processes is ', comm_1d_y%nprocs
                    endif
                    call MPI_Finalize(ierr)
                    stop
                else
                    ! Enlarge the grid size, multiplied by the number of MPI processes
                    ny_lv_aggregation = int(ny/2)
                    mg_sdm(l)%ny = ny/2 * comm_1d_y%nprocs
                    lv_gdm_coarsest_y = l
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Y-grid coarsening at level ', l, ': ',ny,' grids reduced to ',ny/2,' grids.'
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Y-grid aggregation at level ', l, ': ',ny/2,' grids aggregated to ',mg_sdm(l)%ny ,' grids.'    
                endif

                if(mod(nz,2).eq.1) then
                    ! Error if the parent grid size is equal to one: No more processing
                    if(myrank.eq.0) then
                        print '(2(a,i2))', '[Error] Z-grid coarsening impossible at level = ',l,', equal to aggregation level = ', lv_aggregation
                        print '(2(a,i6))', '[Error] The number of grid per process is ', nz, ' and number processes is ', comm_1d_z%nprocs
                    endif
                    call MPI_Finalize(ierr)
                    stop
                else
                    ! Enlarge the grid size, multiplied by the number of MPI processes
                    nz_lv_aggregation = int(nz/2)
                    mg_sdm(l)%nz = nz/2 * comm_1d_z%nprocs
                    lv_gdm_coarsest_z = l
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Z-grid coarsening at level ', l, ': ',nz,' grids reduced to ',nz/2,' grids.'
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Z-grid aggregation at level ', l, ': ',nz/2,' grids aggregated to ',mg_sdm(l)%nz ,' grids.'    
                endif
            else if(l.gt.lv_aggregation) then
                ! Beyond aggregation level: Divide the number of grids in parent level by two (# is even) or keep it (# is odd)
                if(mod(nx,2).EQ.1) then
                    ! Keep the grid size
                    if(myrank.eq.0) then
                        print '(a,i2,a)',  '[MG] No more X-grid coarsening from level ', l,'. Keeping X-grid number.'
                        print '(2(a,i6))', '[MG] The number of grid is ', nx, ' and number processes is ', comm_1d_x%nprocs
                    endif
                    mg_sdm(l)%nx = nx
                else
                    ! Divide the parent grid size by two
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] X-grid coarsening at level ', l, ': ',nx,' grids reduced to ',nx/2,' grids.'
                    mg_sdm(l)%nx = nx/2
                    lv_gdm_coarsest_x = l
                endif

                if(mod(ny,2).EQ.1) then
                    ! Keep the grid size
                    if(myrank.eq.0) then
                        print '(a,i2,a)',  '[MG] No more Y-grid coarsening from level ', l,'. Keeping Y-grid number.'
                        print '(2(a,i6))', '[MG] The number of grid is ', ny, ' and number processes is ', comm_1d_y%nprocs
                    endif
                    mg_sdm(l)%ny = ny
                else
                    ! Divide the parent grid size by two before aggregation level
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Y-grid coarsening at level ', l, ': ',ny,' grids reduced to ',ny/2,' grids.'
                    mg_sdm(l)%ny = ny/2
                    lv_gdm_coarsest_y = l
                endif

                if(mod(nz,2).EQ.1) then
                    ! Keep the grid size
                    if(myrank.eq.0) then
                        print '(a,i2,a)',  '[MG] No more Z-grid coarsening from level ', l,'. Keeping Z-grid number.'
                        print '(2(a,i6))', '[MG] The number of grid is ', nz, ' and number processes is ', comm_1d_z%nprocs
                    endif
                    mg_sdm(l)%nz = nz
                else
                    ! Divide the parent grid size by two before aggregation level
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Z-grid coarsening at level ', l, ': ',nz,' grids reduced to ',nz/2,' grids.'
                    mg_sdm(l)%nz = nz/2
                    lv_gdm_coarsest_z = l
                endif
            endif

            ! Assign the grid sizes of finer level as those of current level
            nx = mg_sdm(l)%nx
            ny = mg_sdm(l)%ny
            nz = mg_sdm(l)%nz

            ! If the number of grids is odd, no more grid coarsening
            if(mod(nx*ny*nz,2).EQ.1) then
                exit
            endif
        enddo

        do l = 1, n_levels
            if(l.lt.lv_aggregation) then
                ! is_aggregated is .true. if the number of process is equal to 1 (same as serial) before aggregation level
                if(comm_1d_x%nprocs.eq.1) then
                    mg_sdm(l)%is_aggregated(0) = .true.
                else
                    mg_sdm(l)%is_aggregated(0) = .false.
                endif
                if(comm_1d_y%nprocs.eq.1) then
                    mg_sdm(l)%is_aggregated(1) = .true.
                else
                    mg_sdm(l)%is_aggregated(1) = .false.
                endif
                if(comm_1d_z%nprocs.eq.1) then
                    mg_sdm(l)%is_aggregated(2) = .true.
                else
                    mg_sdm(l)%is_aggregated(2) = .false.
                endif
            else if(l.ge.lv_aggregation) then
                ! is_aggregated is .true. from aggregation level
                mg_sdm(l)%is_aggregated(0:2) = .true.
            endif
        enddo

        ! Assign the max coarsest level and error check
        lv_gdm_coarsest_max = max(lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z)
        if(myrank.eq.0) then 
            print '(a,4i3)', '[MG] Number of levels = ',n_levels
            print '(a,4i3)', '[MG] Final coarsest levels max, x, y, z directions= ',lv_gdm_coarsest_max, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z
        endif

        ! Check the number of levels and the max coarsest level
        if(n_levels.gt.lv_gdm_coarsest_max) then
            if(myrank.eq.0) then 
                print '(a)', '[MG] Number of levels is greater than the coarsest level in global domain. '
                print '(a)', '[MG] It is not allowed and change the number of levels.'
            endif
            call MPI_Finalize(ierr)
            stop
        endif

        if(myrank.eq.0) then 
            print '(a,4i3)', '[MG] Final gather level    = ',lv_aggregation
        endif

        if(lv_aggregation.gt.n_levels) then
            if(myrank.eq.0) print '(a)', '[MG] Gather level is greater than the coarsest level and no gather level exists'
        endif

        if(myrank.eq.0) print '(a)', '[MG] Generating grid dimension '

        ! Setup the mesh size, grid spacing and grid coordinates in x-direction.
        ! Begin with the partitioned subdomain as parent (finer) MG subdomain
        dxm_f => sdm%dxm
        dxg_f => sdm%dxg
        xg_f  => sdm%xg

        do l = 1, n_levels

            nx = mg_sdm(l)%nx
            allocate(mg_sdm(l)%dxm(0:nx+1))
            allocate(mg_sdm(l)%dxg(0:nx+1))
            allocate(mg_sdm(l)%xg(0:nx+1))

            ! Setup the current MG domain from parent (finer) MG domain.
            call multigrid_subdomain_aggregation_make_grid(mg_sdm(l)%dxm, mg_sdm(l)%dxg, mg_sdm(l)%xg, nx, &
                                                            dxm_f, dxg_f, xg_f, sdm%ox, &
                                                            l, lv_gdm_coarsest_x, lv_aggregation, &
                                                            nx_lv_aggregation, comm_1d_x,'x')

            nullify(dxm_f, dxg_f, xg_f)

            ! Change the current MG subdomain to the next parent (finer) MG subdomain
            dxm_f => mg_sdm(l)%dxm
            dxg_f => mg_sdm(l)%dxg
            xg_f  => mg_sdm(l)%xg
        enddo

        nullify(dxm_f, dxg_f, xg_f)

        ! Setup the mesh size, grid spacing and grid coordinates in y-direction.
        ! Begin with the partitioned subdomain as parent (finer) MG subdomain
        dym_f => sdm%dym
        dyg_f => sdm%dyg
        yg_f  => sdm%yg

        do l = 1, n_levels

            ny = mg_sdm(l)%ny
            allocate(mg_sdm(l)%dym(0:ny+1))
            allocate(mg_sdm(l)%dyg(0:ny+1))
            allocate(mg_sdm(l)%yg(0:ny+1))
            mg_sdm(l)%dym = 0.0d0
            mg_sdm(l)%dyg = 0.0d0
            mg_sdm(l)%yg = 0.0d0

            ! Setup the current MG domain from parent (finer) MG domain.
            call multigrid_subdomain_aggregation_make_grid(mg_sdm(l)%dym, mg_sdm(l)%dyg, mg_sdm(l)%yg, ny, &
                                                            dym_f, dyg_f, yg_f, sdm%oy, &
                                                            l, lv_gdm_coarsest_y, lv_aggregation, &
                                                            ny_lv_aggregation, comm_1d_y,'y')
            nullify(dym_f, dyg_f, yg_f)

            ! Change the current MG subdomain to the next parent (finer) MG subdomain
            dym_f => mg_sdm(l)%dym
            dyg_f => mg_sdm(l)%dyg
            yg_f  => mg_sdm(l)%yg
        enddo

        nullify(dym_f, dyg_f, yg_f)

        ! Setup the mesh size, grid spacing and grid coordinates in z-direction.
        ! Begin with the partitioned subdomain as parent (finer) MG subdomain
        dzm_f => sdm%dzm
        dzg_f => sdm%dzg
        zg_f  => sdm%zg

        do l = 1, n_levels

            nz = mg_sdm(l)%nz
            allocate(mg_sdm(l)%dzm(0:nz+1))
            allocate(mg_sdm(l)%dzg(0:nz+1))
            allocate(mg_sdm(l)%zg(0:nz+1))
            mg_sdm(l)%dzm = 0.0d0
            mg_sdm(l)%dzg = 0.0d0
            mg_sdm(l)%zg = 0.0d0

            ! Setup the current MG domain from parent (finer) MG domain.
            call multigrid_subdomain_aggregation_make_grid(mg_sdm(l)%dzm, mg_sdm(l)%dzg, mg_sdm(l)%zg, nz, &
                                                            dzm_f, dzg_f, zg_f, sdm%oz, &
                                                            l, lv_gdm_coarsest_z, lv_aggregation, &
                                                            nz_lv_aggregation, comm_1d_z,'z')

            nullify(dzm_f, dzg_f, zg_f)

            ! Change the current MG subdomain to the next parent (finer) MG subdomain
            dzm_f => mg_sdm(l)%dzm
            dzg_f => mg_sdm(l)%dzg
            zg_f  => mg_sdm(l)%zg

        enddo

        nullify(dzm_f, dzg_f, zg_f)

    end subroutine multigrid_subdomain_create_CGA

    !>
    !> @brief       Make meshes and grids of the current MG subdomain form the parent (finer) MG subdomain
    !>              in the defined direction
    !> @detail      Mesh and grid creation proceeds until lv_coarsest using communication information and direction
    !> @param       dxm                 Mesh size in current MG subdomain
    !> @param       dxg                 Grid spacing in current MG subdomain
    !> @param       xg                  Grid coordinate in current MG subdomain
    !> @param       nx                  Number of grid in current MG subdomain
    !> @param       dxm_f               Mesh size in finer MG subdomain
    !> @param       dxg_f               Grid spacing in finer MG subdomain
    !> @param       xg_f                Grid coordinate in finer MG subdomain
    !> @param       ox                  Starting point of subdomain
    !> @param       lv_cur              Level of current MG subdomain
    !> @param       lv_coarsest         Coarsest level in the target direction, dir
    !> @param       lv_aggregation      Aggregation level in the target direction
    !> @param       nx_lv_aggregation   Number of grids at aggregation level in the target direction
    !> @param       comm_1d             1D Subcommunicator in the target direction, dir
    !> @param       dir                 Target direction
    !>
    subroutine multigrid_subdomain_aggregation_make_grid(dxm, dxg, xg, nx, dxm_f, dxg_f, xg_f, &
                                                        ox, lv_cur, lv_coarsest, lv_aggregation, &
                                                        nx_lv_aggregation, comm_1d, dir)
        use mpi_topology, only : cart_comm_1d, myrank

        implicit none

        real(kind=8), intent(inout)         :: dxm(0:)
        real(kind=8), intent(inout)         :: dxg(0:)
        real(kind=8), intent(inout)         :: xg(0:)
        integer(kind=4), intent(in)         :: nx
        real(kind=8), pointer, intent(in)   :: dxm_f(:)
        real(kind=8), pointer, intent(in)   :: dxg_f(:)
        real(kind=8), pointer, intent(in)   :: xg_f(:)
        real(kind=8), intent(in)            :: ox
        integer(kind=4), intent(in)         :: lv_cur, lv_coarsest, lv_aggregation
        integer(kind=4), intent(in)         :: nx_lv_aggregation
        type(cart_comm_1d), intent(in)      :: comm_1d
        character, intent(in)               :: dir

        integer(kind = 4)           :: i, ierr
        integer(kind=4)             :: request1(2), request2(2)
        real(kind=8), allocatable   :: dxm_gl(:)

        if(lv_cur.eq.lv_aggregation) then
            if(myrank.eq.0) print '(a,i2,a,i5,i5,a,a)', '[MG] level(aggregation) : ',lv_cur,', size : ', size(dxm_f), size(dxm), ', dir :', dir
        else
            if(myrank.eq.0) print '(a,i15,a,i5,i5,a,a)', '[MG] level : ',lv_cur,', size : ', size(dxm_f), size(dxm), ', dir :', dir
        endif
        
        dxm = 0.0d0
        dxg = 0.0d0
        xg = 0.0d0

        if(lv_cur.le.lv_coarsest) then
            ! Grid coarsening if lv_cur <= lv_coarsest
            ! Mese size of current level by summing two neighboring mesh sizes in finer level
            if(lv_cur.lt.lv_aggregation) then
                ! Grid coarsening before aggregation level

                ! Mese size of current level by summing two neighboring mesh sizes in finer level
                do i = 1, nx
                    dxm(i) = dxm_f(2*i-1) + dxm_f(2*i)
                enddo

                ! Update the mesh size of ghostcell
                if(comm_1d%west_rank.ne.MPI_PROC_NULL) then
                    ! If not boundary at x,y,z=0
                    call MPI_Isend(dxm(1), 1, MPI_REAL8, comm_1d%west_rank, 1, comm_1d%mpi_comm, request1(1), ierr)
                    call MPI_Irecv(dxm(0), 1, MPI_REAL8, comm_1d%west_rank, 2, comm_1d%mpi_comm, request1(2), ierr)
                else
                    ! If boundary at x,y,z=0
                    dxm(0) = dxm_f(0)
                endif
                if(comm_1d%east_rank.ne.MPI_PROC_NULL) then
                    ! If not boundary at x,y,z=L
                    call MPI_Isend(dxm(nx+0), 1, MPI_REAL8, comm_1d%east_rank, 2, comm_1d%mpi_comm, request2(1), ierr)
                    call MPI_Irecv(dxm(nx+1), 1, MPI_REAL8, comm_1d%east_rank, 1, comm_1d%mpi_comm, request2(2), ierr)
                else
                    ! If boundary at x,y,z=L
                    dxm(nx+1) = dxm_f(2*nx+1)
                endif

                if(comm_1d%west_rank.ne.MPI_PROC_NULL) then
                    call MPI_Waitall(2, request1, MPI_STATUSES_IGNORE, ierr)
                endif
                if(comm_1d%east_rank.ne.MPI_PROC_NULL) then
                    call MPI_Waitall(2, request2, MPI_STATUSES_IGNORE, ierr)
                endif

                ! Update grid spacing and coordinate using mesh size information
                xg(0) = ox - 0.5d0 * dxm(0)
                do i = 1, nx
                    dxg(i) = 0.5d0 * (dxm(i) + dxm(i-1))
                    xg(i) = xg(i-1) + dxg(i)
                enddo
                dxg(nx+1) = 0.5d0 * (dxm(nx+1) + dxm(nx))
                xg(nx+1) = xg(nx) + dxg(nx+1)
            else if(lv_cur.eq.lv_aggregation) then
                ! Grid coarsening at aggregation level, grid size is enlarged.
                allocate(dxm_gl(nx_lv_aggregation))

                ! Mese size of current level by summing two neighboring mesh sizes in finer level
                do i = 1, nx_lv_aggregation
                    dxm_gl(i) = dxm_f(2*i-1) + dxm_f(2*i)
                enddo

                ! Aggregate the grid information
                call MPI_Allgather(dxm_gl(1), nx_lv_aggregation, MPI_REAL8, dxm(1), nx_lv_aggregation, MPI_REAL8, comm_1d%mpi_comm, ierr)

                ! Ghostcell information broadcasting
                dxm(0) = dxm_f(0)
                dxm(nx+1) = dxm_f(2*nx_lv_aggregation+1)
                call MPI_Bcast(dxm(0), 1, MPI_REAL8, 0, comm_1d%mpi_comm, ierr)
                call MPI_Bcast(dxm(nx+1), 1, MPI_REAL8, comm_1d%nprocs-1, comm_1d%mpi_comm, ierr)

                ! Update grid spacing and coordinate using mesh size information
                xg(0) = ox - 0.5d0 * dxm(0)
                call MPI_Bcast(xg(0), 1, MPI_REAL8, 0, comm_1d%mpi_comm, ierr)
                do i = 1, nx
                    dxg(i) = 0.5d0 * (dxm(i) + dxm(i-1))
                    xg(i) = xg(i-1) + dxg(i)
                enddo
                dxg(nx+1) = 0.5d0 * (dxm(nx+1) + dxm(nx))
                xg(nx+1) = xg(nx) + dxg(nx+1)
                deallocate(dxm_gl)
            else if(lv_cur.gt.lv_aggregation) then
                ! Grid coarsening beyond aggregation level, no communication is required on ghostcell
                ! Mese size of current level by summing two neighboring mesh sizes in finer level
                do i = 1, nx
                    dxm(i) = dxm_f(2*i-1) + dxm_f(2*i)
                enddo
                dxm(0) = dxm_f(0)
                dxm(nx+1) = dxm_f(2*nx+1)

                ! Update grid spacing and coordinate using mesh size information
                xg(0) = (xg_f(0) + 0.5d0 * dxm_f(0)) - 0.5d0 * dxm(0)
                do i = 1, nx
                    dxg(i) = 0.5d0 * (dxm(i) + dxm(i-1))
                    xg(i) = xg(i-1) + dxg(i)
                enddo
                dxg(nx+1) = 0.5d0 * (dxm(nx+1) + dxm(nx))
                xg(nx+1) = xg(nx) + dxg(nx+1)
            endif
        else
            ! No grid coarsening if lv_cur > lv_coarsest
            dxm(0:nx+1) = dxm_f(0:nx+1)
            dxg(0:nx+1) = dxg_f(0:nx+1)
            xg(0:nx+1)  = xg_f(0:nx+1)
        endif

    end subroutine multigrid_subdomain_aggregation_make_grid

    !>
    !> @brief       Create MG subdomains for CGPSA
    !> @detail      Define MG domain and assign grid dimensions such as mesh size, grid spacing 
    !>              and grid coordinates in each MG levels.
    !> @param       sdm                 Subdomain
    !>
    subroutine multigrid_subdomain_create_CGPSA(sdm)

        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z
        use mpi_topology, only : myrank

        implicit none

        type(subdomain), target, intent(in)     :: sdm

        integer(kind=4)             :: l, ierr
        integer(kind=4)             :: nx, ny, nz
        !> @{ Number of grid in aggregation levels, nx * number of MPI processes
        integer(kind=4)             :: nx_lv_aggregation, ny_lv_aggregation, nz_lv_aggregation
        !> @}
        !> @{ Pointer for mesh size in current MG level
        real(kind=8), pointer       :: dxm_f(:), dym_f(:), dzm_f(:)
        !> @}
        !> @{ Pointer for grid spacing in current MG level
        real(kind=8), pointer       :: dxg_f(:), dyg_f(:), dzg_f(:)
        !> @}
        !> @{ Pointer for grid coordinates in current MG level
        real(kind=8), pointer       :: xg_f(:),  yg_f(:),  zg_f(:)
        !> @}

        if(myrank.eq.0) print '(a)', '[MG] Grid coarsening with adaptive aggretation. lv_aggregation is neglected'

        nx = sdm%nx
        ny = sdm%ny
        nz = sdm%nz

        ! Define the number of grids in MG domains and coarsest level in each coordinate dimensions.
        do l = 1, n_levels
            if(myrank.eq.0) print '(a,i2)', '[MG] Grid coarsening at level ', l

            if(comm_1d_x%nprocs.eq.1) then
                ! No aggregation if the MPI procress is equal to 1 in x-direction
                lv_aggregation_x = 0
                if(myrank.eq.0) print '(a,i2)', '[MG] Aggretaion level in x-direction is set to zero : Npx = 1.'
            endif

            if(mod(nx,2).eq.1) then
                ! Keep the grid size: no more grid coarsening
                mg_sdm(l)%nx = nx
                if(myrank.eq.0) then
                    print '(a,i2,a)',  '[MG] No more x-grid coarsening at level ', l,' in serial. Keeping x-grid number'
                    print '(a,i2)',    '[MG] Aggregation level in x-direction is ', lv_aggregation_x
                    print '(2(a,i6))', '[MG] The number of grid is ', mg_sdm(l)%nx , ' and number processes is ', comm_1d_x%nprocs
                endif
                if(lv_aggregation_x.ge.l) then
                    ! Error if the number of grids is odd before aggregation level
                    print '(2(a,i2))', '[Error] X-grid coarsening impossible at level = ',l,', <= X-aggregation level X= ', lv_aggregation_x
                    print '(2(a,i6))', '[Error] The number of grid per process is ', nx, ' and number processes is ', comm_1d_x%nprocs
                    call MPI_Abort(MPI_COMM_WORLD,301,ierr)
                endif
            else
                if(l.eq.lv_aggregation_x) then
                    ! At aggregation level : Enlarge the grid size, multiplied by the number of MPI processes, and coarsen the grids
                    nx_lv_aggregation = int(nx/2)
                    mg_sdm(l)%nx = nx/2 * comm_1d_x%nprocs
                    lv_gdm_coarsest_x = l
                    if(myrank.eq.0) then
                        print '(a,i2)',    '[MG] Aggregation level in x-direction is ', lv_aggregation_x
                        print '(2(a,i6))', '[MG] The number of grid is ', mg_sdm(l)%nx, ' and number processes is ', comm_1d_x%nprocs
                    endif
                else
                    ! Before and beyond aggregation level : Divide the number of grids in parent level by two (# is even) 
                    mg_sdm(l)%nx = nx/2
                    lv_gdm_coarsest_x = l
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] x-grid coarsening at level ', l, ': ',nx ,' grids reduced to ',mg_sdm(l)%nx ,' grids.'
                endif
            endif

            if(comm_1d_y%nprocs.eq.1) then
                ! No aggregation if the MPI procress is equal to 1 in y-direction
                lv_aggregation_y = 0
                if(myrank.eq.0) print '(a,i2)', '[MG] Aggretaion level in y-direction is set to zero : Npy = 1.'
            endif

            if(mod(ny,2).eq.1) then
                ! Keep the grid size: no more grid coarsening
                mg_sdm(l)%ny = ny
                if(myrank.eq.0) then
                    print '(a,i2,a)',  '[MG] No more y-grid coarsening at level ', l,' in serial. Keeping y-grid number'
                    print '(a,i2)',    '[MG] Aggregation level in y-direction is ', lv_aggregation_y
                    print '(2(a,i6))', '[MG] The number of grid is ', mg_sdm(l)%ny , ' and number processes is ', comm_1d_y%nprocs
                endif
                if(lv_aggregation_y.ge.l) then
                    ! Error if the number of grids is odd before aggregation level
                    print '(2(a,i2))', '[Error] Y-grid coarsening impossible at level = ',l,', <= Y-aggregation level X= ', lv_aggregation_y
                    print '(2(a,i6))', '[Error] The number of grid per process is ', ny, ' and number processes is ', comm_1d_y%nprocs
                    call MPI_Abort(MPI_COMM_WORLD,301,ierr)
                endif
            else
                if(l.eq.lv_aggregation_y) then
                    ! At aggregation level : Enlarge the grid size, multiplied by the number of MPI processes, and coarsen the grids
                    ny_lv_aggregation = int(ny/2)
                    mg_sdm(l)%ny = ny/2 * comm_1d_y%nprocs
                    lv_gdm_coarsest_y = l
                    if(myrank.eq.0) then
                        print '(a,i2)',    '[MG] Aggregation level in y-direction is ', lv_aggregation_y
                        print '(2(a,i6))', '[MG] The number of grid is ', mg_sdm(l)%ny, ' and number processes is ', comm_1d_y%nprocs
                    endif
                else
                    ! Before and beyond aggregation level : Divide the number of grids in parent level by two (# is even) 
                    mg_sdm(l)%ny = ny/2
                    lv_gdm_coarsest_y = l
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] y-grid coarsening at level ', l, ': ',ny ,' grids reduced to ',mg_sdm(l)%ny ,' grids.'
                endif
            endif

            if(comm_1d_z%nprocs.eq.1) then
                ! No aggregation if the MPI procress is equal to 1 in z-direction
                lv_aggregation_z = 0
                if(myrank.eq.0) print '(a,i2)', '[MG] Aggretaion level in z-direction is set to zero : Npz = 1.'
            endif

            if(mod(nz,2).eq.1) then
                ! Keep the grid size: no more grid coarsening
                mg_sdm(l)%nz = nz
                if(myrank.eq.0) then
                    print '(a,i2,a)',  '[MG] No more Z-grid coarsening at level ', l,' in serial. Keeping z-grid number'
                    print '(a,i2)',    '[MG] Aggregation level in z-direction is ', lv_aggregation_z
                    print '(2(a,i6))', '[MG] The number of grid is ', mg_sdm(l)%nz , ' and number processes is ', comm_1d_z%nprocs
                endif
                if(lv_aggregation_z.ge.l) then
                    ! Error if the number of grids is odd before aggregation level
                    print '(2(a,i2))', '[Error] Z-grid coarsening impossible at level = ',l,', <= Y-aggregation level X= ', lv_aggregation_z
                    print '(2(a,i6))', '[Error] The number of grid per process is ', nz, ' and number processes is ', comm_1d_z%nprocs
                    call MPI_Abort(MPI_COMM_WORLD,301,ierr)
                endif
            else
                if(l.eq.lv_aggregation_z) then
                    ! At aggregation level : Enlarge the grid size, multiplied by the number of MPI processes, and coarsen the grids
                    nz_lv_aggregation = int(nz/2)
                    mg_sdm(l)%nz = nz/2 * comm_1d_z%nprocs
                    lv_gdm_coarsest_z = l
                    if(myrank.eq.0) then
                        print '(a,i2)',    '[MG] Aggregation level in z-direction is ', lv_aggregation_z
                        print '(2(a,i6))', '[MG] The number of grid is ', mg_sdm(l)%nz, ' and number processes is ', comm_1d_z%nprocs
                    endif
                else
                    ! Before and beyond aggregation level : Divide the number of grids in parent level by two (# is even) 
                    mg_sdm(l)%nz = nz/2
                    lv_gdm_coarsest_z = l
                    if(myrank.eq.0) print '(3(a,i6),a)', '[MG] Z-grid coarsening at level ', l, ': ',nz ,' grids reduced to ',mg_sdm(l)%nz ,' grids.'
                endif
            endif

            ! Assign the grid sizes of finer level as those of current level
            nx = mg_sdm(l)%nx
            ny = mg_sdm(l)%ny
            nz = mg_sdm(l)%nz

            ! If the number of grids is odd, no more grid coarsening
            if(mod(nx*ny*nz,2).EQ.1) then
                exit
            endif
        enddo

        do l = 1, n_levels
            ! is_aggregated is .true. from aggregation level
            if(l.lt.lv_aggregation_x) then
                mg_sdm(l)%is_aggregated(0) = .false.
            else if(l.ge.lv_aggregation_x) then
                mg_sdm(l)%is_aggregated(0) = .true.
            endif
            if(l.lt.lv_aggregation_y) then
                mg_sdm(l)%is_aggregated(1) = .false.
            else if(l.ge.lv_aggregation_y) then
                mg_sdm(l)%is_aggregated(1) = .true.
            endif
            if(l.lt.lv_aggregation_z) then
                mg_sdm(l)%is_aggregated(2) = .false.
            else if(l.ge.lv_aggregation_z) then
                mg_sdm(l)%is_aggregated(2) = .true.
            endif
        enddo

        ! Assign the max coarsest level and error check
        lv_gdm_coarsest_max = max(lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z)
        lv_aggregation_max = max(lv_aggregation_x, lv_aggregation_y, lv_aggregation_z)
        if(myrank.eq.0) then 
            print '(a,4i3)', '[MG] Number of levels = ',n_levels
            print '(a,4i3)', '[MG] Final aggregation levels max x, y, z directions= ',lv_aggregation_max, lv_aggregation_x, lv_aggregation_y, lv_aggregation_z
            print '(a,4i3)', '[MG] Final coarsest levels max, x, y, z directions= ',lv_gdm_coarsest_max, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z
        endif

        ! Check the number of levels and the max coarsest level
        if(n_levels.ne.lv_gdm_coarsest_max) then
            if(myrank.eq.0) then 
                print '(a)', '[MG] Number of levels should be equal to the coarsest level in global domain. '
                print '(a)', '[MG] It is not allowed and change the number of levels.'
            endif
            call MPI_Finalize(ierr)
            stop
        endif

        ! Setup the mesh size, grid spacing and grid coordinates in x-direction.
        ! Begin with the partitioned subdomain as parent (finer) MG subdomain
        dxm_f => sdm%dxm
        dxg_f => sdm%dxg
        xg_f  => sdm%xg

        do l = 1, n_levels

            nx = mg_sdm(l)%nx
            allocate(mg_sdm(l)%dxm(0:nx+1))
            allocate(mg_sdm(l)%dxg(0:nx+1))
            allocate(mg_sdm(l)%xg(0:nx+1))

            ! Setup the current MG domain from parent (finer) MG domain.
            call multigrid_subdomain_aggregation_make_grid(mg_sdm(l)%dxm, mg_sdm(l)%dxg, mg_sdm(l)%xg, nx, &
                                                                    dxm_f, dxg_f, xg_f, sdm%ox, &
                                                                    l, lv_gdm_coarsest_x, lv_aggregation_x, &
                                                                    nx_lv_aggregation, comm_1d_x,'x')

            nullify(dxm_f, dxg_f, xg_f)

            ! Change the current MG subdomain to the next parent (finer) MG subdomain
            dxm_f => mg_sdm(l)%dxm
            dxg_f => mg_sdm(l)%dxg
            xg_f  => mg_sdm(l)%xg
        enddo

        nullify(dxm_f, dxg_f, xg_f)

        ! Setup the mesh size, grid spacing and grid coordinates in y-direction.
        ! Begin with the partitioned subdomain as parent (finer) MG subdomain
        dym_f => sdm%dym
        dyg_f => sdm%dyg
        yg_f  => sdm%yg

        do l = 1, n_levels

            ny = mg_sdm(l)%ny
            allocate(mg_sdm(l)%dym(0:ny+1))
            allocate(mg_sdm(l)%dyg(0:ny+1))
            allocate(mg_sdm(l)%yg(0:ny+1))
            mg_sdm(l)%dym = 0.0d0
            mg_sdm(l)%dyg = 0.0d0
            mg_sdm(l)%yg = 0.0d0

            ! Setup the current MG domain from parent (finer) MG domain.
            call multigrid_subdomain_aggregation_make_grid(mg_sdm(l)%dym, mg_sdm(l)%dyg, mg_sdm(l)%yg, ny, &
                                                            dym_f, dyg_f, yg_f, sdm%oy, &
                                                            l, lv_gdm_coarsest_y, lv_aggregation_y, &
                                                            ny_lv_aggregation, comm_1d_y,'y')
            nullify(dym_f, dyg_f, yg_f)

            ! Change the current MG subdomain to the next parent (finer) MG subdomain
            dym_f => mg_sdm(l)%dym
            dyg_f => mg_sdm(l)%dyg
            yg_f  => mg_sdm(l)%yg
        enddo

        nullify(dym_f, dyg_f, yg_f)

        ! Setup the mesh size, grid spacing and grid coordinates in z-direction.
        ! Begin with the partitioned subdomain as parent (finer) MG subdomain
        dzm_f => sdm%dzm
        dzg_f => sdm%dzg
        zg_f  => sdm%zg

        do l = 1, n_levels

            nz = mg_sdm(l)%nz
            allocate(mg_sdm(l)%dzm(0:nz+1))
            allocate(mg_sdm(l)%dzg(0:nz+1))
            allocate(mg_sdm(l)%zg(0:nz+1))
            mg_sdm(l)%dzm = 0.0d0
            mg_sdm(l)%dzg = 0.0d0
            mg_sdm(l)%zg = 0.0d0

            ! Setup the current MG domain from parent (finer) MG domain.
            call multigrid_subdomain_aggregation_make_grid(mg_sdm(l)%dzm, mg_sdm(l)%dzg, mg_sdm(l)%zg, nz, &
                                                            dzm_f, dzg_f, zg_f, sdm%oz, &
                                                            l, lv_gdm_coarsest_z, lv_aggregation_z, &
                                                            nz_lv_aggregation, comm_1d_z,'z')

            nullify(dzm_f, dzg_f, zg_f)

            ! Change the current MG subdomain to the next parent (finer) MG subdomain
            dzm_f => mg_sdm(l)%dzm
            dzg_f => mg_sdm(l)%dzg
            zg_f  => mg_sdm(l)%zg

        enddo

        nullify(dzm_f, dzg_f, zg_f)

    end subroutine multigrid_subdomain_create_CGPSA

    !>
    !> @brief       Allocate grid variables, x, b, and r, and derived datatypes in all MG subdomains
    !>
    subroutine multigrid_allocate_subdomain_variables

        implicit none

        integer(kind=4)     :: l

        do l = 1, n_levels
            allocate(mg_sdm(l)%b(0:mg_sdm(l)%nx+1, 0:mg_sdm(l)%ny+1, 0:mg_sdm(l)%nz+1 ))
            allocate(mg_sdm(l)%r(0:mg_sdm(l)%nx+1, 0:mg_sdm(l)%ny+1, 0:mg_sdm(l)%nz+1 ))
            allocate(mg_sdm(l)%x(0:mg_sdm(l)%nx+1, 0:mg_sdm(l)%ny+1, 0:mg_sdm(l)%nz+1 ))

            mg_sdm(l)%b(:,:,:) = 0.0d0
            mg_sdm(l)%r(:,:,:) = 0.0d0
            mg_sdm(l)%x(:,:,:) = 0.0d0

        enddo

        do l = 1, n_levels
            call geometry_subdomain_ddt_create(mg_sdm(l))
        enddo

    end subroutine multigrid_allocate_subdomain_variables


    !>
    !> @brief       Destroy all variables in all MG subdomains
    !>
    subroutine multigrid_destroy

        implicit none

        integer(kind=4)     :: l

        do l = 1, n_levels
            call matrix_heptadiagonal_destroy(mg_a_poisson(l))
            call geometry_subdomain_destroy(mg_sdm(l))
            call geometry_subdomain_ddt_destroy(mg_sdm(l))
        enddo

        deallocate(mg_sdm)

    end subroutine multigrid_destroy

    !>
    !> @brief       Restriction subroutine for CGS, CGA, and CGPSA
    !> @detail      This subroutine is commonly used by using stride_f and offest_f variables
    !> @param       val_c       Variable in coarser grid
    !> @param       val_f       Variable in finer grid
    !> @param       dm_c        MG subdomain in coarser level
    !> @param       dm_f        MG subdomain in finer level
    !> @param       level       MG level where restriction is conducted
    !>
    subroutine multigrid_restriction(val_c, val_f, dm_c, dm_f, level)

        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z, myrank
        implicit none

        real(kind=8), intent(inout)     :: val_c(0:,0:,0:)
        real(kind=8), intent(in)        :: val_f(0:,0:,0:)
        type(subdomain), intent(in)     :: dm_c, dm_f
        integer(kind=4), intent(in)     :: level

        integer(kind=4)                 :: i, j, k
        integer(kind=4)                 :: i_c, j_c, k_c
        integer(kind=4)                 :: iz_f, jz_f, kz_f
        integer(kind=4)                 :: ip_f, jp_f, kp_f
        integer(kind=4)                 :: i_stride_f, i_offset_f
        integer(kind=4)                 :: j_stride_f, j_offset_f
        integer(kind=4)                 :: k_stride_f, k_offset_f
        integer(kind=4)                 :: nx_c, ny_c, nz_c
        real(kind=8)                    :: vol_f(2,2,2), vol_c
        real(kind=8), allocatable       :: dxf(:), dyf(:), dzf(:)

        val_c = 0.0d0

        allocate(dxf(0:dm_f%nx+1))
        allocate(dyf(0:dm_f%ny+1))
        allocate(dzf(0:dm_f%nz+1))

        dxf = 0.0d0
        dyf = 0.0d0
        dzf = 0.0d0

        select case(aggregation_type)
            case(0)
                nx_c = dm_c%nx
                ny_c = dm_c%ny
                nz_c = dm_c%nz
            case(1)
                if(level.eq.lv_aggregation-1) then
                    ! The actual number of grid for restriction at aggregation level.
                    nx_c = dm_c%nx / comm_1d_x%nprocs
                    ny_c = dm_c%ny / comm_1d_y%nprocs
                    nz_c = dm_c%nz / comm_1d_z%nprocs
                else
                    nx_c = dm_c%nx
                    ny_c = dm_c%ny
                    nz_c = dm_c%nz
                endif
            case(2)
                if(level.eq.lv_aggregation_x-1) then
                    ! The actual number of grid for restriction at aggregation level.
                    nx_c = dm_c%nx / comm_1d_x%nprocs
                else
                    nx_c = dm_c%nx
                endif
                if(level.eq.lv_aggregation_y-1) then
                    ! The actual number of grid for restriction at aggregation level.
                    ny_c = dm_c%ny / comm_1d_y%nprocs
                else
                    ny_c = dm_c%ny
                endif
                if(level.eq.lv_aggregation_z-1) then
                    ! The actual number of grid for restriction at aggregation level.
                    nz_c = dm_c%nz / comm_1d_z%nprocs
                else
                    nz_c = dm_c%nz
                endif
            case default
                if(myrank.eq.0) print '(a,i2)', '[Error] Aggregation method should be 0, 1, or 2.', aggregation_type
        end select

        if(level.ge.lv_gdm_coarsest_x) then
            ! Restriction beyond coarsest level, index = i * 1 (i_stride_f) - 0 (i_offset_f)
            dxf(:) = dm_f%dxm(:)
            i_stride_f = 1
            i_offset_f = 0
        else
            ! Normal restriction until coarsest level, index = i * 2 (i_stride_f) - 1 (i_offset_f)
            do i = 1, dm_f%nx
                i_c = int((i-1)/2)+1
                dxf(i) = dm_c%dxm(i_c)-dm_f%dxm(i)
            enddo
            i_stride_f = 2
            i_offset_f = 1
        endif

        if(level.ge.lv_gdm_coarsest_y) then
            ! Restriction beyond coarsest level, index = j * 1 (j_stride_f) - 0 (j_offset_f)
            dyf(:) = dm_f%dym(:)
            j_stride_f = 1
            j_offset_f = 0
        else
            ! Normal restriction until coarsest level, index = j * 2 (j_stride_f) - 1 (j_offset_f)
            do j = 1, dm_f%ny
                j_c = int((j-1)/2)+1
                dyf(j) = dm_c%dym(j_c)-dm_f%dym(j)
            enddo
            j_stride_f = 2
            j_offset_f = 1
        endif

        if(level.ge.lv_gdm_coarsest_z) then
            ! Restriction beyond coarsest level, index = k * 1 (k_stride_f) - 0 (k_offset_f)
            dzf(:) = dm_f%dzm(:)
            k_stride_f = 1
            k_offset_f = 0
        else
            ! Normal restriction until coarsest level, index = k * 2 (k_stride_f) - 1 (k_offset_f)
            do k = 1, dm_f%nz
                k_c = int((k-1)/2)+1
                dzf(k) = dm_c%dzm(k_c)-dm_f%dzm(k)
            enddo
            k_stride_f = 2
            k_offset_f = 1
        endif

        ! Restriction base on volume. Beyond coarsest level, val_c = val_f (vol_c = 8 * vol_f and vol_f are same)
!$omp parallel do private(kp_f, kz_f, jp_f, jz_f, ip_f, iz_f, vol_f, vol_c) &
!$omp firstprivate(k_stride_f,k_offset_f,j_stride_f,j_offset_f,i_stride_f,i_offset_f) &
!$omp default(shared)
        do k = 1, nz_c
            kp_f = k * k_stride_f
            kz_f = kp_f - k_offset_f
            do j = 1, ny_c
                jp_f = j * j_stride_f
                jz_f = jp_f - j_offset_f
                do i = 1, nx_c
                    ip_f = i * i_stride_f
                    iz_f = ip_f - i_offset_f

                    vol_f(1,1,1) = dxf(iz_f) * dyf(jz_f) * dzf(kz_f)
                    vol_f(2,1,1) = dxf(ip_f) * dyf(jz_f) * dzf(kz_f)
                    vol_f(1,2,1) = dxf(iz_f) * dyf(jp_f) * dzf(kz_f)
                    vol_f(2,2,1) = dxf(ip_f) * dyf(jp_f) * dzf(kz_f)
                    vol_f(1,1,2) = dxf(iz_f) * dyf(jz_f) * dzf(kp_f)
                    vol_f(2,1,2) = dxf(ip_f) * dyf(jz_f) * dzf(kp_f)
                    vol_f(1,2,2) = dxf(iz_f) * dyf(jp_f) * dzf(kp_f)
                    vol_f(2,2,2) = dxf(ip_f) * dyf(jp_f) * dzf(kp_f)
                    vol_c        = sum(vol_f)
                        
                    val_c(i, j, k) = (vol_f(1,1,1) * val_f(ip_f, jp_f, kp_f ) &
                                    + vol_f(2,1,1) * val_f(iz_f, jp_f, kp_f ) &
                                    + vol_f(1,2,1) * val_f(ip_f, jz_f, kp_f ) &
                                    + vol_f(2,2,1) * val_f(iz_f, jz_f, kp_f ) &
                                    + vol_f(1,1,2) * val_f(ip_f, jp_f, kz_f ) &
                                    + vol_f(2,1,2) * val_f(iz_f, jp_f, kz_f ) &
                                    + vol_f(1,2,2) * val_f(ip_f, jz_f, kz_f ) &
                                    + vol_f(2,2,2) * val_f(iz_f, jz_f, kz_f ) ) / vol_c
                enddo
            enddo
        enddo
        deallocate(dxf , dyf, dzf)

        if(myrank.eq.0) print '(a,i2,a,i2)', '[MG] Restricton from level ',level,' to level ',level+1

#ifdef DEBUG_RESTRICTION
        call multigrid_debug_print_restriction_info(val_c, val_f, vol_c, vol_f, level, nx_c, ny_c, nz_c, &
                                                    i_stride_f, j_stride_f, k_stride_f, i_offset_f, j_offset_f, k_offset_f)
#endif
    
    end subroutine multigrid_restriction

    !>
    !> @brief       Prolongation subroutine for CGS, CGA, and CGPSA
    !> @detail      This subroutine is commonly used by using stride_f and offest_f variables
    !>              Constant profile in assumed in finer level
    !> @param       val_f       Variable in finer grid
    !> @param       val_c       Variable in coarser grid
    !> @param       dm_f        MG subdomain in finer level
    !> @param       dm_c        MG subdomain in coarser level
    !> @param       level       MG level where restriction is conducted
    !>
    subroutine multigrid_prolongation_linear_on_nonuniform_grid(val_f, val_c, dm_f, dm_c, level)

        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z, myrank

        implicit none

        real(kind=8), intent(out)       :: val_f(0:,0:,0:)
        real(kind=8), intent(in)        :: val_c(0:,0:,0:)
        type(subdomain), intent(in)     :: dm_c
        type(subdomain), intent(in)     :: dm_f
        integer(kind=4), intent(in)     :: level

        integer(kind=4)                 :: nx_c, ny_c, nz_c
        integer(kind=4)                 :: i, j, k
        integer(kind=4)                 :: iz_f, jz_f, kz_f
        integer(kind=4)                 :: ip_f, jp_f, kp_f
        integer(kind=4)                 :: i_stride_f, i_offset_f
        integer(kind=4)                 :: j_stride_f, j_offset_f
        integer(kind=4)                 :: k_stride_f, k_offset_f
        integer(kind=4)                 :: im_c, jm_c, km_c
        integer(kind=4)                 :: ip_c, jp_c, kp_c
        real(kind=8), allocatable       :: dxp(:), dyp(:), dzp(:)
        real(kind=8), allocatable       :: dxn(:), dyn(:), dzn(:)

        allocate(dxp(0:dm_f%nx+1))
        allocate(dyp(0:dm_f%ny+1))
        allocate(dzp(0:dm_f%nz+1))

        allocate(dxn(0:dm_f%nx+1))
        allocate(dyn(0:dm_f%ny+1))
        allocate(dzn(0:dm_f%nz+1))

        dxp = 0.0d0
        dyp = 0.0d0
        dzp = 0.0d0

        dxn = 0.0d0
        dyn = 0.0d0
        dzn = 0.0d0

        select case(aggregation_type)
            case(0)
                nx_c = dm_c%nx
                ny_c = dm_c%ny
                nz_c = dm_c%nz
            case(1)
                if(level.eq.lv_aggregation-1) then
                    ! The actual number of grid for restriction at aggregation level.
                    nx_c = dm_c%nx / comm_1d_x%nprocs
                    ny_c = dm_c%ny / comm_1d_y%nprocs
                    nz_c = dm_c%nz / comm_1d_z%nprocs
                else
                    nx_c = dm_c%nx
                    ny_c = dm_c%ny
                    nz_c = dm_c%nz
                endif
            case(2)
                if(level.eq.lv_aggregation_x-1) then
                    ! The actual number of grid for restriction at aggregation level.
                    nx_c = dm_c%nx / comm_1d_x%nprocs
                else
                    nx_c = dm_c%nx
                endif
                if(level.eq.lv_aggregation_y-1) then
                    ny_c = dm_c%ny / comm_1d_y%nprocs
                else
                    ny_c = dm_c%ny
                endif
                if(level.eq.lv_aggregation_z-1) then
                    ! The actual number of grid for restriction at aggregation level.
                    nz_c = dm_c%nz / comm_1d_z%nprocs
                else
                    nz_c = dm_c%nz
                endif
            case default
                if(myrank.eq.0) print '(a,i2)', '[Error] Aggregation method should be 0, 1, or 2.', aggregation_type
        end select

        if(level.ge.lv_gdm_coarsest_x) then
            ! Prolongtaion beyond coarsest level, index = i * 1 (i_stride_f) - 0 (i_offset_f)
            i_stride_f = 1
            i_offset_f = 0
            dxp(:) = dm_c%dxm(:)
            dxn(:) = dm_c%dxm(:)
        else
            ! Normal prolongation until coarsest level, index = i * 2 (i_stride_f) - 1 (i_offset_f)
            i_stride_f = 2
            i_offset_f = 1
            do i = 1, nx_c
                ip_f = 2 * i
                dxp(ip_f-1)= 1.0d0*dm_f%dxm(ip_f-2) + 0.5d0*dm_f%dxm(ip_f-1)
                dxp(ip_f)  = 0.5d0*dm_f%dxm(ip_f-1)
                dxn(ip_f-1)= 0.5d0*dm_f%dxm(ip_f  )
                dxn(ip_f)  = 1.0d0*dm_f%dxm(ip_f+1) + 0.5d0*dm_f%dxm(ip_f  )
            enddo
        endif

        if(level.ge.lv_gdm_coarsest_y) then
            ! Prolongtaion beyond coarsest level, index = j * 1 (j_stride_f) - 0 (j_offset_f)
            j_stride_f = 1
            j_offset_f = 0
            dyp(:) = dm_c%dym(:)
            dyn(:) = dm_c%dym(:)
        else
            ! Normal prolongation until coarsest level, index = j * 2 (j_stride_f) - 1 (j_offset_f)
            j_stride_f = 2
            j_offset_f = 1
            do j = 1, ny_c
                jp_f = 2 * j
                dyp(jp_f-1) = 1.0d0*dm_f%dym(jp_f-2) + 0.5d0*dm_f%dym(jp_f-1)
                dyp(jp_f)   = 0.5d0*dm_f%dym(jp_f-1)
                dyn(jp_f-1) = 0.5d0*dm_f%dym(jp_f  )
                dyn(jp_f)   = 1.0d0*dm_f%dym(jp_f+1) + 0.5d0*dm_f%dym(jp_f  )
            enddo
        endif

        if(level.ge.lv_gdm_coarsest_z) then
            ! Prolongtaion beyond coarsest level, index = k * 1 (k_stride_f) - 0 (k_offset_f)
            k_stride_f = 1
            k_offset_f = 0
            dzp(:) = dm_c%dzm(:)
            dzn(:) = dm_c%dzm(:)
        else
            ! Normal prolongation until coarsest level, index = k * 2 (k_stride_f) - 1 (k_offset_f)
            k_stride_f = 2
            k_offset_f = 1
            do k = 1, nz_c
                kp_f = 2 * k
                dzp(kp_f-1) = 1.0d0*dm_f%dzm(kp_f-2) + 0.5d0*dm_f%dzm(kp_f-1)
                dzp(kp_f)   = 0.5d0*dm_f%dzm(kp_f-1)
                dzn(kp_f-1) = 0.5d0*dm_f%dzm(kp_f  )
                dzn(kp_f)   = 1.0d0*dm_f%dzm(kp_f+1) + 0.5d0*dm_f%dzm(kp_f  )
            enddo
        endif

        ! Prolongation based on volume. Beyond coarsest level, val_c = val_f (vol_c = 8 * vol_f and vol_f are same)
!$omp parallel do private(kp_f, kz_f, km_c, kp_c, jp_f, jz_f, jm_c, jp_c, ip_f, iz_f, im_c, ip_c) &
!$omp firstprivate(k_stride_f,k_offset_f,j_stride_f,j_offset_f,i_stride_f,i_offset_f) &
!$omp default(shared)
        do k = 1, nz_c

            kp_f = k * k_stride_f
            kz_f = kp_f - k_offset_f
            km_c = k - k_offset_f
            kp_c = k + k_offset_f

            do j = 1, ny_c
            
                jp_f = j * j_stride_f
                jz_f = jp_f - j_offset_f
                jm_c = j - j_offset_f
                jp_c = j + j_offset_f

                do i = 1, nx_c

                    ip_f = i * i_stride_f
                    iz_f = ip_f - i_offset_f
                    im_c = i - i_offset_f
                    ip_c = i + i_offset_f

                    val_f(iz_f, jz_f, kz_f) =   (dxp(iz_f)*dyp(jz_f)*dzp(kz_f)*val_c(i, j, k) &
                                                +dxp(ip_f)*dyp(jz_f)*dzp(kz_f)*val_c(im_c, j, k) &
                                                +dxp(iz_f)*dyp(jp_f)*dzp(kz_f)*val_c(i, jm_c, k) &
                                                +dxp(iz_f)*dyp(jz_f)*dzp(kp_f)*val_c(i, j, km_c) &
                                                +dxp(ip_f)*dyp(jp_f)*dzp(kz_f)*val_c(im_c, jm_c, k) &
                                                +dxp(ip_f)*dyp(jz_f)*dzp(kp_f)*val_c(im_c, j, km_c) &
                                                +dxp(iz_f)*dyp(jp_f)*dzp(kp_f)*val_c(i, jm_c, km_c) &
                                                +dxp(ip_f)*dyp(jp_f)*dzp(kp_f)*val_c(im_c, jm_c, km_c) ) &
                                            / ( (dxp(iz_f)+dxp(ip_f))*(dyp(jz_f)+dyp(jp_f))*(dzp(kz_f)+dzp(kp_f)) )

                    val_f(ip_f, jz_f, kz_f) =   (dxn(ip_f)*dyp(jz_f)*dzp(kz_f)*val_c(i, j, k) &
                                                +dxn(iz_f)*dyp(jz_f)*dzp(kz_f)*val_c(ip_c, j, k) &
                                                +dxn(ip_f)*dyp(jp_f)*dzp(kz_f)*val_c(i, jm_c, k) &
                                                +dxn(ip_f)*dyp(jz_f)*dzp(kp_f)*val_c(i, j, km_c) &
                                                +dxn(iz_f)*dyp(jp_f)*dzp(kz_f)*val_c(ip_c, jm_c, k) &
                                                +dxn(iz_f)*dyp(jz_f)*dzp(kp_f)*val_c(ip_c, j, km_c) &
                                                +dxn(ip_f)*dyp(jp_f)*dzp(kp_f)*val_c(i, jm_c, km_c) &
                                                +dxn(iz_f)*dyp(jp_f)*dzp(kp_f)*val_c(ip_c, jm_c, km_c) ) &
                                            / ( (dxn(iz_f)+dxn(ip_f))*(dyp(jz_f)+dyp(jp_f))*(dzp(kz_f)+dzp(kp_f)) )

                    val_f(iz_f, jp_f, kz_f)  =  (dxp(iz_f)*dyn(jp_f)*dzp(kz_f)*val_c(i, j, k) &
                                                +dxp(ip_f)*dyn(jp_f)*dzp(kz_f)*val_c(im_c, j, k) &
                                                +dxp(iz_f)*dyn(jz_f)*dzp(kz_f)*val_c(i, jp_c, k) &
                                                +dxp(iz_f)*dyn(jp_f)*dzp(kp_f)*val_c(i, j, km_c) &
                                                +dxp(ip_f)*dyn(jz_f)*dzp(kz_f)*val_c(im_c, jp_c, k) &
                                                +dxp(ip_f)*dyn(jp_f)*dzp(kp_f)*val_c(im_c, j, km_c) &
                                                +dxp(iz_f)*dyn(jz_f)*dzp(kp_f)*val_c(i, jp_c, km_c) &
                                                +dxp(ip_f)*dyn(jz_f)*dzp(kp_f)*val_c(im_c, jp_c, km_c) ) &
                                            / ( (dxp(iz_f)+dxp(ip_f))*(dyn(jz_f)+dyn(jp_f))*(dzp(kz_f)+dzp(kp_f)) )

                    val_f(ip_f, jp_f, kz_f)  =  (dxn(ip_f)*dyn(jp_f)*dzp(kz_f)*val_c(i, j, k) &
                                                +dxn(iz_f)*dyn(jp_f)*dzp(kz_f)*val_c(ip_c, j, k) &
                                                +dxn(ip_f)*dyn(jz_f)*dzp(kz_f)*val_c(i, jp_c, k) &
                                                +dxn(ip_f)*dyn(jp_f)*dzp(kp_f)*val_c(i, j, km_c) &
                                                +dxn(iz_f)*dyn(jz_f)*dzp(kz_f)*val_c(ip_c, jp_c, k) &
                                                +dxn(iz_f)*dyn(jp_f)*dzp(kp_f)*val_c(ip_c, j, km_c) &
                                                +dxn(ip_f)*dyn(jz_f)*dzp(kp_f)*val_c(i, jp_c, km_c) &
                                                +dxn(iz_f)*dyn(jz_f)*dzp(kp_f)*val_c(ip_c, jp_c, km_c) ) &
                                            / ( (dxn(iz_f)+dxn(ip_f))*(dyn(jz_f)+dyn(jp_f))*(dzp(kz_f)+dzp(kp_f)) )

                    val_f(iz_f, jz_f, kp_f)  =  (dxp(iz_f)*dyp(jz_f)*dzn(kp_f)*val_c(i, j, k) &
                                                +dxp(ip_f)*dyp(jz_f)*dzn(kp_f)*val_c(im_c, j, k) &
                                                +dxp(iz_f)*dyp(jp_f)*dzn(kp_f)*val_c(i, jm_c, k) &
                                                +dxp(iz_f)*dyp(jz_f)*dzn(kz_f)*val_c(i, j, kp_c) &
                                                +dxp(ip_f)*dyp(jp_f)*dzn(kp_f)*val_c(im_c, jm_c, k) &
                                                +dxp(ip_f)*dyp(jz_f)*dzn(kz_f)*val_c(im_c, j, kp_c) &
                                                +dxp(iz_f)*dyp(jp_f)*dzn(kz_f)*val_c(i, jm_c, kp_c) &
                                                +dxp(ip_f)*dyp(jp_f)*dzn(kz_f)*val_c(im_c, jm_c, kp_c) ) &
                                            / ( (dxp(iz_f)+dxp(ip_f))*(dyp(jz_f)+dyp(jp_f))*(dzn(kz_f)+dzn(kp_f)) )

                    val_f(ip_f, jz_f, kp_f)  =  (dxn(ip_f)*dyp(jz_f)*dzn(kp_f)*val_c(i, j, k) &
                                                +dxn(iz_f)*dyp(jz_f)*dzn(kp_f)*val_c(ip_c, j, k) &
                                                +dxn(ip_f)*dyp(jp_f)*dzn(kp_f)*val_c(i, jm_c, k) &
                                                +dxn(ip_f)*dyp(jz_f)*dzn(kz_f)*val_c(i, j, kp_c) &
                                                +dxn(iz_f)*dyp(jp_f)*dzn(kp_f)*val_c(ip_c, jm_c, k) &
                                                +dxn(iz_f)*dyp(jz_f)*dzn(kz_f)*val_c(ip_c, j, kp_c) &
                                                +dxn(ip_f)*dyp(jp_f)*dzn(kz_f)*val_c(i, jm_c, kp_c) &
                                                +dxn(iz_f)*dyp(jp_f)*dzn(kz_f)*val_c(ip_c, jm_c, kp_c) ) &
                                            / ( (dxn(iz_f)+dxn(ip_f))*(dyp(jz_f)+dyp(jp_f))*(dzn(kz_f)+dzn(kp_f)) )

                    val_f(iz_f, jp_f, kp_f)  =  (dxp(iz_f)*dyn(jp_f)*dzn(kp_f)*val_c(i, j, k) &
                                                +dxp(ip_f)*dyn(jp_f)*dzn(kp_f)*val_c(im_c, j, k) &
                                                +dxp(iz_f)*dyn(jz_f)*dzn(kp_f)*val_c(i, jp_c, k) &
                                                +dxp(iz_f)*dyn(jp_f)*dzn(kz_f)*val_c(i, j, kp_c) &
                                                +dxp(ip_f)*dyn(jz_f)*dzn(kp_f)*val_c(im_c, jp_c, k) &
                                                +dxp(ip_f)*dyn(jp_f)*dzn(kz_f)*val_c(im_c, j, kp_c) &
                                                +dxp(iz_f)*dyn(jz_f)*dzn(kz_f)*val_c(i, jp_c, kp_c) &
                                                +dxp(ip_f)*dyn(jz_f)*dzn(kz_f)*val_c(im_c, jp_c, kp_c) ) &
                                            / ( (dxp(iz_f)+dxp(ip_f))*(dyn(jz_f)+dyn(jp_f))*(dzn(kz_f)+dzn(kp_f)) )

                    val_f(ip_f, jp_f, kp_f) =   (dxn(ip_f)*dyn(jp_f)*dzn(kp_f)*val_c(i, j, k) &
                                                +dxn(iz_f)*dyn(jp_f)*dzn(kp_f)*val_c(ip_c, j, k) &
                                                +dxn(ip_f)*dyn(jz_f)*dzn(kp_f)*val_c(i, jp_c, k) &
                                                +dxn(ip_f)*dyn(jp_f)*dzn(kz_f)*val_c(i, j, kp_c) &
                                                +dxn(iz_f)*dyn(jz_f)*dzn(kp_f)*val_c(ip_c, jp_c, k) &
                                                +dxn(iz_f)*dyn(jp_f)*dzn(kz_f)*val_c(ip_c, j, kp_c) &
                                                +dxn(ip_f)*dyn(jz_f)*dzn(kz_f)*val_c(i, jp_c, kp_c) &
                                                +dxn(iz_f)*dyn(jz_f)*dzn(kz_f)*val_c(ip_c, jp_c, kp_c) ) &
                                            / ( (dxn(iz_f)+dxn(ip_f))*(dyn(jz_f)+dyn(jp_f))*(dzn(kz_f)+dzn(kp_f)) )

                enddo
            enddo
        enddo

        deallocate(dxn, dxp)
        deallocate(dyn, dyp)
        deallocate(dzn, dzp)

        if(myrank.eq.0) print '(a,i2,a,i2)', '[MG] Prolongation from level ',level+1,' to level ',level

#ifdef DEBUG_PROLONGATION
        call multigrid_debug_print_prolongation_info(val_c, val_f, level, dm_f%nx, dm_f%ny, dm_f%nz, i_gl_c, j_gl_c, k_gl_c, &
                                                    i_stride_f, j_stride_f, k_stride_f, i_offset_f, j_offset_f, k_offset_f)
#endif

    end subroutine multigrid_prolongation_linear_on_nonuniform_grid

    !>
    !> @brief       Residual calculation. rsd = rhs - a_poisson * x
    !> @param       rsd             Residual
    !> @param       a_poisson       A matrix
    !> @param       x               x vector
    !> @param       rhs             RHS
    !> @param       dm              MG subdomain in target level
    !> @param       is_aggregated   Aggregation information on whether domain is aggregated (.true.) or not (.false.)
    !>
    subroutine multigrid_residual(rsd, a_poisson, x, rhs, dm, is_aggregated)

        implicit none

        real(kind=8), intent(out)               :: rsd(0:,0:,0:)
        type(matrix_heptadiagonal), intent(in)  :: a_poisson
        real(kind=8), intent(inout)             :: x(0:,0:,0:)
        real(kind=8), intent(in)                :: rhs(0:,0:,0:)
        type(subdomain), intent(in)             :: dm
        logical, intent(in)                     :: is_aggregated(0:2)

        integer(kind=4)                 :: i, j, k

        rsd(:,:,:) = 0.0d0
        call mv_mul_poisson_matrix(rsd, a_poisson, x, dm, is_aggregated)

!$omp parallel do shared(rsd,rhs)
        do k = 1, dm%nz
            do j = 1, dm%ny
                do i = 1, dm%nx
                    rsd(i,j,k) = rhs(i,j,k) - rsd(i,j,k)
                enddo
            enddo
        enddo

    end subroutine multigrid_residual

    !>
    !> @brief       Coarsest level solver. Sover a_poisson * x = rhs
    !> @param       x               Solution vector
    !> @param       a_poisson       A matrix
    !> @param       rhs             RHS
    !> @param       maxiteration    Maximum number of iterations
    !> @param       tolerance       Convergence criteria
    !> @param       omega           Relexation factor
    !> @param       is_aggregated   Boolean whether domain is aggregated (.true.) or not (.false.)
    !>
    subroutine multigrid_solve_coarset_level(x, a_poisson, rhs, dm, maxiteration, tolerance, omega, is_aggregated)

        use mpi_topology, only  : myrank

        implicit none

        real(kind=8), intent(inout)             :: x(0:,0:,0:)
        type(matrix_heptadiagonal), intent(in)  :: a_poisson
        real(kind=8), intent(in)                :: rhs(0:,0:,0:)
        type(subdomain), intent(in)             :: dm
        integer(kind=4), intent(in)             :: maxiteration
        real(kind=8), intent(in)                :: tolerance
        real(kind=8), intent(in)                :: omega
        logical, intent(in)                     :: is_aggregated(0:2)
        integer(kind=4)                         :: i, j, k

        x = 0.0d0

        if((dm%nx*dm%ny*dm%nz.eq.1).and.(all(dm%is_aggregated).eq..true.)) then
            ! When the number of grid is equal to 1.
            if(myrank.eq.0) print '(a)', '[MG] Obtain the solution on a single grid in the coarsest level.'
            x(1,1,1) = rhs(1,1,1) / a_poisson%coeff(0,1,1,1)
        else
            ! When the number of grid is not equal to 1, use solver.
            if(myrank.eq.0) print '(a)', '[MG] Obtain the solution on a grid in the coarsest level with RGBS solver.'
            call rbgs_solver_poisson_matrix(x, &
                                            a_poisson, &
                                            rhs, &
                                            dm, maxiteration, tolerance, omega, &
                                            is_aggregated)
        endif

    end subroutine multigrid_solve_coarset_level

    !>
    !> @brief       V-Cycle MG solver: call subroutines according to aggregation type
    !> @param       sol             Result solution having a shape of 3D matrix
    !> @param       rsd             Residual
    !> @param       a_poisson       Heptadiagonal poisson matrix
    !> @param       rhs             RHS vector having a shape of 3D matrix
    !> @param       sdm             Subdomain
    !> @param       maxiteration    Maximum number of iterations
    !> @param       tolerance       Convergence criteria
    !> @param       omega_sor       Relexation factor
    !>
    subroutine multigrid_solve_vcycle(sol, rsd, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor)

        use mpi
        use mpi_topology, only : myrank

        implicit none

        real(kind=8), intent(inout)             :: sol(0:,0:,0:), rsd(0:,0:,0:)
        type(matrix_heptadiagonal), intent(in)  :: a_poisson
        real(kind=8), intent(in)                :: rhs(0:,0:,0:)
        type(subdomain), intent(in)             :: sdm
        integer(kind=4), intent(in)             :: maxiteration
        real(kind=8), intent(in)                :: tolerance
        real(kind=8), intent(in)                :: omega_sor

        integer(kind=4)                         :: ierr

        select case(aggregation_type)
        case(0)
            call multigrid_CGS_vcycle_solver(sol, rsd, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor)
        case(1)
            call multigrid_CGA_vcycle_solver(sol, rsd, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor)
        case(2)
            call multigrid_CGPSA_vcycle_solver(sol, rsd, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor)
        case default
            if(myrank.eq.0) print '(a,i2)', '[Error] Aggregation method should be 0, 1, or 2.', aggregation_type
            call MPI_Finalize(ierr)
            stop
        end select

    end subroutine multigrid_solve_vcycle

    !>
    !> @brief       V-Cycle MG solver with CGS
    !> @param       sol             Result solution having a shape of 3D matrix
    !> @param       rsd             Residual
    !> @param       a_poisson       Heptadiagonal poisson matrix
    !> @param       rhs             RHS vector having a shape of 3D matrix
    !> @param       sdm             Subdomain
    !> @param       maxiteration    Maximum number of iterations
    !> @param       tolerance       Convergence criteria
    !> @param       omega_sor       Relexation factor
    !>
    subroutine multigrid_CGS_vcycle_solver(sol, rsd, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor)

        use mpi_topology, only  : myrank

        implicit none

        real(kind=8), intent(inout)                 :: sol(0:,0:,0:), rsd(0:,0:,0:)
        type(matrix_heptadiagonal), intent(in)      :: a_poisson
        real(kind=8), intent(in)                    :: rhs(0:,0:,0:)
        type(subdomain), intent(in)                 :: sdm
        integer(kind=4), intent(in)                 :: maxiteration
        real(kind=8), intent(in)                    :: tolerance
        real(kind=8), intent(in)                    :: omega_sor

        real(kind=8)                                :: rsd_val, res0tol
        integer(kind=4)     :: l, cyc

        call multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm%is_aggregated)
        call vv_dot_3d_matrix(res0tol, rsd, rsd, sdm%nx, sdm%ny, sdm%nz, sdm%is_aggregated)

        ! V-cycle Iterations
        do cyc = 1, n_vcycles

            ! Smooting and restriction
            call rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm%is_aggregated)
            call multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm%is_aggregated)
            call multigrid_restriction(mg_sdm(1)%b, rsd, mg_sdm(1), sdm, 0)

            ! Smooting and restriction
            do l = 1, n_levels-1
                mg_sdm(l)%x = 0.0d0
                call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                call multigrid_residual(mg_sdm(l)%r, mg_a_poisson(l), mg_sdm(l)%x, mg_sdm(l)%b, mg_sdm(l), mg_sdm(l)%is_aggregated)
                call multigrid_restriction(mg_sdm(l+1)%b, mg_sdm(l)%r, mg_sdm(l+1), mg_sdm(l),l)
            enddo

            ! Solve at the coarset level
            call multigrid_solve_coarset_level( mg_sdm(n_levels)%x, &
                                                    mg_a_poisson(n_levels), &
                                                    mg_sdm(n_levels)%b, &
                                                    mg_sdm(n_levels), &
                                                    1000, tolerance, omega_sor, &
                                                    mg_sdm(n_levels)%is_aggregated)

            call multigrid_residual(mg_sdm(n_levels)%r, mg_a_poisson(n_levels), mg_sdm(n_levels)%x, mg_sdm(n_levels)%b, mg_sdm(n_levels), mg_sdm(n_levels)%is_aggregated)
            if(myrank.eq.0) print '(a,e18.10,a,e18.10)','[MG] Solution in the coarset level : x(1,1,1) = ',mg_sdm(n_levels)%x(1,1,1),', residue = ',mg_sdm(n_levels)%r(1,1,1)
#ifdef DEBUG_COARSEST
            call multigrid_debug_print_coarsest_level_solution(cyc, mg_sdm(n_levels))
#endif
            ! Prolongation and smoothing: Ghostcell update is required for prolongation
            do l = n_levels-1, 1, -1
                call multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm(l)%r, mg_sdm(l+1)%x, mg_sdm(l), mg_sdm(l+1),l)
                mg_sdm(l)%x = mg_sdm(l)%x + mg_sdm(l)%r
                call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                call geometry_halocell_update_selectively(mg_sdm(l)%x, mg_sdm(l), mg_sdm(1)%is_aggregated)
            enddo

            ! Prolongation and smoothing
            call multigrid_prolongation_linear_on_nonuniform_grid(rsd, mg_sdm(1)%x, sdm, mg_sdm(1), 0)
            sol = sol + rsd
            call rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, mg_sdm(l)%is_aggregated)

            ! Convergence check
            call multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm%is_aggregated)
            call vv_dot_3d_matrix(rsd_val, rsd, rsd, sdm%nx, sdm%ny, sdm%nz, sdm%is_aggregated)
            if(myrank.eq.0) print '(a,i4,a,3(e15.7,x))', '[MG] cycle = ',cyc,', Error:',sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol)
            if (sqrt(rsd_val/res0tol) .LT. tolerance) exit 
        enddo

    end subroutine multigrid_CGS_vcycle_solver

    !>
    !> @brief       V-Cycle MG solver with CGA
    !> @detail      MPI_Gather/MPI_Scatterv with derived datatypes are used for aggregation.
    !> @param       sol             Result solution having a shape of 3D matrix
    !> @param       rsd             Residual
    !> @param       a_poisson       Heptadiagonal poisson matrix
    !> @param       rhs             RHS vector having a shape of 3D matrix
    !> @param       sdm             Subdomain
    !> @param       maxiteration    Maximum number of iterations
    !> @param       tolerance       Convergence criteria
    !> @param       omega_sor       Relexation factor
    !>
    subroutine multigrid_CGA_vcycle_solver(sol, rsd, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor)

        use mpi
        use mpi_topology, only  : myrank, nprocs, mpi_world_cart, comm_1d_x, comm_1d_y, comm_1d_z

        implicit none

        real(kind=8), intent(inout)                 :: sol(0:,0:,0:), rsd(0:,0:,0:)
        type(matrix_heptadiagonal), intent(in)      :: a_poisson
        real(kind=8), intent(in)                    :: rhs(0:,0:,0:)
        type(subdomain), intent(in)                 :: sdm
        integer(kind=4), intent(in)                 :: maxiteration
        real(kind=8), intent(in)                    :: tolerance
        real(kind=8), intent(in)                    :: omega_sor

        real(kind=8)                                :: rsd_val, res0tol
        integer(kind=4)                             :: l, i, cyc
        integer(kind=4)                             :: ierr

        ! MPI_Gatherv/MPI_Scatterv used for aggregation
        integer(kind=4)                             :: ddtype_temp1
        integer(kind=4)                             :: ddtype_scatterv  !> Derived datatype for grid disaggregation
        integer(kind=4)                             :: ddtype_gatherv   !> Derived datatype for grid aggregation
        integer(kind=4)                             :: sizes(0:2), subsizes(0:2), starts(0:2)
        integer(kind=4)                             :: r8size
        integer(kind=4)                             :: nx_aggr, ny_aggr, nz_aggr
        integer(kind=4)                             :: nx, ny, nz
        integer(kind=MPI_ADDRESS_KIND)              :: extent, lb
        integer(kind=4), allocatable, dimension(:)  :: cnts
        integer(kind=4), allocatable, dimension(:)  :: disps
        integer(kind=4), allocatable, dimension(:,:):: cart_coord
#ifndef MPI_INPLACE
        real(kind=8), allocatable               :: tmp(:,:,:)
#endif

        ! Start of MPI communication setup
        allocate( cart_coord(0:2,0:nprocs-1) )
        do i = 0, nprocs-1
            call MPI_Cart_coords(mpi_world_cart, i, 3, cart_coord(:,i), ierr )
        enddo

        nx_aggr = mg_sdm(lv_aggregation)%nx
        ny_aggr = mg_sdm(lv_aggregation)%ny
        nz_aggr = mg_sdm(lv_aggregation)%nz

        nx = nx_aggr / comm_1d_x%nprocs
        ny = ny_aggr / comm_1d_y%nprocs
        nz = nz_aggr / comm_1d_z%nprocs

        ! Derivide datatype for grid aggregation using MPI_Gatherv
        sizes    = (/nx_aggr+2, ny_aggr+2, nz_aggr+2/)
        subsizes = (/nx, ny, nz/)
        starts   = (/1, 1, 1/)

        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, ddtype_temp1, ierr)
        CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
        lb = 0
        extent = r8size
        call MPI_Type_create_resized(ddtype_temp1, lb, extent, ddtype_gatherv, ierr)
        call MPI_Type_commit(ddtype_gatherv, ierr)

        ! Derivide datatype for grid disaggregation using MPI_Scatterv
        sizes    = (/nx_aggr+2, ny_aggr+2, nz_aggr+2/)
        subsizes = (/nx+2, ny+2, nz+2/)
        starts   = (/0, 0, 0/)

        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, ddtype_temp1, ierr)
        CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
        lb = 0
        extent = r8size
        call MPI_Type_create_resized(ddtype_temp1, lb, extent, ddtype_scatterv, ierr)
        call MPI_Type_commit(ddtype_scatterv, ierr)

        ! DDT counts and displacements
        allocate(cnts(0:nprocs-1))
        allocate(disps(0:nprocs-1))
        do i = 0, nprocs-1
            cnts(i) = 1
            disps(i) = nx * cart_coord(0,i)  &
                    + ny * cart_coord(1,i) * (nx_aggr + 2)  &
                    + nz * cart_coord(2,i) * (nx_aggr + 2) * (ny_aggr + 2)
        enddo
        ! End of MPI communication setup

        call multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm%is_aggregated)
        call vv_dot_3d_matrix(res0tol, rsd, rsd, sdm%nx, sdm%ny, sdm%nz, sdm%is_aggregated)

        ! V-cycle Iterations
        do cyc = 1, n_vcycles

            ! Smooting and restriction until aggregation level
            call rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm%is_aggregated)
            call multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm%is_aggregated)
            call multigrid_restriction(mg_sdm(1)%b, rsd, mg_sdm(1), sdm, 0)

            do l = 1, lv_aggregation-1
                mg_sdm(l)%x = 0.0d0
                call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                call multigrid_residual(mg_sdm(l)%r, mg_a_poisson(l), mg_sdm(l)%x, mg_sdm(l)%b, mg_sdm(l), mg_sdm(l)%is_aggregated)
                call multigrid_restriction(mg_sdm(l+1)%b, mg_sdm(l)%r, mg_sdm(l+1), mg_sdm(l),l)
            enddo

            ! Grid aggregation 
#ifdef MPI_INPLACE
            if(myrank.eq.0) then
                call MPI_Gatherv(MPI_IN_PLACE, 0, ddtype_gatherv, mg_sdm(lv_aggregation)%b, cnts, disps, ddtype_gatherv, 0, MPI_COMM_WORLD, ierr)
            else
                call MPI_Gatherv(mg_sdm(lv_aggregation)%b, 1, ddtype_gatherv, mg_sdm(lv_aggregation)%b, cnts, disps, ddtype_gatherv, 0, MPI_COMM_WORLD, ierr)
            endif
#else
            allocate(tmp(0:size(mg_sdm(lv_aggregation)%b, 1)-1,0:size(mg_sdm(lv_aggregation)%b,2)-1,0:size(mg_sdm(lv_aggregation)%b,3)-1))
            call MPI_Gatherv(mg_sdm(lv_aggregation)%b, 1, ddtype_gatherv, tmp, cnts, disps, ddtype_gatherv, 0, MPI_COMM_WORLD, ierr)
            mg_sdm(lv_aggregation)%b = tmp
#endif

            ! Smooting and restriction from aggregation level: Only myrank0 conducts MG steps by grid disaggregation
            if(myrank.eq.0) then
                do l = lv_aggregation, n_levels-1
                    mg_sdm(l)%x = 0.0d0
                    call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                    call multigrid_residual(mg_sdm(l)%r, mg_a_poisson(l), mg_sdm(l)%x, mg_sdm(l)%b, mg_sdm(l), mg_sdm(l)%is_aggregated)
                    call multigrid_restriction(mg_sdm(l+1)%b, mg_sdm(l)%r, mg_sdm(l+1), mg_sdm(l),l)
                enddo

                ! Solve at the coarset level
                call multigrid_solve_coarset_level( mg_sdm(n_levels)%x, &
                                                        mg_a_poisson(n_levels), &
                                                        mg_sdm(n_levels)%b, &
                                                        mg_sdm(n_levels), &
                                                        1000, tolerance, omega_sor, &
                                                        mg_sdm(n_levels)%is_aggregated)

                call multigrid_residual(mg_sdm(n_levels)%r, mg_a_poisson(n_levels), mg_sdm(n_levels)%x, mg_sdm(n_levels)%b, mg_sdm(n_levels), mg_sdm(n_levels)%is_aggregated)
                print '(a,e18.10,a,e18.10)','[MG] Solution in the coarset level : x(1,1,1) = ',mg_sdm(n_levels)%x(1,1,1),', residue = ',mg_sdm(n_levels)%r(1,1,1)
#ifdef DEBUG_COARSEST
                call multigrid_debug_print_coarsest_level_solution(cyc, mg_sdm(n_levels))
#endif
                ! Prolongation and smoothing : Ghostcell update is required for prolongation
                do l = n_levels-1, lv_aggregation, -1
                    call multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm(l)%r, mg_sdm(l+1)%x, mg_sdm(l), mg_sdm(l+1),l)
                    mg_sdm(l)%x = mg_sdm(l)%x + mg_sdm(l)%r
                    call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                enddo
            endif

            ! Grid disaggregation: From here, all ranks conducts MG steps
#ifdef MPI_INPLACE
            if(myrank.eq.0) then
                call MPI_Scatterv(mg_sdm(lv_aggregation)%x, cnts, disps, ddtype_scatterv, MPI_IN_PLACE, 0, ddtype_scatterv, 0, MPI_COMM_WORLD, ierr)
            else
                call MPI_Scatterv(mg_sdm(lv_aggregation)%x, cnts, disps, ddtype_scatterv, mg_sdm(lv_aggregation)%x, 1, ddtype_scatterv, 0, MPI_COMM_WORLD, ierr)
            endif
#else
            call MPI_Scatterv(mg_sdm(lv_aggregation)%x, cnts, disps, ddtype_scatterv, tmp, 1, ddtype_scatterv, 0, MPI_COMM_WORLD, ierr)
            mg_sdm(lv_aggregation)%x = tmp
            deallocate(tmp)
#endif

            ! Prolongation and smoothing: Ghostcell update is required for prolongation
            do l = lv_aggregation-1, 1, -1
                call multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm(l)%r, mg_sdm(l+1)%x, mg_sdm(l), mg_sdm(l+1),l)
                mg_sdm(l)%x = mg_sdm(l)%x + mg_sdm(l)%r
                call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                call geometry_halocell_update_selectively(mg_sdm(l)%x, mg_sdm(l), mg_sdm(l)%is_aggregated)
            enddo

            call multigrid_prolongation_linear_on_nonuniform_grid(rsd, mg_sdm(1)%x, sdm, mg_sdm(1), 0)
            sol = sol + rsd
            call rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, mg_sdm(l)%is_aggregated)

            ! Convergence check
            call multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm%is_aggregated)
            call vv_dot_3d_matrix(rsd_val, rsd, rsd, sdm%nx, sdm%ny, sdm%nz, sdm%is_aggregated)
            if(myrank.eq.0) print '(a,i4,a,3(e15.7,x))', '[MG] cycle = ',cyc,', Error:',sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol)
            if (sqrt(rsd_val/res0tol) .LT. tolerance) exit 
        enddo

        deallocate( cnts )
        deallocate( disps )
        deallocate( cart_coord )

    end subroutine multigrid_CGA_vcycle_solver

    !>
    !> @brief       V-Cycle MG solver with CGPSA
    !> @detail      MPI_Gather/MPI_Scatterv with derived datatypes are used for aggregation.
    !> @param       sol             Result solution having a shape of 3D matrix
    !> @param       rsd             Residual
    !> @param       a_poisson       Heptadiagonal poisson matrix
    !> @param       rhs             RHS vector having a shape of 3D matrix
    !> @param       sdm             Subdomain
    !> @param       maxiteration    Maximum number of iterations
    !> @param       tolerance       Convergence criteria
    !> @param       omega_sor       Relexation factor
    !>
    subroutine multigrid_CGPSA_vcycle_solver(sol, rsd, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor)

        use mpi
        use mpi_topology, only  : myrank, cart_comm_1d, comm_1d_x, comm_1d_y, comm_1d_z

        implicit none

        real(kind=8), intent(inout)                 :: sol(0:,0:,0:), rsd(0:,0:,0:)
        type(matrix_heptadiagonal), intent(in)      :: a_poisson
        real(kind=8), intent(in)                    :: rhs(0:,0:,0:)
        type(subdomain), intent(in)                 :: sdm
        integer(kind=4), intent(in)                 :: maxiteration
        real(kind=8), intent(in)                    :: tolerance
        real(kind=8), intent(in)                    :: omega_sor

        real(kind=8)                                :: rsd_val, res0tol
        real(kind=8), allocatable                   :: tmp(:,:,:)
        integer(kind=4)                             :: i, l, cyc
        integer(kind=4)                             :: ierr
        integer(kind=4), pointer                    :: lv_aggr_max_ptr, lv_aggr_med_ptr, lv_aggr_min_ptr
        type(cart_comm_1d), pointer                 :: comm_max_ptr, comm_med_ptr, comm_min_ptr
        character(len=1)                            :: max, med, min


        ! MPI_Gatherv/MPI_Scatterv used for aggregation
        integer(kind=4)                             :: sizes(0:2), subsizes(0:2), starts(0:2)
        integer(kind=4)                             :: r8size
        integer(kind=4)                             :: nx_aggr, ny_aggr, nz_aggr
        integer(kind=4)                             :: n_part
        integer(kind=MPI_ADDRESS_KIND)              :: extent, lb
        integer(kind=4)                             :: ddtype_temp1
        integer(kind=4), target                     :: ddtype_gatherv_x         !> Derived datatype for grid aggregation
        integer(kind=4), target                     :: ddtype_gatherv_y         !> Derived datatype for grid aggregation
        integer(kind=4), target                     :: ddtype_gatherv_z         !> Derived datatype for grid aggregation
        integer(kind=4), pointer                    :: ddtype_gatherv_max_ptr   !> Pointer to change aggregation order
        integer(kind=4), pointer                    :: ddtype_gatherv_med_ptr   !> Pointer to change aggregation order
        integer(kind=4), pointer                    :: ddtype_gatherv_min_ptr   !> Pointer to change aggregation order
        integer(kind=4), target                     :: ddtype_scatterv_x        !> Derived datatype for grid disaggregation
        integer(kind=4), target                     :: ddtype_scatterv_y        !> Derived datatype for grid disaggregation
        integer(kind=4), target                     :: ddtype_scatterv_z        !> Derived datatype for grid disaggregation
        integer(kind=4), pointer                    :: ddtype_scatterv_max_ptr  !> Pointer to change disaggregation order
        integer(kind=4), pointer                    :: ddtype_scatterv_med_ptr  !> Pointer to change disaggregation order
        integer(kind=4), pointer                    :: ddtype_scatterv_min_ptr  !> Pointer to change disaggregation order
        integer(kind=4), allocatable, dimension(:), target  :: cnt_x, disps_x
        integer(kind=4), allocatable, dimension(:), target  :: cnt_y, disps_y
        integer(kind=4), allocatable, dimension(:), target  :: cnt_z, disps_z
        integer(kind=4), dimension(:), pointer              :: cnt_max_ptr, disps_max_ptr
        integer(kind=4), dimension(:), pointer              :: cnt_med_ptr, disps_med_ptr
        integer(kind=4), dimension(:), pointer              :: cnt_min_ptr, disps_min_ptr

        ! Start of MPI communication setup
        if(lv_aggregation_x .gt. 0) then

            nx_aggr = mg_sdm(lv_aggregation_x)%nx
            ny_aggr = mg_sdm(lv_aggregation_x)%ny
            nz_aggr = mg_sdm(lv_aggregation_x)%nz

            n_part = nx_aggr / comm_1d_x%nprocs

            ! Derivide datatype for grid aggregation using MPI_Gatherv in x-direction
            sizes    = (/nx_aggr+2, ny_aggr+2, nz_aggr+2/)
            subsizes = (/n_part, ny_aggr, nz_aggr/)
            starts   = (/1, 1, 1/)

            call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_REAL8, ddtype_temp1, ierr)
            CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
            lb = 0
            extent = r8size
            call MPI_Type_create_resized(ddtype_temp1, lb, extent, ddtype_gatherv_x, ierr)
            call MPI_Type_commit(ddtype_gatherv_x, ierr)

            ! Derivide datatype for grid aggregation using MPI_Scatterv in x-direction
            sizes    = (/nx_aggr+2, ny_aggr+2, nz_aggr+2/)
            subsizes = (/n_part+2, ny_aggr+2, nz_aggr+2/)
            starts   = (/0, 0, 0/)

            call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_REAL8, ddtype_temp1, ierr)
            CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
            lb = 0
            extent = r8size
            call MPI_Type_create_resized(ddtype_temp1, lb, extent, ddtype_scatterv_x, ierr)
            call MPI_Type_commit(ddtype_scatterv_x, ierr)

            ! DDT counts and displacements
            allocate(cnt_x(0:comm_1d_x%nprocs-1))
            allocate(disps_x(0:comm_1d_x%nprocs-1))

            do i = 0, comm_1d_x%nprocs-1
                cnt_x(i) = 1
                disps_x(i) = n_part * i
            enddo

        endif


        if(lv_aggregation_y .gt. 0) then

            nx_aggr = mg_sdm(lv_aggregation_y)%nx
            ny_aggr = mg_sdm(lv_aggregation_y)%ny
            nz_aggr = mg_sdm(lv_aggregation_y)%nz

            n_part = ny_aggr / comm_1d_y%nprocs

            ! Derivide datatype for grid aggregation using MPI_Gatherv in y-direction
            sizes    = (/nx_aggr+2, ny_aggr+2, nz_aggr+2/)
            subsizes = (/nx_aggr, n_part, nz_aggr/)
            starts   = (/1, 1, 1/)

            call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_REAL8, ddtype_temp1, ierr)
            CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
            lb = 0
            extent = r8size
            call MPI_Type_create_resized(ddtype_temp1, lb, extent, ddtype_gatherv_y, ierr)
            call MPI_Type_commit(ddtype_gatherv_y, ierr)

            ! Derivide datatype for grid aggregation using MPI_Scatterv in y-direction
            sizes    = (/nx_aggr+2, ny_aggr+2, nz_aggr+2/)
            subsizes = (/nx_aggr+2, n_part+2, nz_aggr+2/)
            starts   = (/0, 0, 0/)

            call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_REAL8, ddtype_temp1, ierr)
            CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
            lb = 0
            extent = r8size
            call MPI_Type_create_resized(ddtype_temp1, lb, extent, ddtype_scatterv_y, ierr)
            call MPI_Type_commit(ddtype_scatterv_y, ierr)

            ! DDT counts and displacements
            allocate(cnt_y(0:comm_1d_y%nprocs-1))
            allocate(disps_y(0:comm_1d_y%nprocs-1))

            do i = 0, comm_1d_y%nprocs-1
                cnt_y(i) = 1
                disps_y(i) = (nx_aggr+2) * n_part * i
            enddo
        endif

        if(lv_aggregation_z .gt. 0) then

            nx_aggr = mg_sdm(lv_aggregation_z)%nx
            ny_aggr = mg_sdm(lv_aggregation_z)%ny
            nz_aggr = mg_sdm(lv_aggregation_z)%nz

            n_part = nz_aggr / comm_1d_z%nprocs

            ! Derivide datatype for grid aggregation using MPI_Gatherv in z-direction
            sizes    = (/nx_aggr+2, ny_aggr+2, nz_aggr+2/)
            subsizes = (/nx_aggr, ny_aggr, n_part/)
            starts   = (/1, 1, 1/)

            call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_REAL8, ddtype_temp1, ierr)
            CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
            lb = 0
            extent = r8size
            call MPI_Type_create_resized(ddtype_temp1, lb, extent, ddtype_gatherv_z, ierr)
            call MPI_Type_commit(ddtype_gatherv_z, ierr)

            ! Derivide datatype for grid aggregation using MPI_Scatterv in z-direction
            sizes    = (/nx_aggr+2, ny_aggr+2, nz_aggr+2/)
            subsizes = (/nx_aggr+2, ny_aggr+2, n_part+2/)
            starts   = (/0, 0, 0/)

            call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_REAL8, ddtype_temp1, ierr)
            CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
            lb = 0
            extent = r8size
            call MPI_Type_create_resized(ddtype_temp1, lb, extent, ddtype_scatterv_z, ierr)
            call MPI_Type_commit(ddtype_scatterv_z, ierr)

            ! DDT counts and displacements
            allocate(cnt_z(0:comm_1d_z%nprocs-1))
            allocate(disps_z(0:comm_1d_z%nprocs-1))

            do i = 0, comm_1d_z%nprocs-1
                cnt_z(i) = 1
                disps_z(i) = (nx_aggr+2) * (ny_aggr+2) * n_part * i
            enddo
        endif

        ! Change the order of aggregation in each direction according to aggregation level
        if(lv_aggregation_x.gt.lv_aggregation_y) then
            if(lv_aggregation_x.gt.lv_aggregation_z) then
                if(lv_aggregation_y.gt.lv_aggregation_z) then
                    ! x -> y -> z
                    lv_aggr_max_ptr => lv_aggregation_x
                    lv_aggr_med_ptr => lv_aggregation_y
                    lv_aggr_min_ptr => lv_aggregation_z
                    comm_max_ptr => comm_1d_x
                    comm_med_ptr => comm_1d_y
                    comm_min_ptr => comm_1d_z

                    ddtype_gatherv_max_ptr => ddtype_gatherv_x
                    ddtype_gatherv_med_ptr => ddtype_gatherv_y
                    ddtype_gatherv_min_ptr => ddtype_gatherv_z

                    ddtype_scatterv_max_ptr => ddtype_scatterv_x
                    ddtype_scatterv_med_ptr => ddtype_scatterv_y
                    ddtype_scatterv_min_ptr => ddtype_scatterv_z

                    cnt_max_ptr => cnt_x
                    cnt_med_ptr => cnt_y
                    cnt_min_ptr => cnt_z
                    disps_max_ptr => disps_x
                    disps_med_ptr => disps_y
                    disps_min_ptr => disps_z

                    max = 'x'
                    med = 'y'
                    min = 'z'
                else
                    ! x -> z -> y
                    lv_aggr_max_ptr => lv_aggregation_x
                    lv_aggr_med_ptr => lv_aggregation_z
                    lv_aggr_min_ptr => lv_aggregation_y
                    comm_max_ptr => comm_1d_x
                    comm_med_ptr => comm_1d_z
                    comm_min_ptr => comm_1d_y

                    ddtype_gatherv_max_ptr => ddtype_gatherv_x
                    ddtype_gatherv_med_ptr => ddtype_gatherv_z
                    ddtype_gatherv_min_ptr => ddtype_gatherv_y

                    ddtype_scatterv_max_ptr => ddtype_scatterv_x
                    ddtype_scatterv_med_ptr => ddtype_scatterv_z
                    ddtype_scatterv_min_ptr => ddtype_scatterv_y

                    cnt_max_ptr => cnt_x
                    cnt_med_ptr => cnt_z
                    cnt_min_ptr => cnt_y
                    disps_max_ptr => disps_x
                    disps_med_ptr => disps_z
                    disps_min_ptr => disps_y

                    max = 'x'
                    med = 'z'
                    min = 'y'
                endif
            else
                ! z -> x -> y
                lv_aggr_max_ptr => lv_aggregation_z
                lv_aggr_med_ptr => lv_aggregation_x
                lv_aggr_min_ptr => lv_aggregation_y
                comm_max_ptr => comm_1d_z
                comm_med_ptr => comm_1d_x
                comm_min_ptr => comm_1d_y

                ddtype_gatherv_max_ptr => ddtype_gatherv_z
                ddtype_gatherv_med_ptr => ddtype_gatherv_x
                ddtype_gatherv_min_ptr => ddtype_gatherv_y

                ddtype_scatterv_max_ptr => ddtype_scatterv_z
                ddtype_scatterv_med_ptr => ddtype_scatterv_x
                ddtype_scatterv_min_ptr => ddtype_scatterv_y

                cnt_max_ptr => cnt_z
                cnt_med_ptr => cnt_x
                cnt_min_ptr => cnt_y
                disps_max_ptr => disps_z
                disps_med_ptr => disps_x
                disps_min_ptr => disps_y

                max = 'z'
                med = 'x'
                min = 'y'
            endif
        else
            if(lv_aggregation_y.gt.lv_aggregation_z) then
                if(lv_aggregation_x.gt.lv_aggregation_z) then
                    ! y -> x -> z
                    lv_aggr_max_ptr => lv_aggregation_y
                    lv_aggr_med_ptr => lv_aggregation_x
                    lv_aggr_min_ptr => lv_aggregation_z
                    comm_max_ptr => comm_1d_y
                    comm_med_ptr => comm_1d_x
                    comm_min_ptr => comm_1d_z

                    ddtype_gatherv_max_ptr => ddtype_gatherv_y
                    ddtype_gatherv_med_ptr => ddtype_gatherv_x
                    ddtype_gatherv_min_ptr => ddtype_gatherv_z

                    ddtype_scatterv_max_ptr => ddtype_scatterv_y
                    ddtype_scatterv_med_ptr => ddtype_scatterv_x
                    ddtype_scatterv_min_ptr => ddtype_scatterv_z

                    cnt_max_ptr => cnt_y
                    cnt_med_ptr => cnt_x
                    cnt_min_ptr => cnt_z
                    disps_max_ptr => disps_y
                    disps_med_ptr => disps_x
                    disps_min_ptr => disps_z

                    max = 'y'
                    med = 'x'
                    min = 'z'
                else
                    ! y -> z -> x
                    lv_aggr_max_ptr => lv_aggregation_y
                    lv_aggr_med_ptr => lv_aggregation_z
                    lv_aggr_min_ptr => lv_aggregation_x
                    comm_max_ptr => comm_1d_y
                    comm_med_ptr => comm_1d_z
                    comm_min_ptr => comm_1d_x

                    ddtype_gatherv_max_ptr => ddtype_gatherv_y
                    ddtype_gatherv_med_ptr => ddtype_gatherv_z
                    ddtype_gatherv_min_ptr => ddtype_gatherv_x

                    ddtype_scatterv_max_ptr => ddtype_scatterv_y
                    ddtype_scatterv_med_ptr => ddtype_scatterv_z
                    ddtype_scatterv_min_ptr => ddtype_scatterv_x

                    cnt_max_ptr => cnt_y
                    cnt_med_ptr => cnt_z
                    cnt_min_ptr => cnt_x
                    disps_max_ptr => disps_y
                    disps_med_ptr => disps_z
                    disps_min_ptr => disps_x

                    max = 'y'
                    med = 'z'
                    min = 'x'
                endif
            else
                ! z -> y -> x
                lv_aggr_max_ptr => lv_aggregation_z
                lv_aggr_med_ptr => lv_aggregation_y
                lv_aggr_min_ptr => lv_aggregation_x
                comm_max_ptr => comm_1d_z
                comm_med_ptr => comm_1d_y
                comm_min_ptr => comm_1d_x

                ddtype_gatherv_max_ptr => ddtype_gatherv_z
                ddtype_gatherv_med_ptr => ddtype_gatherv_y
                ddtype_gatherv_min_ptr => ddtype_gatherv_x

                ddtype_scatterv_max_ptr => ddtype_scatterv_z
                ddtype_scatterv_med_ptr => ddtype_scatterv_y
                ddtype_scatterv_min_ptr => ddtype_scatterv_x

                cnt_max_ptr => cnt_z
                cnt_med_ptr => cnt_y
                cnt_min_ptr => cnt_x
                disps_max_ptr => disps_z
                disps_med_ptr => disps_y
                disps_min_ptr => disps_x

                max = 'z'
                med = 'y'
                min = 'x'
            endif
        endif
        ! End of MPI communication setup

        call multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm%is_aggregated)
        call vv_dot_3d_matrix(res0tol, rsd, rsd, sdm%nx, sdm%ny, sdm%nz, sdm%is_aggregated)

        ! V-cycle Iterations
        do cyc = 1, n_vcycles

            ! Smooting and restriction until aggregation level min -> med -> max
            call rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm%is_aggregated)
            call multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm%is_aggregated)
            call multigrid_restriction(mg_sdm(1)%b, rsd, mg_sdm(1), sdm, 0)

            do l = 1, max0(1, lv_aggr_min_ptr)-1
                mg_sdm(l)%x = 0.0d0
                call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                call multigrid_residual(mg_sdm(l)%r, mg_a_poisson(l), mg_sdm(l)%x, mg_sdm(l)%b, mg_sdm(l), mg_sdm(l)%is_aggregated)
                call multigrid_restriction(mg_sdm(l+1)%b, mg_sdm(l)%r, mg_sdm(l+1), mg_sdm(l),l)
            enddo

            ! Grid aggregation in the direction where aggregation level is minimum
            if(lv_aggr_min_ptr.gt.0) then
                allocate(tmp(0:size(mg_sdm(lv_aggr_min_ptr)%b, 1)-1,0:size(mg_sdm(lv_aggr_min_ptr)%b,2)-1,0:size(mg_sdm(lv_aggr_min_ptr)%b,3)-1))
                tmp = 0.0d0
                call MPI_Gatherv(mg_sdm(lv_aggr_min_ptr)%b(0,0,0), 1, ddtype_gatherv_min_ptr, tmp(0,0,0), cnt_min_ptr, disps_min_ptr, ddtype_gatherv_min_ptr, 0, comm_min_ptr%mpi_comm, ierr)
                mg_sdm(lv_aggr_min_ptr)%b = tmp
                deallocate(tmp)
            endif

            ! Rank0 in the communicator of minimum aggregation level, take parts in MG steps
            if(comm_min_ptr%myrank.eq.0) then

                do l = max0(1, lv_aggr_min_ptr), lv_aggr_med_ptr-1
                    mg_sdm(l)%x = 0.0d0
                    call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                    call multigrid_residual(mg_sdm(l)%r, mg_a_poisson(l), mg_sdm(l)%x, mg_sdm(l)%b, mg_sdm(l), mg_sdm(l)%is_aggregated)
                    call multigrid_restriction(mg_sdm(l+1)%b, mg_sdm(l)%r, mg_sdm(l+1), mg_sdm(l),l)
                enddo
    
                ! Grid aggregation in the direction where aggregation level is median
                if(lv_aggr_med_ptr.gt.0) then
                    allocate(tmp(0:size(mg_sdm(lv_aggr_med_ptr)%b, 1)-1,0:size(mg_sdm(lv_aggr_med_ptr)%b,2)-1,0:size(mg_sdm(lv_aggr_med_ptr)%b,3)-1))
                    tmp = 0.0d0
                    call MPI_Gatherv(mg_sdm(lv_aggr_med_ptr)%b(0,0,0), 1, ddtype_gatherv_med_ptr, tmp(0,0,0), cnt_med_ptr, disps_med_ptr, ddtype_gatherv_med_ptr, 0, comm_med_ptr%mpi_comm, ierr)
                    mg_sdm(lv_aggr_med_ptr)%b = tmp
                    deallocate(tmp)
                endif

                ! Rank0 in the communicator of median aggregation level, take parts in MG steps
                if(comm_med_ptr%myrank.eq.0) then

                    do l = max0(1, lv_aggr_med_ptr), lv_aggr_max_ptr-1
                        mg_sdm(l)%x = 0.0d0
                        call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                        call multigrid_residual(mg_sdm(l)%r, mg_a_poisson(l), mg_sdm(l)%x, mg_sdm(l)%b, mg_sdm(l), mg_sdm(l)%is_aggregated)
                        call multigrid_restriction(mg_sdm(l+1)%b, mg_sdm(l)%r, mg_sdm(l+1), mg_sdm(l),l)
                    enddo
        
                    ! Grid aggregation in the direction where aggregation level is maximum
                    if(lv_aggr_max_ptr.gt.0) then
                        allocate(tmp(0:size(mg_sdm(lv_aggr_max_ptr)%b, 1)-1,0:size(mg_sdm(lv_aggr_max_ptr)%b,2)-1,0:size(mg_sdm(lv_aggr_max_ptr)%b,3)-1))
                        tmp = 0.0d0
                        call MPI_Gatherv(mg_sdm(lv_aggr_max_ptr)%b(0,0,0), 1, ddtype_gatherv_max_ptr, tmp(0,0,0), cnt_max_ptr, disps_max_ptr, ddtype_gatherv_max_ptr, 0, comm_max_ptr%mpi_comm, ierr)
                        mg_sdm(lv_aggr_max_ptr)%b = tmp
                        deallocate(tmp)
                    endif
        
                    ! Rank0 in the communicator of maximum aggregation level, take parts in MG steps
                    if(comm_max_ptr%myrank.eq.0) then
                        do l = max0(1, lv_aggr_max_ptr), n_levels-1
                            mg_sdm(l)%x = 0.0d0
                            call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                            call multigrid_residual(mg_sdm(l)%r, mg_a_poisson(l), mg_sdm(l)%x, mg_sdm(l)%b, mg_sdm(l), mg_sdm(l)%is_aggregated)
                            call multigrid_restriction(mg_sdm(l+1)%b, mg_sdm(l)%r, mg_sdm(l+1), mg_sdm(l),l)
                        enddo

                        ! Solve at the coarset level only for rank0 in the communicator of min., med., and max..
                        call multigrid_solve_coarset_level( mg_sdm(n_levels)%x, &
                                                                mg_a_poisson(n_levels), &
                                                                mg_sdm(n_levels)%b, &
                                                                mg_sdm(n_levels), &
                                                                1000, tolerance, omega_sor, &
                                                                mg_sdm(n_levels)%is_aggregated)

                        call multigrid_residual(mg_sdm(n_levels)%r, mg_a_poisson(n_levels), mg_sdm(n_levels)%x, mg_sdm(n_levels)%b, mg_sdm(n_levels), mg_sdm(n_levels)%is_aggregated)
                        print '(a,e18.10,a,e18.10)','[MG] Solution in the coarset level : x(1,1,1) = ',mg_sdm(n_levels)%x(1,1,1),', residue = ',mg_sdm(n_levels)%r(1,1,1)
#ifdef DEBUG_COARSEST
                        call multigrid_debug_print_coarsest_level_solution(cyc, mg_sdm(n_levels))
#endif
                        ! Prolongation and smoothing : Ghostcell update is required for prolongation
                        do l = n_levels-1, max0(1, lv_aggr_max_ptr), -1
                            call multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm(l)%r, mg_sdm(l+1)%x, mg_sdm(l), mg_sdm(l+1),l)
                            mg_sdm(l)%x = mg_sdm(l)%x + mg_sdm(l)%r
                            call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                        enddo
                    endif

                    ! Grid disaggregation in the direction where aggregation level is maximum
                    if(lv_aggr_max_ptr.gt.0) then
                        allocate(tmp(0:size(mg_sdm(lv_aggr_max_ptr)%x, 1)-1,0:size(mg_sdm(lv_aggr_max_ptr)%x,2)-1,0:size(mg_sdm(lv_aggr_max_ptr)%x,3)-1))
                        tmp = 0.0d0
                        call MPI_Scatterv(mg_sdm(lv_aggr_max_ptr)%x(0,0,0), cnt_max_ptr, disps_max_ptr, ddtype_scatterv_max_ptr, tmp, 1, ddtype_scatterv_max_ptr, 0, comm_max_ptr%mpi_comm, ierr)
                        mg_sdm(lv_aggr_max_ptr)%x = tmp
                        deallocate(tmp)
                    endif

                    ! Prolongation and smoothing : Ghostcell update is required for prolongation
                    do l = lv_aggr_max_ptr-1, max0(1, lv_aggr_med_ptr), -1
                        call multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm(l)%r, mg_sdm(l+1)%x, mg_sdm(l), mg_sdm(l+1),l)
                        mg_sdm(l)%x = mg_sdm(l)%x + mg_sdm(l)%r
                        call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                        call geometry_halocell_update_selectively(mg_sdm(l)%x, mg_sdm(l), mg_sdm(l)%is_aggregated)
                    enddo
                endif

                ! Grid disaggregation in the direction where aggregation level is median
                if(lv_aggr_med_ptr.gt.0) then
                    allocate(tmp(0:size(mg_sdm(lv_aggr_med_ptr)%x, 1)-1,0:size(mg_sdm(lv_aggr_med_ptr)%x,2)-1,0:size(mg_sdm(lv_aggr_med_ptr)%x,3)-1))
                    tmp = 0.0d0
                    call MPI_Scatterv(mg_sdm(lv_aggr_med_ptr)%x(0,0,0), cnt_med_ptr, disps_med_ptr, ddtype_scatterv_med_ptr, tmp, 1, ddtype_scatterv_med_ptr, 0, comm_med_ptr%mpi_comm, ierr)
                    mg_sdm(lv_aggr_med_ptr)%x = tmp
                    deallocate(tmp)
                endif

                ! Prolongation and smoothing : Ghostcell update is required for prolongation
                do l = lv_aggr_med_ptr-1, max0(1, lv_aggr_min_ptr), -1
                    call multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm(l)%r, mg_sdm(l+1)%x, mg_sdm(l), mg_sdm(l+1),l)
                    mg_sdm(l)%x = mg_sdm(l)%x + mg_sdm(l)%r
                    call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                    call geometry_halocell_update_selectively(mg_sdm(l)%x, mg_sdm(l), mg_sdm(l)%is_aggregated)
                enddo
            endif

            ! Grid disaggregation in the direction where aggregation level is minimum
            if(lv_aggr_min_ptr.gt.0) then
                allocate(tmp(0:size(mg_sdm(lv_aggr_min_ptr)%x, 1)-1,0:size(mg_sdm(lv_aggr_min_ptr)%x,2)-1,0:size(mg_sdm(lv_aggr_min_ptr)%x,3)-1))
                tmp = 0.0d0
                call MPI_Scatterv(mg_sdm(lv_aggr_min_ptr)%x(0,0,0), cnt_min_ptr, disps_min_ptr, ddtype_scatterv_min_ptr, tmp, 1, ddtype_scatterv_min_ptr, 0, comm_min_ptr%mpi_comm, ierr)
                mg_sdm(lv_aggr_min_ptr)%x = tmp
                deallocate(tmp)
            endif

            ! Prolongation and smoothing : Ghostcell update is required for prolongation
            do l = lv_aggr_min_ptr-1, 1, -1
                call multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm(l)%r, mg_sdm(l+1)%x, mg_sdm(l), mg_sdm(l+1),l)
                mg_sdm(l)%x = mg_sdm(l)%x + mg_sdm(l)%r
                call rbgs_iterator_poisson_matrix(mg_sdm(l)%x, mg_a_poisson(l), mg_sdm(l)%b, mg_sdm(l), maxiteration, omega_sor, mg_sdm(l)%is_aggregated)
                call geometry_halocell_update_selectively(mg_sdm(l)%x, mg_sdm(l), mg_sdm(l)%is_aggregated)
            enddo

            call multigrid_prolongation_linear_on_nonuniform_grid(rsd, mg_sdm(1)%x, sdm, mg_sdm(1), 0)
            sol = sol + rsd
            call rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm%is_aggregated)

            ! Convergence check
            call multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm%is_aggregated)
            call vv_dot_3d_matrix(rsd_val, rsd, rsd, sdm%nx, sdm%ny, sdm%nz, sdm%is_aggregated)
            if(myrank.eq.0) print '(a,i4,a,3(e15.7,x))', '[MG] cycle = ',cyc,', Error:',sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol)
            if (sqrt(rsd_val/res0tol) .LT. tolerance) exit
        enddo

        if(lv_aggregation_x .gt. 0) then
            deallocate(cnt_x  )
            deallocate(disps_x)
        endif
        if(lv_aggregation_y .gt. 0) then
            deallocate(cnt_y  )
            deallocate(disps_y)
        endif
        if(lv_aggregation_z .gt. 0) then
            deallocate(cnt_z  )
            deallocate(disps_z)
        endif

    end subroutine multigrid_CGPSA_vcycle_solver

end module multigrid