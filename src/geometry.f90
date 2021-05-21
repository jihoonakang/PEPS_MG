!======================================================================================================================
!> @file        mpi_subdomain.f90
!> @brief       This file contains a module that defines a geometry for the example problem of PEPS_MG.
!> @details     In the geometry module, global domain and partitioned domains(subdomain) are defined and their grid 
!>              information is allocated. Derived Datatype for ghostcell communication is also defined.
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
!> @brief       Module for building ubdomains from the physical global domain.
!> @details     This module containes domain and grid information of global domain and subdomains. 
!>              Derived Datatype for ghostcell communication between the subdomains is included.
!>
module geometry

    implicit none

    !> @brief   Global domain information
    !> @details It contains the geometry of global domain in three-dimensional Cartesian cooridnates. 
    !>          Dimensions of grids and meshes are included. Cell-centered grid system is employed.
    type, public    :: domain

        !> @{ Number of grids in the global domain
        integer(kind=4)     :: nx, ny, nz
        !> @}
        !> @{ Physical length of the global domain
        real(kind=8)        :: lx, ly, lz
        !> @}
        !> @{ Starting coordinates of the global domain. (x0, y0, z0)
        real(kind=8)        :: ox, oy, oz
        !> @}
        !> @{ Mesh (cell) size in each dimension
        real(kind=8), allocatable, dimension(:)         :: dxm, dym, dzm
        !> @}
        !> @{ Space between grid points
        real(kind=8), allocatable, dimension(:)         :: dxg, dyg, dzg
        !> @}
        !> @{ Coordinates of grid pointss
        real(kind=8), allocatable, dimension(:)         :: xg, yg, zg
        !> @}
        !> @{ Periodic boundary condition: .true. if periodic, .false. if not
        logical                                         :: is_periodic(0:2)
        !> @}
    end type domain

    !> @brief   Partitioned subdomain information
    !> @details It contains the geometry of subdomain in three-dimensional Cartesian cooridnates. 
    !>          Dimensions of grids and meshes are included. Cell-centered grid system is employed.
    !>          Grid variables, solution, RHS and residual, is defined in this type variable.
    type, public    :: subdomain
        !> @{ Grid numbers in the partitioned subdomain
        integer(kind=4)     :: nx, ny, nz
        !> @}
        !> @{ Physical length of the subdomain
        real(kind=8)        :: lx, ly, lz
        !> @}
        !> @{ Starting coordinates of the subdomain domain.
        real(kind=8)        :: ox, oy, oz
        !> @}
        !> @{ Mesh (cell) size in each dimension of subdomain
        real(kind=8), allocatable, dimension(:)         :: dxm, dym, dzm
        !> @}
        !> @{ Space between grid points of subdomain
        real(kind=8), allocatable, dimension(:)         :: dxg, dyg, dzg
        !> @}
        !> @{ Coordinates of grid pointss of subdomain
        real(kind=8), allocatable, dimension(:)         :: xg, yg, zg
        !> @}
        !> @{ Periodic boundary condition: .true. if periodic, .false. if not
        logical                                         :: is_periodic(0:2)
        !> @}
        real(kind=8), allocatable, dimension(:,:,:)     :: x    !> Solution variable in subdomain
        real(kind=8), allocatable, dimension(:,:,:)     :: b    !> RHS variable in subdomain
        real(kind=8), allocatable, dimension(:,:,:)     :: r    !> Residual variable in subdomain

        !> @{ Indicate whether subdomain is aggregated in each direction: .true. if aggregated, .false. if not
        logical                 :: is_aggregated(0:2)
        !> @}
        !> @{ Indices of first and last grids of subdomain in x-direction, with regards to global domain
        integer(kind=4), public :: ista, iend
        !> @}
        !> @{ Indices of first and last grids of subdomain in y-direction, with regards to global domain
        integer(kind=4), public :: jsta, jend
        !> @}
        !> @{ Indices of first and last grids of subdomain in z-direction, with regards to global domain
        integer(kind=4), public :: ksta, kend
        !> @}
        !> @{ Derived datatypes for ghostcell communication in x-direction
        integer(kind=4), public :: ddt_yz_plane_x0, ddt_yz_plane_x1, ddt_yz_plane_xn, ddt_yz_plane_xn1
        !> @}
        !> @{ Derived datatypes for ghostcell communication in y-direction
        integer(kind=4), public :: ddt_xz_plane_y0, ddt_xz_plane_y1, ddt_xz_plane_yn, ddt_xz_plane_yn1
        !> @}
        !> @{ Derived datatypes for ghostcell communication in z-direction
        integer(kind=4), public :: ddt_xy_plane_z0, ddt_xy_plane_z1, ddt_xy_plane_zn, ddt_xy_plane_zn1
        !> @}
        !> @{ Derived datatypes for inner domain without ghostcell
        integer(kind=4), public :: ddt_inner_domain
        !> @}
        !> @{ Indicates whether the subdomain is x-boundary: .true. if at boundary, .false. if not
        logical                 :: is_x0_boundary, is_x1_boundary
        !> @}
        !> @{ Indicates whether the subdomain is y-boundary: .true. if at boundary, .false. if not
        logical                 :: is_y0_boundary, is_y1_boundary
        !> @}
        !> @{ Indicates whether the subdomain is z-boundary: .true. if at boundary, .false. if not
        logical                 :: is_z0_boundary, is_z1_boundary
        !> @}
        end type subdomain

    contains

    !>
    !> @brief       Prepare the global domain and grid information.
    !> @param       gdm         Global domain
    !> @param       nx          Number of grid in x-direction
    !> @param       ny          Number of grid in y-direction
    !> @param       nz          Number of grid in z-direction
    !> @param       ox          Starting coordinate of global domain in x-direction
    !> @param       oy          Starting coordinate of global domain in y-direction
    !> @param       oz          Starting coordinate of  global domain in z-direction
    !> @param       lx          Physical length of global domain in x-direction
    !> @param       ly          Physical length of global domain in y-direction
    !> @param       lz          Physical length of global domain in z-direction
    !> @param       ax          Mesh stretch ratio in x-direction
    !> @param       ay          Mesh stretch ratio in y-direction
    !> @param       az          Mesh stretch ratio in z-direction
    !> @param       period      Periodicity of global domain
    !>
    subroutine geometry_domain_create(gdm, nx, ny, nz, ox, oy, oz, lx, ly, lz, ax, ay, az, period)
        
        implicit none

        type(domain), intent(inout)     :: gdm
        integer(kind=4), intent(in)     :: nx, ny, nz
        real(kind=8),    intent(in)     :: ox, oy, oz
        real(kind=8),    intent(in)     :: lx, ly, lz
        real(kind=8),    intent(in)     :: ax, ay, az
        logical,         intent(in)     :: period(0:2)

        integer(kind=4) :: i

        ! Assign basic parameters of global domain
        gdm%nx = nx
        gdm%ny = ny
        gdm%nz = nz
        gdm%ox = ox
        gdm%oy = oy
        gdm%oz = oz
        gdm%lx = lx
        gdm%ly = ly
        gdm%lz = lz
        gdm%is_periodic(0:2) = period(0:2)

        allocate( gdm%dxm(0:gdm%nx+1) ) 
        allocate( gdm%dym(0:gdm%ny+1) ) 
        allocate( gdm%dzm(0:gdm%nz+1) ) 
    
        allocate( gdm%dxg(0:gdm%nx+1) ) 
        allocate( gdm%dyg(0:gdm%ny+1) ) 
        allocate( gdm%dzg(0:gdm%nz+1) ) 
    
        allocate( gdm%xg(0:gdm%nx+1) ) 
        allocate( gdm%yg(0:gdm%ny+1) ) 
        allocate( gdm%zg(0:gdm%nz+1) ) 
    
        ! Setup mesh size, grid spacing and grid coordinates in x-direction
        if(ax.eq.1.0d0) then
            ! Uniform grid
            gdm%dxm(0) = 0.0d0
            gdm%dxg(0) = 0.0d0
            gdm%xg(0) = gdm%ox
            do i = 1, gdm%nx
                gdm%dxm(i) = gdm%lx / dble(gdm%nx)
                gdm%dxg(i) = 0.5d0 * (gdm%dxm(i) + gdm%dxm(i-1))
                gdm%xg(i) = gdm%xg(i-1) + gdm%dxg(i)
            enddo
            gdm%dxm(gdm%nx+1) = 0.0d0
            gdm%dxg(gdm%nx+1) = 0.5d0 * (gdm%dxm(gdm%nx+1) + gdm%dxm(gdm%nx))
            gdm%xg(gdm%nx+1) = gdm%xg(gdm%nx) + gdm%dxg(gdm%nx+1)
        else
            ! Non-uniform grid with stretch ratio of ax
            gdm%dxm(0) = 0.0d0
            gdm%dxm(1) = 0.5d0 * lx * (ax - 1.0d0) / (ax**(gdm%nx/2) - 1.0d0)
            do i = 2, gdm%nx / 2
                gdm%dxm(i) = gdm%dxm(i-1) * ax
                gdm%dxm(gdm%nx-i+1) = gdm%dxm(i)
            enddo
            gdm%dxm(gdm%nx) = gdm%dxm(1)
            gdm%dxm(gdm%nx+1) = 0.0d0

            gdm%dxg(0) = 0.0d0
            gdm%xg(0) = gdm%ox
            do i = 1, gdm%nx
                gdm%dxg(i) = 0.5d0 * (gdm%dxm(i) + gdm%dxm(i-1))
                gdm%xg(i) = gdm%xg(i-1) + gdm%dxg(i)
            enddo
            gdm%dxg(gdm%nx+1) = 0.5d0 * (gdm%dxm(gdm%nx+1) + gdm%dxm(gdm%nx))
            gdm%xg(gdm%nx+1)  = gdm%xg(gdm%nx) + gdm%dxg(gdm%nx+1)
        endif

        ! Setup mesh size, grid spacing and grid coordinates in y-direction
        if(ay.eq.1.0d0) then
            ! Uniform grid
            gdm%dym(0) = 0.0d0
            gdm%dyg(0) = 0.0d0
            gdm%yg(0) = gdm%oy
            do i = 1, gdm%ny
                gdm%dym(i) = gdm%ly / dble(gdm%ny)
                gdm%dyg(i) = 0.5d0 * (gdm%dym(i) + gdm%dym(i-1))
                gdm%yg(i) = gdm%yg(i-1) + gdm%dyg(i)
            enddo
            gdm%dym(gdm%ny+1) = 0.0d0
            gdm%dyg(gdm%ny+1) = 0.5d0 * (gdm%dym(gdm%ny+1) + gdm%dym(gdm%ny))
            gdm%yg(gdm%ny+1) = gdm%yg(gdm%ny) + gdm%dyg(gdm%ny+1)
        else
            ! Non-uniform grid with stretch ratio of ay
            gdm%dym(0) = 0.0d0
            gdm%dym(1) = 0.5d0 * ly * (ay - 1.0d0) / (ay**(gdm%ny/2) - 1.0d0)
            do i = 2, gdm%ny / 2
                gdm%dym(i) = gdm%dym(i-1) * ay
                gdm%dym(gdm%ny-i+1) = gdm%dym(i)
            enddo
            gdm%dym(gdm%ny) = gdm%dym(1)
            gdm%dym(gdm%ny+1) = 0.0d0

            gdm%dyg(0) = 0.0d0
            gdm%yg(0) = gdm%oy
            do i = 1, gdm%ny
                gdm%dyg(i) = 0.5d0 * (gdm%dym(i) + gdm%dym(i-1))
                gdm%yg(i) = gdm%yg(i-1) + gdm%dyg(i)
            enddo
            gdm%dyg(gdm%ny+1) = 0.5d0 * (gdm%dym(gdm%ny+1) + gdm%dym(gdm%ny))
            gdm%yg(gdm%ny+1)  = gdm%yg(gdm%ny) + gdm%dyg(gdm%ny+1)
        endif

        ! Setup mesh size, grid spacing and grid coordinates in z-direction
        if(az.eq.1.0d0) then
            ! Uniform grid
            gdm%dzm(0) = 0.0d0
            gdm%dzg(0) = 0.0d0
            gdm%zg(0) = gdm%oz
            do i = 1, gdm%nz
                gdm%dzm(i) = gdm%lz / dble(gdm%nz)
                gdm%dzg(i) = 0.5d0 * (gdm%dzm(i) + gdm%dzm(i-1))
                gdm%zg(i) = gdm%zg(i-1) + gdm%dzg(i)
            enddo
            gdm%dzm(gdm%nz+1) = 0.0d0
            gdm%dzg(gdm%nz+1) = 0.5d0 * (gdm%dzm(gdm%nz+1) + gdm%dzm(gdm%nz))
            gdm%zg(gdm%nz+1) = gdm%zg(gdm%nz) + gdm%dzg(gdm%nz+1)
        else
            ! Non-uniform grid with stretch ratio of az
            gdm%dzm(0) = 0.0d0
            gdm%dzm(1) = 0.5d0 * lz * (az - 1.0d0) / (az**(gdm%nz/2) - 1.0d0)
            do i = 2, gdm%nz / 2
                gdm%dzm(i) = gdm%dzm(i-1) * az
                gdm%dzm(gdm%nz-i+1) = gdm%dzm(i)
            enddo
            gdm%dzm(gdm%nz) = gdm%dzm(1)
            gdm%dzm(gdm%nz+1) = 0.0d0

            gdm%dzg(0) = 0.0d0
            gdm%zg(0) = gdm%oz
            do i = 1, gdm%nz
                gdm%dzg(i) = 0.5d0 * (gdm%dzm(i) + gdm%dzm(i-1))
                gdm%zg(i) = gdm%zg(i-1) + gdm%dzg(i)
            enddo
            gdm%dzg(gdm%nz+1) = 0.5d0 * (gdm%dzm(gdm%nz+1) + gdm%dzm(gdm%nz))
            gdm%zg(gdm%nz+1)  = gdm%zg(gdm%nz) + gdm%dzg(gdm%nz+1)
        endif

    end subroutine geometry_domain_create

    !>
    !> @brief       Free the variables of global domain
    !> @param       gdm         Global domain
    !>
    subroutine geometry_domain_destroy(gdm)

        implicit none

        type(domain), intent(inout)     :: gdm

        deallocate(gdm%dxm, gdm%dym, gdm%dzm)
        deallocate(gdm%dxg, gdm%dyg, gdm%dzg)
        deallocate(gdm%xg, gdm%yg, gdm%zg)
    
    end subroutine geometry_domain_destroy

    !>
    !> @brief       Prepare the partitioned subdomain from global domain information.
    !> @param       gdm         Global domain
    !> @param       sdm         Partitioned subdomain
    !>
    subroutine geometry_subdomain_create(sdm, gdm)

        use mpi
        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z
        implicit none

        type(subdomain), intent(inout)  :: sdm
        type(domain), intent(in)        :: gdm

        ! Periodicity of subdomain is equal to that of global domain
        sdm%is_periodic(0:2)    = gdm%is_periodic(0:2)

        ! Domain is not partitioned and aggregated when the number of process is equal to 1
        ! In x-direction
        if(comm_1d_x%nprocs.eq.1) then
            sdm%is_aggregated(0)  = .true.
        else
            sdm%is_aggregated(0)  = .false.
        endif

        ! In y-direction
        if(comm_1d_y%nprocs.eq.1) then
            sdm%is_aggregated(1) = .true.
        else
            sdm%is_aggregated(1) = .false.
        endif

        ! In z-direction
        if(comm_1d_z%nprocs.eq.1) then
            sdm%is_aggregated(2) = .true.
        else
            sdm%is_aggregated(2) = .false.
        endif
        
        ! Divide global domain into the partitioned subdomain
        ! 1. Determine the number of grids in subdomain for all coordinate dimensions
        call para_range(1, gdm%nx, comm_1d_x%nprocs, comm_1d_x%myrank, sdm%ista, sdm%iend)
        sdm%nx = sdm%iend - sdm%ista + 1
        call para_range(1, gdm%ny, comm_1d_y%nprocs, comm_1d_y%myrank, sdm%jsta, sdm%jend)
        sdm%ny = sdm%jend - sdm%jsta + 1
        call para_range(1, gdm%nz, comm_1d_z%nprocs, comm_1d_z%myrank, sdm%ksta, sdm%kend)
        sdm%nz = sdm%kend - sdm%ksta + 1

        allocate( sdm%dxm(0:sdm%nx+1) )
        allocate( sdm%dym(0:sdm%ny+1) )
        allocate( sdm%dzm(0:sdm%nz+1) )

        allocate( sdm%dxg(0:sdm%nx+1) )
        allocate( sdm%dyg(0:sdm%ny+1) )
        allocate( sdm%dzg(0:sdm%nz+1) )

        allocate( sdm%xg(0:sdm%nx+1) )
        allocate( sdm%yg(0:sdm%ny+1) )
        allocate( sdm%zg(0:sdm%nz+1) )

        ! 2. Setup the mesh size, grid spacing and grid coordinates of subdomain for assigned ranges
        sdm%dxm(0:sdm%nx+1) = gdm%dxm(sdm%ista-1:sdm%iend+1)
        sdm%dym(0:sdm%ny+1) = gdm%dym(sdm%jsta-1:sdm%jend+1)
        sdm%dzm(0:sdm%nz+1) = gdm%dzm(sdm%ksta-1:sdm%kend+1)

        sdm%dxg(0:sdm%nx+1) = gdm%dxg(sdm%ista-1:sdm%iend+1)
        sdm%dyg(0:sdm%ny+1) = gdm%dyg(sdm%jsta-1:sdm%jend+1)
        sdm%dzg(0:sdm%nz+1) = gdm%dzg(sdm%ksta-1:sdm%kend+1)

        sdm%xg(0:sdm%nx+1) = gdm%xg(sdm%ista-1:sdm%iend+1)
        sdm%yg(0:sdm%ny+1) = gdm%yg(sdm%jsta-1:sdm%jend+1)
        sdm%zg(0:sdm%nz+1) = gdm%zg(sdm%ksta-1:sdm%kend+1)

        ! 3. Setup the starting coordinates of subdomain
        sdm%ox = sdm%xg(1) - 0.5d0 * sdm%dxm(1)
        sdm%oy = sdm%yg(1) - 0.5d0 * sdm%dym(1)
        sdm%oz = sdm%zg(1) - 0.5d0 * sdm%dzm(1)

        ! 4. Determine whether the subdomain is at boundary or not
        if(comm_1d_x%west_rank.eq.MPI_PROC_NULL) then
            sdm%is_x0_boundary = .true.
        else
            sdm%is_x0_boundary = .false.
        endif

        if(comm_1d_x%east_rank.eq.MPI_PROC_NULL) then
            sdm%is_x1_boundary = .true.
        else
            sdm%is_x1_boundary = .false.
        endif

        if(comm_1d_y%west_rank.eq.MPI_PROC_NULL) then
            sdm%is_y0_boundary = .true.
        else
            sdm%is_y0_boundary = .false.
        endif

        if(comm_1d_y%east_rank.eq.MPI_PROC_NULL) then
            sdm%is_y1_boundary = .true.
        else
            sdm%is_y1_boundary = .false.
        endif

        if(comm_1d_z%west_rank.eq.MPI_PROC_NULL) then
            sdm%is_z0_boundary = .true.
        else
            sdm%is_z0_boundary = .false.
        endif

        if(comm_1d_z%east_rank.eq.MPI_PROC_NULL) then
            sdm%is_z1_boundary = .true.
        else
            sdm%is_z1_boundary = .false.
        endif

        ! Allocate and initialize grid variables of subdomain
        allocate(sdm%x(0:sdm%nx+1, 0:sdm%ny+1, 0:sdm%nz+1))
        allocate(sdm%b(0:sdm%nx+1, 0:sdm%ny+1, 0:sdm%nz+1))
        allocate(sdm%r(0:sdm%nx+1, 0:sdm%ny+1, 0:sdm%nz+1))

        sdm%x(:,:,:) = 0.0d0
        sdm%b(:,:,:) = 0.0d0
        sdm%r(:,:,:) = 0.0d0

    end subroutine geometry_subdomain_create

    !>
    !> @brief       Free the variables of subdomain
    !> @param       sdm         Subdomain
    !>
    subroutine geometry_subdomain_destroy(sdm)

        implicit none
        type(subdomain), intent(inout)  :: sdm

        deallocate( sdm%dxm, sdm%dym, sdm%dzm )
        deallocate( sdm%dxg, sdm%dyg, sdm%dzg )
        deallocate( sdm%xg, sdm%yg, sdm%zg )
        deallocate( sdm%x, sdm%b, sdm%r)

    end subroutine geometry_subdomain_destroy

    !>
    !> @brief       Build derived datatypes for subdomain communication using ghostcells
    !> @param       sdm         Subdomain
    !>
    subroutine geometry_subdomain_ddt_create(sdm)

        use mpi

        implicit none

        type(subdomain), intent(inout)     :: sdm

        integer(kind=4) :: sizes(0:2), subsizes(0:2), starts(0:2), ierr     ! Local variables for MPI_Type_create_subarray

        ! Derived datatype of inner domain without ghostcell
        sizes    = (/sdm%nx+2,sdm%ny+2,sdm%nz+2/)
        subsizes = (/sdm%nx,  sdm%ny  ,sdm%nz  /)
        starts   = (/1,      1,      1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_inner_domain, ierr)
        call MPI_Type_commit(sdm%ddt_inner_domain,ierr)

        ! Derived datatype for sending data to east MPI process (x+ neighbor)
        sizes    = (/sdm%nx+2,sdm%ny+2,sdm%nz+2/)
        subsizes = (/      1,sdm%ny+2,sdm%nz+2/)
        starts   = (/sdm%nx,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_yz_plane_xn, ierr)
        call MPI_Type_commit(sdm%ddt_yz_plane_xn,ierr)

        ! Derived datatype for receiving data from west MPI process (x- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_yz_plane_x0, ierr)
        call MPI_Type_commit(sdm%ddt_yz_plane_x0,ierr)

        ! Derived datatype for sending data to west MPI process (x- neighbor)
        starts   = (/      1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_yz_plane_x1, ierr)
        call MPI_Type_commit(sdm%ddt_yz_plane_x1,ierr)

        ! Derived datatype for receiving data from east MPI process (x+ neighbor)
        starts   = (/  sdm%nx+1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_yz_plane_xn1, ierr)
        call MPI_Type_commit(sdm%ddt_yz_plane_xn1,ierr)

        ! Derived datatype for sending data to north MPI process (y+ neighbor)
        sizes    = (/sdm%nx+2,sdm%ny+2,sdm%nz+2/)
        subsizes = (/sdm%nx+2,      1,sdm%nz+2/)
        starts   = (/      0,sdm%ny,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xz_plane_yn, ierr)
        call MPI_Type_commit(sdm%ddt_xz_plane_yn,ierr)

        ! Derived datatype for receiving data from south MPI process (y- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xz_plane_y0, ierr)
        call MPI_Type_commit(sdm%ddt_xz_plane_y0,ierr)

        ! Derived datatype for sending data to south MPI process (y- neighbor)
        starts   = (/      0,      1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xz_plane_y1, ierr)
        call MPI_Type_commit(sdm%ddt_xz_plane_y1,ierr)

        ! Derived datatype for receiving data from north MPI process (y+ neighbor)
        starts   = (/      0,  sdm%ny+1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xz_plane_yn1, ierr)
        call MPI_Type_commit(sdm%ddt_xz_plane_yn1,ierr)

        ! Derived datatype for sending data to forth MPI process (z+ neighbor)
        sizes    = (/sdm%nx+2,sdm%ny+2,sdm%nz+2/)
        subsizes = (/sdm%nx+2,sdm%ny+2,      1/)
        starts   = (/      0,      0,sdm%nz/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xy_plane_zn, ierr)
        call MPI_Type_commit(sdm%ddt_xy_plane_zn,ierr)

        ! Derived datatype for receiving data from back MPI process (z- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xy_plane_z0, ierr)
        call MPI_Type_commit(sdm%ddt_xy_plane_z0,ierr)

        ! Derived datatype for sending data to back MPI process (z- neighbor)
        starts   = (/      0,      0,      1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xy_plane_z1, ierr)
        call MPI_Type_commit(sdm%ddt_xy_plane_z1,ierr)

        ! Derived datatype for receiving data from forth MPI process (z+ neighbor)
        starts   = (/      0,      0,  sdm%nz+1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xy_plane_zn1, ierr)
        call MPI_Type_commit(sdm%ddt_xy_plane_zn1,ierr)

    end subroutine geometry_subdomain_ddt_create

    !>
    !> @brief       Free derived datatypes for subdomain communication using ghostcells
    !> @param       sdm         subdomain
    !>
    subroutine geometry_subdomain_ddt_destroy(sdm)

        use mpi

        implicit none

        type(subdomain), intent(inout)     :: sdm
        integer(kind=4) :: ierr

        call MPI_Type_free(sdm%ddt_inner_domain, ierr)
        call MPI_Type_free(sdm%ddt_yz_plane_x0, ierr)
        call MPI_Type_free(sdm%ddt_yz_plane_x1, ierr)
        call MPI_Type_free(sdm%ddt_yz_plane_xn, ierr)
        call MPI_Type_free(sdm%ddt_yz_plane_xn1, ierr)

        call MPI_Type_free(sdm%ddt_xz_plane_y0, ierr)
        call MPI_Type_free(sdm%ddt_xz_plane_y1, ierr)
        call MPI_Type_free(sdm%ddt_xz_plane_yn, ierr)
        call MPI_Type_free(sdm%ddt_xz_plane_yn1, ierr)

        call MPI_Type_free(sdm%ddt_xy_plane_z0, ierr)
        call MPI_Type_free(sdm%ddt_xy_plane_z1, ierr)
        call MPI_Type_free(sdm%ddt_xy_plane_zn, ierr)
        call MPI_Type_free(sdm%ddt_xy_plane_zn1, ierr)

    end subroutine geometry_subdomain_ddt_destroy

    !>
    !> @brief       Update the values of boundary ghostcells through communication in all directions.
    !> @param       u           Grid variable to be updated
    !> @param       sdm         Subdomain
    !>
    subroutine geometry_halocell_update(u, sdm)

        use mpi
        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z

        implicit none
    
        real(kind=8), intent(inout)     :: u(0:,0:,0:)
        type(subdomain), intent(in)     :: sdm
    
        integer(kind=4)     :: ierr
        integer(kind=4)     :: request(4)
    
        ! Update the ghostcells in x-direction using derived datatypes and subcommunicator
        call MPI_Isend(u, 1, sdm%ddt_yz_plane_xn,  comm_1d_x%east_rank, 111, comm_1d_x%mpi_comm, request(1), ierr)
        call MPI_Irecv(u, 1, sdm%ddt_yz_plane_x0,  comm_1d_x%west_rank, 111, comm_1d_x%mpi_comm, request(2), ierr)
        call MPI_Isend(u, 1, sdm%ddt_yz_plane_x1,  comm_1d_x%west_rank, 222, comm_1d_x%mpi_comm, request(3), ierr)
        call MPI_Irecv(u, 1, sdm%ddt_yz_plane_xn1, comm_1d_x%east_rank, 222, comm_1d_x%mpi_comm, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
        
        ! Update the ghostcells in y-direction using derived datatypes and subcommunicator
        call MPI_Isend(u, 1, sdm%ddt_xz_plane_yn,  comm_1d_y%east_rank, 333, comm_1d_y%mpi_comm, request(1), ierr)
        call MPI_Irecv(u, 1, sdm%ddt_xz_plane_y0,  comm_1d_y%west_rank, 333, comm_1d_y%mpi_comm, request(2), ierr)
        call MPI_Isend(u, 1, sdm%ddt_xz_plane_y1,  comm_1d_y%west_rank, 444, comm_1d_y%mpi_comm, request(3), ierr)
        call MPI_Irecv(u, 1, sdm%ddt_xz_plane_yn1, comm_1d_y%east_rank, 444, comm_1d_y%mpi_comm, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
        
        ! Update the ghostcells in z-direction using derived datatypes and subcommunicator
        call MPI_Isend(u, 1, sdm%ddt_xy_plane_zn,  comm_1d_z%east_rank, 555, comm_1d_z%mpi_comm, request(1) , ierr)
        call MPI_Irecv(u, 1, sdm%ddt_xy_plane_z0,  comm_1d_z%west_rank, 555, comm_1d_z%mpi_comm, request(2), ierr)
        call MPI_Isend(u, 1, sdm%ddt_xy_plane_z1,  comm_1d_z%west_rank, 666, comm_1d_z%mpi_comm, request(3), ierr)
        call MPI_Irecv(u, 1, sdm%ddt_xy_plane_zn1, comm_1d_z%east_rank, 666, comm_1d_z%mpi_comm, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
    
    end subroutine geometry_halocell_update

    !>
    !> @brief       Update the values of boundary ghostcells through communication when the domain is not aggregated.
    !> @param       u           Grid variable to be updated
    !> @param       sdm         Subdomain
    !> @param       is_serial   Boolean whether the subdomain is aggregated(.true.) or not(.false.)
    !>
    subroutine geometry_halocell_update_selectively(u, sdm, is_serial)

        use mpi
        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z

        implicit none
    
        real(kind=8), intent(inout)     :: u(0:,0:,0:)
        type(subdomain), intent(in)     :: sdm
        logical, intent(in)             :: is_serial(0:2)

        integer(kind=4)     :: ierr
        integer(kind=4)     :: request(4)
    
        ! Update the ghostcells in x-direction using derived datatypes and subcommunicator when the domain is not aggregated
        if(is_serial(0) .eq. .false.) then
            call MPI_Isend(u, 1, sdm%ddt_yz_plane_xn,  comm_1d_x%east_rank, 111, comm_1d_x%mpi_comm, request(1), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_yz_plane_x0,  comm_1d_x%west_rank, 111, comm_1d_x%mpi_comm, request(2), ierr)
            call MPI_Isend(u, 1, sdm%ddt_yz_plane_x1,  comm_1d_x%west_rank, 222, comm_1d_x%mpi_comm, request(3), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_yz_plane_xn1, comm_1d_x%east_rank, 222, comm_1d_x%mpi_comm, request(4), ierr)
            call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)    
        endif
        if(is_serial(1) .eq. .false.) then
        ! Update the ghostcells in y-direction using derived datatypes and subcommunicator when the domain is not aggregated
            call MPI_Isend(u, 1, sdm%ddt_xz_plane_yn,  comm_1d_y%east_rank, 333, comm_1d_y%mpi_comm, request(1), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_xz_plane_y0,  comm_1d_y%west_rank, 333, comm_1d_y%mpi_comm, request(2), ierr)
            call MPI_Isend(u, 1, sdm%ddt_xz_plane_y1,  comm_1d_y%west_rank, 444, comm_1d_y%mpi_comm, request(3), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_xz_plane_yn1, comm_1d_y%east_rank, 444, comm_1d_y%mpi_comm, request(4), ierr)
            call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)    
        endif
        if(is_serial(2) .eq. .false.) then
        ! Update the ghostcells in z-direction using derived datatypes and subcommunicator when the domain is not aggregated
            call MPI_Isend(u, 1, sdm%ddt_xy_plane_zn,  comm_1d_z%east_rank, 555, comm_1d_z%mpi_comm, request(1) , ierr)
            call MPI_Irecv(u, 1, sdm%ddt_xy_plane_z0,  comm_1d_z%west_rank, 555, comm_1d_z%mpi_comm, request(2), ierr)
            call MPI_Isend(u, 1, sdm%ddt_xy_plane_z1,  comm_1d_z%west_rank, 666, comm_1d_z%mpi_comm, request(3), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_xy_plane_zn1, comm_1d_z%east_rank, 666, comm_1d_z%mpi_comm, request(4), ierr)
            call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)    
        endif

    end subroutine geometry_halocell_update_selectively

end module geometry
