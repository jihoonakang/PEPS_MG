module geometry

    implicit none

    type, public    :: domain
        integer(kind=4)     :: nx, ny, nz
        real(kind=8)        :: lx, ly, lz
        real(kind=8)        :: ox, oy, oz
        real(kind=8), allocatable, dimension(:)         :: dxm, dym, dzm
        real(kind=8), allocatable, dimension(:)         :: dxg, dyg, dzg
        real(kind=8), allocatable, dimension(:)         :: xg, yg, zg
        logical                                         :: is_periodic(0:2)
    end type domain

    type, public    :: subdomain
        integer(kind=4)     :: nx, ny, nz
        real(kind=8)        :: lx, ly, lz
        real(kind=8)        :: ox, oy, oz
        real(kind=8), allocatable, dimension(:)         :: dxm, dym, dzm
        real(kind=8), allocatable, dimension(:)         :: dxg, dyg, dzg
        real(kind=8), allocatable, dimension(:)         :: xg, yg, zg
        real(kind=8), allocatable, dimension(:,:,:)     :: x, b, r
        logical                                         :: is_periodic(0:2)

        logical                 :: is_aggregated(0:2)
        integer(kind=4), public :: ista, iend, jsta, jend, ksta, kend
        integer(kind=4), public :: ddt_yz_plane_x0, ddt_yz_plane_x1, ddt_yz_plane_xn, ddt_yz_plane_xn1
        integer(kind=4), public :: ddt_xz_plane_y0, ddt_xz_plane_y1, ddt_xz_plane_yn, ddt_xz_plane_yn1
        integer(kind=4), public :: ddt_xy_plane_z0, ddt_xy_plane_z1, ddt_xy_plane_zn, ddt_xy_plane_zn1
        integer(kind=4), public :: ddt_inner_domain
        logical(kind=1)         :: is_x0_boundary, is_x1_boundary, is_y0_boundary, is_y1_boundary, is_z0_boundary, is_z1_boundary
    end type subdomain

    contains

    subroutine geometry_domain_create(gdm, nx, ny, nz, ox, oy, oz, lx, ly, lz, ax, ay, az, period)
        
        implicit none

        type(domain), intent(inout)     :: gdm
        integer(kind=4), intent(in)     :: nx, ny, nz
        real(kind=8),    intent(in)     :: ox, oy, oz
        real(kind=8),    intent(in)     :: lx, ly, lz
        real(kind=8),    intent(in)     :: ax, ay, az
        logical,         intent(in)     :: period(0:2)

        integer(kind=4) :: i

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
    
        if(ax.eq.1.0d0) then
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

        if(ay.eq.1.0d0) then
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

        if(az.eq.1.0d0) then
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

    subroutine geometry_domain_destroy(gdm)

        implicit none

        type(domain), intent(inout)     :: gdm

        deallocate(gdm%dxm, gdm%dym, gdm%dzm)
        deallocate(gdm%dxg, gdm%dyg, gdm%dzg)
        deallocate(gdm%xg, gdm%yg, gdm%zg)
    
    end subroutine geometry_domain_destroy

    subroutine geometry_subdomain_create(sdm, gdm)

        use mpi
        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z
        implicit none

        type(subdomain), intent(inout)  :: sdm
        type(domain), intent(in)        :: gdm

        sdm%is_periodic(0:2)    = gdm%is_periodic(0:2)


        if(comm_1d_x%nprocs.eq.1) then
            sdm%is_aggregated(0)  = .true.
        else
            sdm%is_aggregated(0)  = .false.
        endif

        if(comm_1d_y%nprocs.eq.1) then
            sdm%is_aggregated(1) = .true.
        else
            sdm%is_aggregated(1) = .false.
        endif

        if(comm_1d_z%nprocs.eq.1) then
            sdm%is_aggregated(2) = .true.
        else
            sdm%is_aggregated(2) = .false.
        endif
        
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

        sdm%dxm(0:sdm%nx+1) = gdm%dxm(sdm%ista-1:sdm%iend+1)
        sdm%dym(0:sdm%ny+1) = gdm%dym(sdm%jsta-1:sdm%jend+1)
        sdm%dzm(0:sdm%nz+1) = gdm%dzm(sdm%ksta-1:sdm%kend+1)

        sdm%dxg(0:sdm%nx+1) = gdm%dxg(sdm%ista-1:sdm%iend+1)
        sdm%dyg(0:sdm%ny+1) = gdm%dyg(sdm%jsta-1:sdm%jend+1)
        sdm%dzg(0:sdm%nz+1) = gdm%dzg(sdm%ksta-1:sdm%kend+1)

        sdm%xg(0:sdm%nx+1) = gdm%xg(sdm%ista-1:sdm%iend+1)
        sdm%yg(0:sdm%ny+1) = gdm%yg(sdm%jsta-1:sdm%jend+1)
        sdm%zg(0:sdm%nz+1) = gdm%zg(sdm%ksta-1:sdm%kend+1)

        sdm%ox = sdm%xg(1) - 0.5d0 * sdm%dxm(1)
        sdm%oy = sdm%yg(1) - 0.5d0 * sdm%dym(1)
        sdm%oz = sdm%zg(1) - 0.5d0 * sdm%dzm(1)

        allocate(sdm%x(0:sdm%nx+1, 0:sdm%ny+1, 0:sdm%nz+1))
        allocate(sdm%b(0:sdm%nx+1, 0:sdm%ny+1, 0:sdm%nz+1))
        allocate(sdm%r(0:sdm%nx+1, 0:sdm%ny+1, 0:sdm%nz+1))

        sdm%x(:,:,:) = 0.0d0
        sdm%b(:,:,:) = 0.0d0
        sdm%r(:,:,:) = 0.0d0

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

        ! print *, myrank, 'xbc', comm_1d_x%west_rank, comm_1d_x%east_rank, sdm%is_x0_boundary, sdm%is_x1_boundary
        ! print *, myrank, 'ybc', comm_1d_y%west_rank, comm_1d_y%east_rank, sdm%is_y0_boundary, sdm%is_y1_boundary
        ! print *, myrank, 'zbc', comm_1d_z%west_rank, comm_1d_z%east_rank, sdm%is_z0_boundary, sdm%is_z1_boundary

    end subroutine geometry_subdomain_create

    subroutine geometry_subdomain_destroy(sdm)

        implicit none
        type(subdomain), intent(inout)  :: sdm

        deallocate( sdm%dxm, sdm%dym, sdm%dzm )
        deallocate( sdm%dxg, sdm%dyg, sdm%dzg )
        deallocate( sdm%xg, sdm%yg, sdm%zg )
        deallocate( sdm%x, sdm%b, sdm%r)

    end subroutine geometry_subdomain_destroy

    subroutine geometry_subdomain_ddt_create(sdm)

        use mpi

        implicit none

        type(subdomain), intent(inout)     :: sdm

        integer(kind=4) :: sizes(0:2), subsizes(0:2), starts(0:2), ierr     ! Local variables for MPI_Type_create_subarray

        ! ddtype sending data to east MPI process (x+ neighbor)
        sizes    = (/sdm%nx+2,sdm%ny+2,sdm%nz+2/)
        subsizes = (/sdm%nx,  sdm%ny  ,sdm%nz  /)
        starts   = (/1,      1,      1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_inner_domain, ierr)
        call MPI_Type_commit(sdm%ddt_inner_domain,ierr)

        ! ddtype sending data to east MPI process (x+ neighbor)
        sizes    = (/sdm%nx+2,sdm%ny+2,sdm%nz+2/)
        subsizes = (/      1,sdm%ny+2,sdm%nz+2/)
        starts   = (/sdm%nx,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_yz_plane_xn, ierr)
        call MPI_Type_commit(sdm%ddt_yz_plane_xn,ierr)

        ! ddtype receiving data from west MPI process (x- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_yz_plane_x0, ierr)
        call MPI_Type_commit(sdm%ddt_yz_plane_x0,ierr)

        ! ddtype sending data to west MPI process (x- neighbor)
        starts   = (/      1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_yz_plane_x1, ierr)
        call MPI_Type_commit(sdm%ddt_yz_plane_x1,ierr)

        ! ddtype receiving data from east MPI process (x+ neighbor)
        starts   = (/  sdm%nx+1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_yz_plane_xn1, ierr)
        call MPI_Type_commit(sdm%ddt_yz_plane_xn1,ierr)

        ! ddtype sending data to north MPI process (y+ neighbor)
        sizes    = (/sdm%nx+2,sdm%ny+2,sdm%nz+2/)
        subsizes = (/sdm%nx+2,      1,sdm%nz+2/)
        starts   = (/      0,sdm%ny,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xz_plane_yn, ierr)
        call MPI_Type_commit(sdm%ddt_xz_plane_yn,ierr)

        ! ddtype receiving data from south MPI process (y- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xz_plane_y0, ierr)
        call MPI_Type_commit(sdm%ddt_xz_plane_y0,ierr)

        ! ddtype sending data to south MPI process (y- neighbor)
        starts   = (/      0,      1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xz_plane_y1, ierr)
        call MPI_Type_commit(sdm%ddt_xz_plane_y1,ierr)

        ! ddtype receiving data from north MPI process (y+ neighbor)
        starts   = (/      0,  sdm%ny+1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xz_plane_yn1, ierr)
        call MPI_Type_commit(sdm%ddt_xz_plane_yn1,ierr)

        ! ddtype sending data to forth MPI process (z+ neighbor)
        sizes    = (/sdm%nx+2,sdm%ny+2,sdm%nz+2/)
        subsizes = (/sdm%nx+2,sdm%ny+2,      1/)
        starts   = (/      0,      0,sdm%nz/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xy_plane_zn, ierr)
        call MPI_Type_commit(sdm%ddt_xy_plane_zn,ierr)

        ! ddtype receiving data from back MPI process (z- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xy_plane_z0, ierr)
        call MPI_Type_commit(sdm%ddt_xy_plane_z0,ierr)

        ! ddtype sending data to back MPI process (z- neighbor)
        starts   = (/      0,      0,      1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xy_plane_z1, ierr)
        call MPI_Type_commit(sdm%ddt_xy_plane_z1,ierr)

        ! ddtype receiving data from forth MPI process (z+ neighbor)
        starts   = (/      0,      0,  sdm%nz+1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, sdm%ddt_xy_plane_zn1, ierr)
        call MPI_Type_commit(sdm%ddt_xy_plane_zn1,ierr)

    end subroutine geometry_subdomain_ddt_create

    subroutine geometry_subdomain_ddt_destroy(sdm)

        use mpi

        implicit none

        type(subdomain), intent(inout)     :: sdm
        integer(kind=4) :: ierr

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

    subroutine geometry_halocell_update_selectively(u, sdm, is_serial)

        use mpi
        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z

        implicit none
    
        real(kind=8), intent(inout)     :: u(0:,0:,0:)
        type(subdomain), intent(in)     :: sdm
        logical, intent(in)             :: is_serial(0:2)

        integer(kind=4)     :: ierr
        integer(kind=4)     :: request(4)
    
        if(is_serial(0) .eq. .false.) then
            ! Update the ghostcells in x-direction using derived datatypes and subcommunicator
            call MPI_Isend(u, 1, sdm%ddt_yz_plane_xn,  comm_1d_x%east_rank, 111, comm_1d_x%mpi_comm, request(1), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_yz_plane_x0,  comm_1d_x%west_rank, 111, comm_1d_x%mpi_comm, request(2), ierr)
            call MPI_Isend(u, 1, sdm%ddt_yz_plane_x1,  comm_1d_x%west_rank, 222, comm_1d_x%mpi_comm, request(3), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_yz_plane_xn1, comm_1d_x%east_rank, 222, comm_1d_x%mpi_comm, request(4), ierr)
            call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)    
        endif
        if(is_serial(1) .eq. .false.) then
            ! Update the ghostcells in y-direction using derived datatypes and subcommunicator
            call MPI_Isend(u, 1, sdm%ddt_xz_plane_yn,  comm_1d_y%east_rank, 333, comm_1d_y%mpi_comm, request(1), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_xz_plane_y0,  comm_1d_y%west_rank, 333, comm_1d_y%mpi_comm, request(2), ierr)
            call MPI_Isend(u, 1, sdm%ddt_xz_plane_y1,  comm_1d_y%west_rank, 444, comm_1d_y%mpi_comm, request(3), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_xz_plane_yn1, comm_1d_y%east_rank, 444, comm_1d_y%mpi_comm, request(4), ierr)
            call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)    
        endif
        if(is_serial(2) .eq. .false.) then
            ! Update the ghostcells in z-direction using derived datatypes and subcommunicator
            call MPI_Isend(u, 1, sdm%ddt_xy_plane_zn,  comm_1d_z%east_rank, 555, comm_1d_z%mpi_comm, request(1) , ierr)
            call MPI_Irecv(u, 1, sdm%ddt_xy_plane_z0,  comm_1d_z%west_rank, 555, comm_1d_z%mpi_comm, request(2), ierr)
            call MPI_Isend(u, 1, sdm%ddt_xy_plane_z1,  comm_1d_z%west_rank, 666, comm_1d_z%mpi_comm, request(3), ierr)
            call MPI_Irecv(u, 1, sdm%ddt_xy_plane_zn1, comm_1d_z%east_rank, 666, comm_1d_z%mpi_comm, request(4), ierr)
            call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)    
        endif

    end subroutine geometry_halocell_update_selectively

    subroutine geometry_boundary_value_calculate(u, sdm, gdm, comm_boundary)

        use mpi
        use mpi_topology, only : myrank, boundary_comm

        implicit none

        real(kind=8), intent(inout)     :: u(0:,0:,0:)
        type(subdomain), intent(in)     :: sdm
        type(domain), intent(in)        :: gdm
        type(boundary_comm), intent(in) :: comm_boundary

        real(kind=8), parameter         :: PI = 3.1415926535897932384626d0
        real(kind=8), allocatable, dimension(:,:,:)     :: bval_x, bval_y, bval_z
#ifndef INPLACE
        real(kind=8), allocatable, dimension(:,:,:)     :: bval_xg, bval_yg, bval_zg
#endif

        integer(kind=4) :: i, j, k, ip, jp, kp, irp, jrp, krp
        integer(kind=4) :: nxs, nys, nzs, nxg, nyg, nzg
        integer(kind=4) :: ierr
        real(kind=8)    :: dr, drx, dry, drz, gf_coeff
        real(kind=8)    :: dudx0, dudx1, dudy0, dudy1, dudz0, dudz1
        real(kind=8)    :: alpha, dxg1, dxg2, dyg1, dyg2, dzg1, dzg2, w1, w2

        allocate(bval_x(gdm%ny,gdm%nz,0:1))
        allocate(bval_y(gdm%nx,gdm%nz,0:1))
        allocate(bval_z(gdm%nx,gdm%ny,0:1))

        bval_x = 0.0d0
        bval_y = 0.0d0
        bval_z = 0.0d0

#ifndef INPLACE
        allocate(bval_xg(gdm%ny,gdm%nz,0:1))
        allocate(bval_yg(gdm%nx,gdm%nz,0:1))
        allocate(bval_zg(gdm%nx,gdm%ny,0:1))

        bval_xg = 0.0d0
        bval_yg = 0.0d0
        bval_zg = 0.0d0
#endif

        nxg = gdm%nx
        nyg = gdm%ny
        nzg = gdm%nz

        nxs = sdm%nx
        nys = sdm%ny
        nzs = sdm%nz

        dudx0 = 0.0d0
        dudx1 = 0.0d0
        dudy0 = 0.0d0
        dudy1 = 0.0d0
        dudz0 = 0.0d0
        dudz1 = 0.0d0

        gf_coeff = 0.25d0 / PI

        ! Contributions from boundaries in x-direction
        if(sdm%is_x0_boundary .eq. .true.) then
            dxg1 = sdm%dxg(2)
            dxg2 = sdm%dxg(3)
            alpha = dxg2 / dxg1
            w1 = gf_coeff / dxg1 * (2.0d0+alpha) / (1.0d0+alpha)
            w2 = gf_coeff / dxg2 / (1.0d0+alpha)
            do kp = 1, nzs
                do jp = 1, nys
                    jrp = jp+sdm%jsta-1
                    krp = kp+sdm%ksta-1
                    dudx0 = -((u(2,jp,kp)-u(1,jp,kp))*w1 - (u(3,jp,kp)-u(2,jp,kp))*w2) * sdm%dym(jp) * sdm%dzm(kp)
                    
                    ! Boundary treatment in x-direction
                    do k = 1, nzg
                        do j = 1, nyg
                            ! Lower boundary in x-direction
                            drx = gdm%xg(0) - gdm%xg(1)
                            dry = gdm%yg(j) - gdm%yg(jrp)
                            drz = gdm%zg(k) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,0) = bval_x(j,k,0) + dudx0 / dr

                            ! Upper boundary in x-direction
                            drx = gdm%xg(nxg+1) - gdm%xg(1)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,1) = bval_x(j,k,1) + dudx0 / dr
                        enddo
                    enddo

                    ! Boundary treatment in y-direction
                    do k = 1, nzg
                        do i = 1, nxg
                            ! Lower boundary in y-direction
                            drx = gdm%xg(i) - gdm%xg(1)
                            dry = gdm%yg(0) - gdm%yg(jrp)
                            drz = gdm%zg(k) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,0) = bval_y(i,k,0) + dudx0 / dr

                            ! Upper boundary in y-direction
                            dry = gdm%yg(nyg+1) - gdm%yg(jrp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,1) = bval_y(i,k,1) + dudx0 / dr
                        enddo
                    enddo

                    ! Boundary treatment in z-direction
                    do j = 1, nyg
                        do i = 1, nxg
                            ! Lower boundary in z-direction
                            drx = gdm%xg(i) - gdm%xg(1)
                            dry = gdm%yg(j) - gdm%yg(jrp)
                            drz = gdm%zg(0) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,0) = bval_z(i,j,0) + dudx0 / dr

                            ! Upper boundary in z-direction
                            drz = gdm%zg(nzg+1) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,1) = bval_z(i,j,1) + dudx0 / dr
                        enddo
                    enddo
                enddo
            enddo
        endif

        ! Contributions from boundaries in x-direction
        if(sdm%is_x1_boundary .eq. .true.) then
            dxg1 = sdm%dxg(nxs-1)
            dxg2 = sdm%dxg(nxs)
            alpha = dxg1 / dxg2
            w1 = gf_coeff / dxg2 * (2.0d0+alpha) / (1.0d0+alpha)
            w2 = gf_coeff / dxg1 / (1.0d0+alpha)
            do kp = 1, nzs
                do jp = 1, nys
                    jrp = jp+sdm%jsta-1
                    krp = kp+sdm%ksta-1
                    dudx1 = ((u(nxs,jp,kp)-u(nxs-1,jp,kp))*w1 - (u(nxs-1,jp,kp)-u(nxs-2,jp,kp))*w2) * sdm%dym(jp) * sdm%dzm(kp)
                
                    ! Boundary treatment in x-direction
                    do k = 1, nzg
                        do j = 1, nyg
                            ! Lower boundary in x-direction
                            drx = gdm%xg(0) - gdm%xg(nxg)
                            dry = gdm%yg(j) - gdm%yg(jrp)
                            drz = gdm%zg(k) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,0) = bval_x(j,k,0) + dudx1 / dr

                            ! Upper boundary in x-direction
                            drx = gdm%xg(nxg+1) - gdm%xg(nxg)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,1) = bval_x(j,k,1) + dudx1 / dr
                        enddo
                    enddo

                    ! Boundary treatment in y-direction
                    do k = 1, nzg
                        do i = 1, nxg
                            ! Lower boundary in y-direction
                            drx = gdm%xg(i) - gdm%xg(nxg)
                            dry = gdm%yg(0) - gdm%yg(jrp)
                            drz = gdm%zg(k) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,0) = bval_y(i,k,0) + dudx1 / dr

                            ! Upper boundary in y-direction
                            dry = gdm%yg(nyg+1) - gdm%yg(jrp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,1) = bval_y(i,k,1) + dudx1 / dr
                        enddo
                    enddo

                    ! Boundary treatment in z-direction
                    do j = 1, nyg
                        do i = 1, nxg
                            ! Lower boundary in z-direction
                            drx = gdm%xg(i) - gdm%xg(nxg)
                            dry = gdm%yg(j) - gdm%yg(jrp)
                            drz = gdm%zg(0) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,0) = bval_z(i,j,0) + dudx1 / dr

                            ! Upper boundary in z-direction
                            drz = gdm%zg(nzg+1) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,1) = bval_z(i,j,1) + dudx1 / dr
                        enddo
                    enddo
                enddo
            enddo
        endif

        ! Contributions from boundaries in y-direction
        if(sdm%is_y0_boundary .eq. .true.) then
            dyg1 = sdm%dyg(2)
            dyg2 = sdm%dyg(3)
            alpha = dyg2 / dyg1
            w1 = gf_coeff / dyg1 * (2.0d0+alpha) / (1.0d0+alpha)
            w2 = gf_coeff / dyg2 / (1.0d0+alpha)
            do kp = 1, nzs
                do ip = 1, nxs
                    dudy0 = -((u(ip,2,kp)-u(ip,1,kp)) * w1 - (u(ip,3,kp)-u(ip,2,kp)) * w2) * sdm%dxm(ip) * sdm%dzm(kp)
                    irp = ip+sdm%ista-1
                    krp = kp+sdm%ksta-1

                    ! Boundary treatment in x-direction
                    do k = 1, nzg
                        do j = 1, nyg
                            ! Lower boundary in x-direction
                            drx = gdm%xg(0) - gdm%xg(irp)
                            dry = gdm%yg(j) - gdm%yg(1)
                            drz = gdm%zg(k) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,0) = bval_x(j,k,0) + dudy0 / dr

                            ! Upper boundary in x-direction
                            drx = gdm%xg(nxg+1) - gdm%xg(irp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,1) = bval_x(j,k,1) + dudy0 / dr
                        enddo
                    enddo

                    ! Boundary treatment in y-direction
                    do k = 1, nzg
                        do i = 1, nxg
                            ! Lower boundary in y-direction
                            drx = gdm%xg(i) - gdm%xg(irp)
                            dry = gdm%yg(0) - gdm%yg(1)
                            drz = gdm%zg(k) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,0) = bval_y(i,k,0) + dudy0 / dr

                            ! Upper boundary in y-direction
                            dry = gdm%yg(nyg+1) - gdm%yg(1)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,1) = bval_y(i,k,1) + dudy0 / dr
                        enddo
                    enddo

                    ! Boundary treatment in z-direction
                    do j = 1, nyg
                        do i = 1, nxg
                            ! Lower boundary in z-direction
                            drx = gdm%xg(i) - gdm%xg(irp)
                            dry = gdm%yg(j) - gdm%yg(1)
                            drz = gdm%zg(0) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,0) = bval_z(i,j,0) + dudy0 / dr

                            ! Upper boundary in z-direction
                            drz = gdm%zg(nzg+1) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,1) = bval_z(i,j,1) + dudy0 / dr
                        enddo
                    enddo
                enddo
            enddo
        endif

        ! Contributions from boundaries in y-direction
        if(sdm%is_y1_boundary .eq. .true.) then
            dyg1 = sdm%dyg(nys-1)
            dyg2 = sdm%dyg(nys)
            alpha = dyg1 / dyg2
            w1 = gf_coeff / dyg2 * (2.0d0+alpha) / (1.0d0+alpha)
            w2 = gf_coeff / dyg1 / (1.0d0+alpha)
            do kp = 1, nzs
                do ip = 1, nxs
                    irp = ip+sdm%ista-1
                    krp = kp+sdm%ksta-1
                    dudy1 = ((u(ip,nys,kp)-u(ip,nys-1,kp)) * w1 - (u(ip,nys-1,kp)-u(ip,nys-2,kp)) * w2) * sdm%dxm(ip) * sdm%dzm(kp)

                    ! Boundary treatment in x-direction
                    do k = 1, nzg
                        do j = 1, nyg
                            ! Lower boundary in x-direction
                            drx = gdm%xg(0) - gdm%xg(irp)
                            dry = gdm%yg(j) - gdm%yg(nyg)
                            drz = gdm%zg(k) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,0) = bval_x(j,k,0) + dudy1 / dr

                            ! Upper boundary in x-direction
                            drx = gdm%xg(nxg+1) - gdm%xg(irp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,1) = bval_x(j,k,1) + dudy1 / dr
                        enddo
                    enddo

                    ! Boundary treatment in y-direction
                    do k = 1, nzg
                        do i = 1, nxg
                            ! Lower boundary in y-direction
                            drx = gdm%xg(i) - gdm%xg(irp)
                            dry = gdm%yg(0) - gdm%yg(nyg)
                            drz = gdm%zg(k) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,0) = bval_y(i,k,0) + dudy1 / dr

                            ! Upper boundary in y-direction
                            dry = gdm%yg(nyg+1) - gdm%yg(nyg)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,1) = bval_y(i,k,1) + dudy1 / dr
                        enddo
                    enddo

                    ! Boundary treatment in z-direction
                    do j = 1, nyg
                        do i = 1, nxg
                            ! Lower boundary in z-direction
                            drx = gdm%xg(i) - gdm%xg(irp)
                            dry = gdm%yg(j) - gdm%yg(nyg)
                            drz = gdm%zg(0) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,0) = bval_z(i,j,0) + dudy1 / dr

                            ! Upper boundary in z-direction
                            drz = gdm%zg(nzg+1) - gdm%zg(krp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,1) = bval_z(i,j,1) + dudy1 / dr
                        enddo
                    enddo
                enddo
            enddo
        endif

        ! Contributions from boundaries in z-direction
        if(sdm%is_z0_boundary .eq. .true.) then
            dzg1 = sdm%dzg(2)
            dzg2 = sdm%dzg(3)
            alpha = dzg2 / dzg1
            w1 = gf_coeff / dzg1 * (2.0d0+alpha) / (1.0d0+alpha)
            w2 = gf_coeff / dzg2 / (1.0d0+alpha)
            do jp = 1, nys
                do ip = 1, nxs
                    irp = ip+sdm%ista-1
                    jrp = jp+sdm%jsta-1
                    dudz0 = -((u(ip,jp,2)-u(ip,jp,1)) * w1 - (u(ip,jp,3)-u(ip,jp,2)) * w2) * sdm%dxm(ip) * sdm%dym(jp)

                    ! Boundary treatment in x-direction
                    do k = 1, nzg
                        do j = 1, nyg
                            ! Lower boundary in x-direction
                            drx = gdm%xg(0) - gdm%xg(irp)
                            dry = gdm%yg(j) - gdm%yg(jrp)
                            drz = gdm%zg(k) - gdm%zg(1)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,0) = bval_x(j,k,0) + dudz0 / dr

                            ! Upper boundary in x-direction
                            drx = gdm%xg(nxg+1) - gdm%xg(irp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,1) = bval_x(j,k,1) + dudz0 / dr
                        enddo
                    enddo

                    ! Boundary treatment in y-direction
                    do k = 1, nzg
                        do i = 1, nxg
                            ! Lower boundary in y-direction
                            drx = gdm%xg(i) - gdm%xg(irp)
                            dry = gdm%yg(0) - gdm%yg(jrp)
                            drz = gdm%zg(k) - gdm%zg(1)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,0) = bval_y(i,k,0) + dudz0 / dr

                            ! Upper boundary in y-direction
                            dry = gdm%yg(nyg+1) - gdm%yg(jrp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,1) = bval_y(i,k,1) + dudz0 / dr
                        enddo
                    enddo

                    ! Boundary treatment in z-direction
                    do j = 1, nyg
                        do i = 1, nxg
                            ! Lower boundary in z-direction
                            drx = gdm%xg(i) - gdm%xg(irp)
                            dry = gdm%yg(j) - gdm%yg(jrp)
                            drz = gdm%zg(0) - gdm%zg(1)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,0) = bval_z(i,j,0) + dudz0 / dr

                            ! Upper boundary in z-direction
                            drz = gdm%zg(nzg+1) - gdm%zg(1)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,1) = bval_z(i,j,1) + dudz0 / dr
                        enddo
                    enddo
                enddo
            enddo
        endif

        ! Contributions from boundaries in z-direction
        if(sdm%is_z1_boundary .eq. .true.) then
            dzg1 = sdm%dzg(nzs-1)
            dzg2 = sdm%dzg(nzs)
            alpha = dzg1 / dzg2
            w1 = gf_coeff / dzg2 * (2.0d0+alpha) / (1.0d0+alpha)
            w2 = gf_coeff / dzg1 / (1.0d0+alpha)
            do jp = 1, nys
                do ip = 1, nxs
                    irp = ip+sdm%ista-1
                    jrp = jp+sdm%jsta-1
                    dudz1 = ((u(ip,jp,nzs)-u(ip,jp,nzs-1)) * w1 - (u(ip,jp,nzs-1)-u(ip,jp,nzs-2)) * w2) * sdm%dxm(ip) * sdm%dym(jp)

                    ! Boundary treatment in x-direction
                    do k = 1, nzg
                        do j = 1, nyg
                            ! Lower boundary in x-direction
                            drx = gdm%xg(0) - gdm%xg(irp)
                            dry = gdm%yg(j) - gdm%yg(jrp)
                            drz = gdm%zg(k) - gdm%zg(nzg)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,0) = bval_x(j,k,0) + dudz1 / dr

                            ! Upper boundary in x-direction
                            drx = gdm%xg(nxg+1) - gdm%xg(irp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_x(j,k,1) = bval_x(j,k,1) + dudz1 / dr
                        enddo
                    enddo

                    ! Boundary treatment in y-direction
                    do k = 1, nzg
                        do i = 1, nxg
                            ! Lower boundary in y-direction
                            drx = gdm%xg(i) - gdm%xg(irp)
                            dry = gdm%yg(0) - gdm%yg(jrp)
                            drz = gdm%zg(k) - gdm%zg(nzg)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,0) = bval_y(i,k,0) + dudz1 / dr

                            ! Upper boundary in y-direction
                            dry = gdm%yg(nyg+1) - gdm%yg(jrp)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_y(i,k,1) = bval_y(i,k,1) + dudz1 / dr
                        enddo
                    enddo

                    ! Boundary treatment in z-direction
                    do j = 1, nyg
                        do i = 1, nxg
                            ! Lower boundary in z-direction
                            drx = gdm%xg(i) - gdm%xg(irp)
                            dry = gdm%yg(j) - gdm%yg(jrp)
                            drz = gdm%zg(0) - gdm%zg(nzg)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,0) = bval_z(i,j,0) + dudz1 / dr

                            ! Upper boundary in z-direction
                            drz = gdm%zg(nzg+1) - gdm%zg(nzg)
                            dr = sqrt(drx*drx+dry*dry+drz*drz)
                            bval_z(i,j,1) = bval_z(i,j,1) + dudz1 / dr
                        enddo
                    enddo
                enddo
            enddo
        endif
        
#ifdef MPI_INPLACE
        if(comm_boundary%mpi_comm.ne.MPI_COMM_NULL) then
            call MPI_Allreduce(MPI_IN_PLACE, bval_x, size(bval_x), MPI_REAL8, MPI_SUM, comm_boundary%mpi_comm, ierr)
            call MPI_Allreduce(MPI_IN_PLACE, bval_y, size(bval_y), MPI_REAL8, MPI_SUM, comm_boundary%mpi_comm, ierr)
            call MPI_Allreduce(MPI_IN_PLACE, bval_z, size(bval_z), MPI_REAL8, MPI_SUM, comm_boundary%mpi_comm, ierr)
        endif

        if(sdm%is_x0_boundary .eq. .true.) then
            do k = 1, nzs
                do j = 1, nys
                    u(0,j,k) = -bval_x(j+sdm%jsta-1,k+sdm%ksta-1,0)
                enddo
            enddo
        endif

        if(sdm%is_x1_boundary .eq. .true.) then
            do k = 1, nzs
                do j = 1, nys
                    u(nxs+1,j,k) = -bval_x(j+sdm%jsta-1,k+sdm%ksta-1,1)
                enddo
            enddo
        endif

        if(sdm%is_y0_boundary .eq. .true.) then
            do k = 1, nzs
                do i = 1, nxs
                    u(i,0,k) = -bval_y(i+sdm%ista-1,k+sdm%ksta-1,0)
                enddo
            enddo
        endif

        if(sdm%is_y1_boundary .eq. .true.) then
            do k = 1, nzs
                do i = 1, nxs
                    u(i,nys+1,k) = -bval_y(i+sdm%ista-1,k+sdm%ksta-1,1)
                enddo
            enddo
        endif

        if(sdm%is_z0_boundary .eq. .true.) then
            do j = 1, nys
                do i = 1, nxs
                    u(i,j,0) = -bval_z(i+sdm%ista-1,j+sdm%jsta-1,0)
                enddo
            enddo
        endif

        if(sdm%is_z1_boundary .eq. .true.) then
            do j = 1, nys
                do i = 1, nxs
                    u(i,j,nzs+1) = -bval_z(i+sdm%ista-1,j+sdm%jsta-1,1)
                enddo
            enddo
        endif
        if(myrank.eq.0) print '(a)', '[Boundary calculation] With MPI_IN_PLACE'
#else
        if(comm_boundary%mpi_comm.ne.MPI_COMM_NULL) then
            call MPI_Allreduce(bval_x, bval_xg, size(bval_x), MPI_REAL8, MPI_SUM, comm_boundary%mpi_comm, ierr)
            call MPI_Allreduce(bval_y, bval_yg, size(bval_y), MPI_REAL8, MPI_SUM, comm_boundary%mpi_comm, ierr)
            call MPI_Allreduce(bval_z, bval_zg, size(bval_z), MPI_REAL8, MPI_SUM, comm_boundary%mpi_comm, ierr)
        endif

        if(sdm%is_x0_boundary .eq. .true.) then
            do k = 1, nzs
                do j = 1, nys
                    u(0,j,k) = -bval_xg(j+sdm%jsta-1,k+sdm%ksta-1,0)
                enddo
            enddo
        endif

        if(sdm%is_x1_boundary .eq. .true.) then
            do k = 1, nzs
                do j = 1, nys
                    u(nxs+1,j,k) = -bval_xg(j+sdm%jsta-1,k+sdm%ksta-1,1)
                enddo
            enddo
        endif

        if(sdm%is_y0_boundary .eq. .true.) then
            do k = 1, nzs
                do i = 1, nxs
                    u(i,0,k) = -bval_yg(i+sdm%ista-1,k+sdm%ksta-1,0)
                enddo
            enddo
        endif

        if(sdm%is_y1_boundary .eq. .true.) then
            do k = 1, nzs
                do i = 1, nxs
                    u(i,nys+1,k) = -bval_yg(i+sdm%ista-1,k+sdm%ksta-1,1)
                enddo
            enddo
        endif

        if(sdm%is_z0_boundary .eq. .true.) then
            do j = 1, nys
                do i = 1, nxs
                    u(i,j,0) = -bval_zg(i+sdm%ista-1,j+sdm%jsta-1,0)
                enddo
            enddo
        endif

        if(sdm%is_z1_boundary .eq. .true.) then
            do j = 1, nys
                do i = 1, nxs
                    u(i,j,nzs+1) = -bval_zg(i+sdm%ista-1,j+sdm%jsta-1,1)
                enddo
            enddo
        endif
        if(myrank.eq.0) print '(a)', '[Boundary calculation] Without MPI_IN_PLACE'
        deallocate(bval_xg, bval_yg, bval_zg)

#endif
        deallocate(bval_x, bval_y, bval_z)

    end subroutine geometry_boundary_value_calculate


    subroutine geometry_boundary_normal_derivative(dudx, dudy, dudz, u, sdm)

        implicit none

        real(kind=8), intent(inout) :: dudx(0:,0:,0:), dudy(0:,0:,0:), dudz(0:,0:,0:)
        real(kind=8), intent(in)    :: u(0:,0:,0:)
        type(subdomain), intent(in) :: sdm

        integer(kind=4) :: i, j, k
        integer(kind=4) :: nx, ny, nz
        real(kind=8)    :: alpha, dxg1, dxg2, dyg1, dyg2, dzg1, dzg2

        nx = sdm%nx
        ny = sdm%ny
        nz = sdm%nz

        if(sdm%is_x0_boundary .eq. .true.) then
            dxg1 = sdm%dxg(2)
            dxg2 = sdm%dxg(3)
            alpha = dxg2 / dxg1
    
            do k = 1, nz
                do j = 1, ny
                    dudx(j,k,0) = -((u(2,j,k)-u(1,j,k))/dxg1*(2.0d0+alpha)/(1.0d0+alpha) - (u(3,j,k)-u(2,j,k))/dxg2/(1.0d0+alpha))
                enddo
            enddo
        endif

        if(sdm%is_x1_boundary .eq. .true.) then
            dxg1 = sdm%dxg(nx-1)
            dxg2 = sdm%dxg(nx)
            alpha = dxg1 / dxg2

            do k = 1, nz
                do j = 1, ny
                    dudx(j,k,1) = (u(nx,j,k)-u(nx-1,j,k))/dxg2*(2.0d0+alpha)/(1.0d0+alpha) - (u(nx-1,j,k)-u(nx-2,j,k))/dxg1/(1.0d0+alpha)
                enddo
            enddo
        endif

        if(sdm%is_y0_boundary .eq. .true.) then
            dyg1 = sdm%dyg(2)
            dyg2 = sdm%dyg(3)
            alpha = dyg2 / dyg1

            do k = 1, nz
                do i = 1, nx
                    dudy(i,k,0) = -((u(i,2,k)-u(i,1,k))/dyg1*(2.0d0+alpha)/(1.0d0+alpha) - (u(i,3,k)-u(i,2,k))/dyg2/(1.0d0+alpha))
                enddo
            enddo
        endif

        if(sdm%is_y1_boundary .eq. .true.) then
            dyg1 = sdm%dyg(ny-1)
            dyg2 = sdm%dyg(ny)
            alpha = dyg1 / dyg2

            do k = 1, nz
                do i = 1, nx
                    dudy(i,k,1) = (u(i,ny,k)-u(i,ny-1,k))/dyg2*(2.0d0+alpha)/(1.0d0+alpha) - (u(i,ny-1,k)-u(i,ny-2,k))/dyg1/(1.0d0+alpha)
                enddo
            enddo
        endif

        if(sdm%is_z0_boundary .eq. .true.) then
            dzg1 = sdm%dzg(2)
            dzg2 = sdm%dzg(3)
            alpha = dzg2 / dzg1

            do j = 1, ny
                do i = 1, nx
                    dudz(i,j,0) = -((u(i,j,2)-u(i,j,1))/dzg1*(2.0d0+alpha)/(1.0d0+alpha) - (u(i,j,3)-u(i,j,2))/dzg2/(1.0d0+alpha))
                enddo
            enddo
        endif

        if(sdm%is_z1_boundary .eq. .true.) then
            dzg1 = sdm%dzg(nz-1)
            dzg2 = sdm%dzg(nz)
            alpha = dzg1 / dzg2

            do j = 1, ny
                do i = 1, nx
                    dudz(i,j,1) = (u(i,j,nz)-u(i,j,nz-1))/dzg2*(2.0d0+alpha)/(1.0d0+alpha) - (u(i,j,nz-1)-u(i,j,nz-2))/dzg1/(1.0d0+alpha)
                enddo
            enddo
        endif

    end subroutine geometry_boundary_normal_derivative

    subroutine geometry_boundary_surface_integral(u, dudx, dudy, dudz, sdm, gdm)

        use mpi
        use mpi_topology, only : myrank

        implicit none

        real(kind=8), intent(inout)     :: u(0:,0:,0:)
        real(kind=8), intent(in)        :: dudx(0:,0:,0:), dudy(0:,0:,0:), dudz(0:,0:,0:)
        type(subdomain), intent(in)     :: sdm
        type(domain), intent(inout)     :: gdm

        real(kind=8), parameter         :: PI = 3.1415926535897932384626d0
        real(kind=8), allocatable, dimension(:,:,:)     :: bval_x, bval_y, bval_z, bval_xg, bval_yg, bval_zg
        integer(kind=4) :: i, j, k, ip, jp, kp, irp, jrp, krp
        integer(kind=4) :: nxs, nys, nzs, nxg, nyg, nzg
        integer(kind=4) :: ierr
        real(kind=8)    :: dr, drx, dry, drz
        character(256)  :: myfilename

        allocate(bval_x(gdm%ny,gdm%nz,0:1))
        allocate(bval_y(gdm%nx,gdm%nz,0:1))
        allocate(bval_z(gdm%nx,gdm%ny,0:1))

        allocate(bval_xg(gdm%ny,gdm%nz,0:1))
        allocate(bval_yg(gdm%nx,gdm%nz,0:1))
        allocate(bval_zg(gdm%nx,gdm%ny,0:1))

        bval_x = 0.0d0
        bval_y = 0.0d0
        bval_z = 0.0d0

        bval_xg = 0.0d0
        bval_yg = 0.0d0
        bval_zg = 0.0d0

        nxg = gdm%nx
        nyg = gdm%ny
        nzg = gdm%nz

        nxs = sdm%nx
        nys = sdm%ny
        nzs = sdm%nz

        ! Boundary treatment in x-direction
        do k = 1, nzg
            do j = 1, nyg
                ! Contributions from boundaries in x-direction
                do kp = 1, nzs
                    do jp = 1, nys
                        jrp = jp+sdm%jsta-1
                        krp = kp+sdm%ksta-1
                        ! Lower boundary in x-direction
                        drx = gdm%xg(0) - gdm%xg(1)
                        dry = gdm%yg(j) - gdm%yg(jrp)
                        drz = gdm%zg(k) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,0) = bval_x(j,k,0) + dudx(jp,kp,0) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        drx = gdm%xg(0) - gdm%xg(nxg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,0) = bval_x(j,k,0) + dudx(jp,kp,1) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        ! Upper boundary in x-direction
                        drx = gdm%xg(nxg+1) - gdm%xg(1)
                        dry = gdm%yg(j) - gdm%yg(jrp)
                        drz = gdm%zg(k) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,1) = bval_x(j,k,1) + dudx(jp,kp,0) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        drx = gdm%xg(nxg+1) - gdm%xg(nxg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,1) = bval_x(j,k,1) + dudx(jp,kp,1) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI
                    enddo
                enddo

                ! Contributions from boundaries in y-direction
                do kp = 1, nzs
                    do ip = 1, nxs
                        irp = ip+sdm%ista-1
                        krp = kp+sdm%ksta-1
                        ! Lower boundary in x-direction
                        drx = gdm%xg(0) - gdm%xg(irp)
                        dry = gdm%yg(j) - gdm%yg(1)
                        drz = gdm%zg(k) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,0) = bval_x(j,k,0) + dudy(ip,kp,0) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        dry = gdm%yg(j) - gdm%yg(nyg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,0) = bval_x(j,k,0) + dudy(ip,kp,1) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        ! Upper boundary in x-direction
                        drx = gdm%xg(nxg+1) - gdm%xg(irp)
                        dry = gdm%yg(j) - gdm%yg(1)
                        drz = gdm%zg(k) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,1) = bval_x(j,k,1) + dudy(ip,kp,0) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        dry = gdm%yg(j) - gdm%yg(nyg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,1) = bval_x(j,k,1) + dudy(ip,kp,1) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI
                    enddo
                enddo

                ! Contributions from boundaries in z-direction
                do jp = 1, nys
                    do ip = 1, nxs
                        irp = ip+sdm%ista-1
                        jrp = jp+sdm%jsta-1
                        ! Lower boundary in x-direction
                        drx = gdm%xg(0) - gdm%xg(irp)
                        dry = gdm%yg(j) - gdm%yg(jrp)
                        drz = gdm%zg(k) - gdm%zg(1)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,0) = bval_x(j,k,0) + dudz(ip,jp,0) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI

                        drz = gdm%zg(k) - gdm%zg(nzg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,0) = bval_x(j,k,0) + dudz(ip,jp,1) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI

                        ! Upper boundary in x-direction
                        drx = gdm%xg(nxg+1) - gdm%xg(irp)
                        dry = gdm%yg(j) - gdm%yg(jrp)
                        drz = gdm%zg(k) - gdm%zg(1)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,1) = bval_x(j,k,1) + dudz(ip,jp,0) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI

                        drz = gdm%zg(k) - gdm%zg(nzg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_x(j,k,1) = bval_x(j,k,1) + dudz(ip,jp,1) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI
                    enddo
                enddo
            enddo
        enddo

        ! Boundary treatment in y-direction
        do k = 1, nzg
            do i = 1, nxg
                ! Contributions from boundaries in x-direction
                do kp = 1, nzs
                    do jp = 1, nys
                        jrp = jp+sdm%jsta-1
                        krp = kp+sdm%ksta-1
                        ! Lower boundary in y-direction
                        drx = gdm%xg(i) - gdm%xg(1)
                        dry = gdm%yg(0) - gdm%yg(jrp)
                        drz = gdm%zg(k) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,0) = bval_y(i,k,0) + dudx(jp,kp,0) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        drx = gdm%xg(i) - gdm%xg(nxg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,0) = bval_y(i,k,0) + dudx(jp,kp,1) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        ! Upper boundary in y-direction
                        drx = gdm%xg(i) - gdm%xg(1)
                        dry = gdm%yg(nyg+1) - gdm%yg(jrp)
                        drz = gdm%zg(k) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,1) = bval_y(i,k,1) + dudx(jp,kp,0) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        drx = gdm%xg(i) - gdm%xg(nxg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,1) = bval_y(i,k,1) + dudx(jp,kp,1) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI
                    enddo
                enddo

                ! Contributions from boundaries in y-direction
                do kp = 1, nzs
                    do ip = 1, nxs
                        irp = ip+sdm%ista-1
                        krp = kp+sdm%ksta-1
                        ! Lower boundary in y-direction
                        drx = gdm%xg(i) - gdm%xg(irp)
                        dry = gdm%yg(0) - gdm%yg(1)
                        drz = gdm%zg(k) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,0) = bval_y(i,k,0) + dudy(ip,kp,0) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        dry = gdm%yg(0) - gdm%yg(nyg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,0) = bval_y(i,k,0) + dudy(ip,kp,1) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        ! Upper boundary in y-direction
                        drx = gdm%xg(i) - gdm%xg(irp)
                        dry = gdm%yg(nyg+1) - gdm%yg(1)
                        drz = gdm%zg(k) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,1) = bval_y(i,k,1) + dudy(ip,kp,0) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        dry = gdm%yg(nyg+1) - gdm%yg(nyg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,1) = bval_y(i,k,1) + dudy(ip,kp,1) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI
                    enddo
                enddo

                ! Contributions from boundaries in z-direction
                do jp = 1, nys
                    do ip = 1, nxs
                        irp = ip+sdm%ista-1
                        jrp = jp+sdm%jsta-1
                        ! Lower boundary in y-direction
                        drx = gdm%xg(i) - gdm%xg(irp)
                        dry = gdm%yg(0) - gdm%yg(jrp)
                        drz = gdm%zg(k) - gdm%zg(1)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,0) = bval_y(i,k,0) + dudz(ip,jp,0) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI

                        drz = gdm%zg(k) - gdm%zg(nzg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,0) = bval_y(i,k,0) + dudz(ip,jp,1) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI

                        ! Upper boundary in y-direction
                        drx = gdm%xg(i) - gdm%xg(irp)
                        dry = gdm%yg(nyg+1) - gdm%yg(jrp)
                        drz = gdm%zg(k) - gdm%zg(1)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,1) = bval_y(i,k,1) + dudz(ip,jp,0) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI

                        drz = gdm%zg(k) - gdm%zg(nzg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_y(i,k,1) = bval_y(i,k,1) + dudz(ip,jp,1) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI
                    enddo
                enddo
            enddo
        enddo

        ! Boundary treatment in z-direction
        do j = 1, nyg
            do i = 1, nxg
                ! Contributions from boundaries in x-direction
                do kp = 1, nzs
                    do jp = 1, nys
                        jrp = jp+sdm%jsta-1
                        krp = kp+sdm%ksta-1
                        ! Lower boundary in z-direction
                        drx = gdm%xg(i) - gdm%xg(1)
                        dry = gdm%yg(j) - gdm%yg(jrp)
                        drz = gdm%zg(0) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,0) = bval_z(i,j,0) + dudx(jp,kp,0) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        drx = gdm%xg(i) - gdm%xg(nxg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,0) = bval_z(i,j,0) + dudx(jp,kp,1) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        ! Upper boundary in z-direction
                        drx = gdm%xg(i) - gdm%xg(1)
                        dry = gdm%yg(j) - gdm%yg(jrp)
                        drz = gdm%zg(nzg+1) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,1) = bval_z(i,j,1) + dudx(jp,kp,0) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        drx = gdm%xg(i) - gdm%xg(nxg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,1) = bval_z(i,j,1) + dudx(jp,kp,1) * sdm%dym(jp) * sdm%dzm(kp) / dr / 4.0d0 / PI
                    enddo
                enddo

                ! Contributions from boundaries in y-direction
                do kp = 1, nzs
                    do ip = 1, nxs
                        irp = ip+sdm%ista-1
                        krp = kp+sdm%ksta-1
                        ! Lower boundary in z-direction
                        drx = gdm%xg(i) - gdm%xg(irp)
                        dry = gdm%yg(j) - gdm%yg(1)
                        drz = gdm%zg(0) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,0) = bval_z(i,j,0) + dudy(ip,kp,0) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        dry = gdm%yg(j) - gdm%yg(nyg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,0) = bval_z(i,j,0) + dudy(ip,kp,1) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        ! Upper boundary in z-direction
                        drx = gdm%xg(i) - gdm%xg(irp)
                        dry = gdm%yg(j) - gdm%yg(1)
                        drz = gdm%zg(nzg+1) - gdm%zg(krp)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,1) = bval_z(i,j,1) + dudy(ip,kp,0) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI

                        dry = gdm%yg(j) - gdm%yg(nyg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,1) = bval_z(i,j,1) + dudy(ip,kp,1) * sdm%dxm(ip) * sdm%dzm(kp) / dr / 4.0d0 / PI
                    enddo
                enddo

                ! Contributions from boundaries in z-direction
                do jp = 1, nys
                    do ip = 1, nxs
                        irp = ip+sdm%ista-1
                        jrp = jp+sdm%jsta-1
                        ! Lower boundary in z-direction
                        drx = gdm%xg(i) - gdm%xg(irp)
                        dry = gdm%yg(j) - gdm%yg(jrp)
                        drz = gdm%zg(0) - gdm%zg(1)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,0) = bval_z(i,j,0) + dudz(ip,jp,0) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI

                        drz = gdm%zg(0) - gdm%zg(nzg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,0) = bval_z(i,j,0) + dudz(ip,jp,1) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI

                        ! Upper boundary in z-direction
                        drx = gdm%xg(i) - gdm%xg(irp)
                        dry = gdm%yg(j) - gdm%yg(jrp)
                        drz = gdm%zg(nzg+1) - gdm%zg(1)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,1) = bval_z(i,j,1) + dudz(ip,jp,0) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI

                        drz = gdm%zg(nzg+1) - gdm%zg(nzg)
                        dr = sqrt(drx*drx+dry*dry+drz*drz)
                        bval_z(i,j,1) = bval_z(i,j,1) + dudz(ip,jp,1) * sdm%dxm(ip) * sdm%dym(jp) / dr / 4.0d0 / PI
                    enddo
                enddo
            enddo
        enddo
        call MPI_Allreduce(bval_x, bval_xg, size(bval_x), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(bval_y, bval_yg, size(bval_y), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(bval_z, bval_zg, size(bval_z), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        write(myfilename,'(a,i0.3)') './result/bval.',myrank
        open(myrank, file=myfilename, form='formatted')

        write(myrank,'(a)') 'bval_xg'
        do k=1, nzg
            do j=1, nyg
                write(myrank,102) gdm%yg(j), gdm%zg(k), bval_xg(j,k,0), bval_xg(j,k,1)
            enddo
        enddo

        write(myrank,'(a)') 'bval_yg'
        do k=1, nzg
            do i=1, nxg
                write(myrank,102) gdm%xg(i), gdm%zg(k), bval_yg(i,k,0), bval_yg(i,k,1)
            enddo
        enddo
    
        write(myrank,'(a)') 'bval_zg'
        do j=1, nyg
            do i=1, nxg
                write(myrank,102) gdm%xg(i), gdm%yg(j), bval_zg(i,j,0), bval_zg(i,j,1)
            enddo
        enddo
        102   format(2(f9.5,x),2(e15.7,x))
    
        close(myrank) 

        if(sdm%is_x0_boundary .eq. .true.) then
            do k = 1, nzs
                do j = 1, nys
                    u(0,j,k) = -bval_xg(j+sdm%jsta-1,k+sdm%ksta-1,0)
                enddo
            enddo
        endif

        if(sdm%is_x1_boundary .eq. .true.) then
            do k = 1, nzs
                do j = 1, nys
                    u(nxs+1,j,k) = -bval_xg(j+sdm%jsta-1,k+sdm%ksta-1,1)
                enddo
            enddo
        endif

        if(sdm%is_y0_boundary .eq. .true.) then
            do k = 1, nzs
                do i = 1, nxs
                    u(i,0,k) = -bval_yg(i+sdm%ista-1,k+sdm%ksta-1,0)
                enddo
            enddo
        endif

        if(sdm%is_y1_boundary .eq. .true.) then
            do k = 1, nzs
                do i = 1, nxs
                    u(i,nys+1,k) = -bval_yg(i+sdm%ista-1,k+sdm%ksta-1,1)
                enddo
            enddo
        endif

        if(sdm%is_z0_boundary .eq. .true.) then
            do j = 1, nys
                do i = 1, nxs
                    u(i,j,0) = -bval_zg(i+sdm%ista-1,j+sdm%jsta-1,0)
                enddo
            enddo
        endif

        if(sdm%is_z1_boundary .eq. .true.) then
            do j = 1, nys
                do i = 1, nxs
                    u(i,j,nzs+1) = -bval_zg(i+sdm%ista-1,j+sdm%jsta-1,1)
                enddo
            enddo
        endif

        deallocate(bval_x, bval_y, bval_z)
        deallocate(bval_xg, bval_yg, bval_zg)
    end subroutine geometry_boundary_surface_integral

end module geometry
