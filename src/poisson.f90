program main

    use mpi
    use timer
    use mpi_topology
    use multigrid
    use geometry
    use matrix
    use, intrinsic :: ISO_C_BINDING

    implicit none

    real(kind=8), parameter         :: PI = 3.1415926535897932384626d0

    integer(kind=4)                 :: npx, npy, npz
    integer(kind=4)                 :: nx, ny, nz
    real(kind=8)                    :: ox, oy, oz
    real(kind=8)                    :: lx, ly, lz
    real(kind=8)                    :: ax, ay, az
    real(kind=8)                    :: alpha_x, alpha_y, alpha_z

    type(domain)                    :: g_domain
    type(subdomain)                 :: s_domain
    type(matrix_poisson)            :: a_poisson

    real(kind=8), allocatable       :: ref_sub(:,:,:)
    real(kind=8), allocatable       :: u_x, u_y, u_z

    real(kind=8)                    :: rms, rms_local

    ! aggretation method 
    ! 0 : no aggregation
    ! 1 : single aggregation
    ! 2 : adaptive aggregation (not implemented yet)
    integer(kind=4) :: i, j, k, maxiteration, number_of_vcycles, number_of_levels, aggregation_method, single_aggregation_level, adaptive_aggregation_level(3)
    real(kind=8)    :: tolerance, t0, omega_sor
    ! character(256)  :: myfilename

    integer(kind=4) :: ierr

    call MPI_Init(ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)

    call timer_init

    namelist /meshes/ nx, ny, nz
    namelist /origin/ ox, oy, oz
    namelist /length/ lx, ly, lz
    namelist /mesh_stretch/ ax, ay, az
    namelist /procs/ npx, npy, npz
    namelist /control/ maxiteration, tolerance, number_of_vcycles, number_of_levels, aggregation_method, single_aggregation_level, adaptive_aggregation_level, omega_sor
    namelist /coefficients/ alpha_x, alpha_y, alpha_z

    open (unit = 1, file = "PARA_INPUT.inp")
    read (1, meshes)
    read (1, origin)
    read (1, length)
    read (1, mesh_stretch)
    read (1, procs)
    read (1, control)
    read (1, coefficients)
    close (1)

    np_dim(0)=npx
    np_dim(1)=npy
    np_dim(2)=npz
    period(0)=.false.; period(1)=.false.; period(2)=.false.

    call mpi_topology_create()

    call geometry_domain_create(g_domain, nx, ny, nz, ox, oy, oz, lx, ly, lz, ax, ay, az, period)
    if(myrank.eq.0) print '(a)', '[Poisson] Geometry and matrix size initialized.'

    call geometry_subdomain_create(s_domain, g_domain)
    call geometry_subdomain_ddt_create(s_domain)

    allocate( ref_sub(0:s_domain%nx+1, 0:s_domain%ny+1, 0:s_domain%nz+1) )
    ref_sub(:,:,:) = 0.0d0

    ! Problem 1 : cosine function over -l/2 <= x,y,z <=l/2
    do k = 1, s_domain%nz
        do j = 1, s_domain%ny
            do i = 1, s_domain%nx
                s_domain%x(i, j, k) = 0.0d0
                s_domain%b(i, j, k) =-cos(s_domain%xg(i)/lx*PI) &
                                    *cos(s_domain%yg(j)/ly*PI) &
                                    *cos(s_domain%zg(k)/lz*PI) &
                                    *PI*PI*(1.0d0/(lx*lx)+1.0d0/(ly*ly)+1.0d0/(lz*lz))
                ref_sub(i, j, k) =  cos(s_domain%xg(i)/lx*PI)*cos(s_domain%yg(j)/ly*PI)*cos(s_domain%zg(k)/lz*PI)
            enddo
        enddo
    enddo

    ! Problem 2 : delta function at x = y = z = 0.0
    ! do k = 1, s_domain%nz
    !     kg = k + s_domain%ksta - 1
    !     do j = 1, s_domain%ny
    !         jg = j + s_domain%jsta - 1
    !         do i = 1, s_domain%nx
    !             ig = i + s_domain%ista - 1
    !             if( (ig.eq.g_domain%nx/4).and.(jg.eq.g_domain%ny/4).and.(kg.eq.g_domain%nz/4) ) then
    !                 s_domain%b(i,j,k) = 4.0d0 *PI / 8.0d0 / (s_domain%dxm(i)*s_domain%dym(j)*s_domain%dzm(k) )
    !                 print *, '1', myrank, i, j, k, s_domain%b(i,j,k)
    !             endif
    !             if( (ig.eq.g_domain%nx/4+1).and.(jg.eq.g_domain%ny/4).and.(kg.eq.g_domain%nz/4) ) then
    !                 s_domain%b(i,j,k) = 4.0d0 *PI / 8.0d0  / (s_domain%dxm(i)*s_domain%dym(j)*s_domain%dzm(k) )
    !                 print *, '2', myrank, i, j, k, s_domain%b(i,j,k)
    !             endif
    !             if( (ig.eq.g_domain%nx/4).and.(jg.eq.g_domain%ny/4+1).and.(kg.eq.g_domain%nz/4) ) then
    !                 s_domain%b(i,j,k) = 4.0d0 *PI / 8.0d0  / (s_domain%dxm(i)*s_domain%dym(j)*s_domain%dzm(k) )
    !                 print *, '3', myrank, i, j, k, s_domain%b(i,j,k)
    !             endif
    !             if( (ig.eq.g_domain%nx/4+1).and.(jg.eq.g_domain%ny/4+1).and.(kg.eq.g_domain%nz/4) ) then
    !                 s_domain%b(i,j,k) = 4.0d0 *PI  / 8.0d0 / (s_domain%dxm(i)*s_domain%dym(j)*s_domain%dzm(k) )
    !                 print *, '4', myrank, i, j, k, s_domain%b(i,j,k)
    !             endif
    !             if( (ig.eq.g_domain%nx/4).and.(jg.eq.g_domain%ny/4).and.(kg.eq.g_domain%nz/4+1) ) then
    !                 s_domain%b(i,j,k) = 4.0d0 *PI / 8.0d0  / (s_domain%dxm(i)*s_domain%dym(j)*s_domain%dzm(k) )
    !                 print *, '5', myrank, i, j, k, s_domain%b(i,j,k)
    !             endif
    !             if( (ig.eq.g_domain%nx/4+1).and.(jg.eq.g_domain%ny/4).and.(kg.eq.g_domain%nz/4+1) ) then
    !                 s_domain%b(i,j,k) = 4.0d0 *PI / 8.0d0  / (s_domain%dxm(i)*s_domain%dym(j)*s_domain%dzm(k) )
    !                 print *, '6', myrank, i, j, k, s_domain%b(i,j,k)
    !             endif
    !             if( (ig.eq.g_domain%nx/4).and.(jg.eq.g_domain%ny/4+1).and.(kg.eq.g_domain%nz/4+1) ) then
    !                 s_domain%b(i,j,k) = 4.0d0 *PI / 8.0d0  / (s_domain%dxm(i)*s_domain%dym(j)*s_domain%dzm(k) )
    !                 print *, '7', myrank, i, j, k, s_domain%b(i,j,k)
    !             endif
    !             if( (ig.eq.g_domain%nx/4+1).and.(jg.eq.g_domain%ny/4+1).and.(kg.eq.g_domain%nz/4+1) ) then
    !                 s_domain%b(i,j,k) = 4.0d0 *PI / 8.0d0  / (s_domain%dxm(i)*s_domain%dym(j)*s_domain%dzm(k) )
    !                 print *, '8', myrank, i, j, k, s_domain%b(i,j,k)
    !             endif
    !         enddo
    !     enddo
    ! enddo
    ! do k = 0, s_domain%nz+1
    !     do j = 0, s_domain%ny+1
    !         do i = 0, s_domain%nx+1
    !             r01=sqrt(s_domain%xg(i)**2+s_domain%yg(j)**2+s_domain%zg(k)**2)
    !             ref_sub(i, j, k) =  - 1.0d0 / r01
    !         enddo
    !     enddo
    ! enddo

    ! Problem 3 : Exponential function, -1.0 <= x <= 1.0
    ! do k = 0, s_domain%nz+1
    !     do j = 0, s_domain%ny+1
    !         do i = 0, s_domain%nx+1
    !             u_x = ( 1.0d0 - exp( alpha_x*((s_domain%xg(i))**2 - 1.0d0) ) ) / ( 1.0d0 - exp(-alpha_x) )
    !             u_y = ( 1.0d0 - exp( alpha_y*((s_domain%yg(j))**2 - 1.0d0) ) ) / ( 1.0d0 - exp(-alpha_y) )
    !             u_z = ( 1.0d0 - exp( alpha_z*((s_domain%zg(k))**2 - 1.0d0) ) ) / ( 1.0d0 - exp(-alpha_z) )
                
    !             s_domain%x(i, j, k) = 0.0d0
    !             s_domain%b(i, j, k) = ( 2.0d0*alpha_x + 4.0d0*alpha_x**2 * (s_domain%xg(i))**2 ) * (u_x - 1.0d0/(1.0d0 - exp(-alpha_x))) * u_y * u_z &
    !                                 + ( 2.0d0*alpha_y + 4.0d0*alpha_y**2 * (s_domain%yg(j))**2 ) * (u_y - 1.0d0/(1.0d0 - exp(-alpha_y))) * u_z * u_x &
    !                                 + ( 2.0d0*alpha_z + 4.0d0*alpha_z**2 * (s_domain%zg(k))**2 ) * (u_z - 1.0d0/(1.0d0 - exp(-alpha_z))) * u_x * u_y
    !             ref_sub(i, j, k) =  u_x * u_y * u_z
    !         enddo
    !     enddo
    ! enddo

    if(myrank.eq.0) print '(a)', '[Poisson] Geometry and rhs constructed.'

    call matrix_poisson_create(a_poisson, s_domain)
    if(myrank.eq.0) print '(a)', '[Poisson] Poisson matrix constructed.'

    if(myrank.eq.0) print '(a)', '[Poisson] Start solving equations.'

    t0=MPI_Wtime()

    call multigrid_create(s_domain, number_of_levels, number_of_vcycles, aggregation_method, single_aggregation_level, adaptive_aggregation_level)
    call multigrid_solve_vcycle(s_domain%x, s_domain%r, a_poisson, s_domain%b, s_domain, maxiteration, tolerance, omega_sor)
    call multigrid_destroy
    rms = 0.0d0
    rms_local = 0.0d0
    do k=1, s_domain%nz
        do j=1, s_domain%ny
            do i=1, s_domain%nx
                rms_local = rms_local + (s_domain%x(i,j,k)-ref_sub(i,j,k))**2
            enddo
        enddo
    enddo
    call MPI_Allreduce(rms_local, rms, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    if(myrank.eq.0) print '(a,e20.10)','[Poisson] RMS = ',rms/g_domain%nx/g_domain%ny/g_domain%nz
    if(myrank.eq.0) print '(a,f12.6)', '[Poisson] Solution obtained. Execution time = ',MPI_Wtime()-t0

    call timer_reduction
    call timer_output

    ! call file_write
    ! call binary_file_check

    if(myrank.eq.0) print '(a)', '[Poisson] Solution printed.'

    deallocate(ref_sub)

    call matrix_poisson_destroy(a_poisson)
    call mpi_topology_destroy
    call geometry_subdomain_ddt_destroy(s_domain)
    call geometry_subdomain_destroy(s_domain)
    call geometry_domain_destroy(g_domain)

    if(myrank.eq.0) print '(a)', '[Poisson] Memory deallocated.'

    call MPI_Finalize(ierr)

    contains

    subroutine file_write

        use mpi
        use mpi_topology

        implicit none

        !---- MPI write ----
        real(kind=8), allocatable, dimension(:, :, :)   :: sol_all, global_x
        integer(kind=4), allocatable, dimension(:,:)    :: cart_coord
        integer(kind=4), allocatable, dimension(:)      :: n1subm_cnt,  n2subm_cnt,  n3subm_cnt
        integer(kind=4), allocatable, dimension(:)      :: n1subm_disp, n2subm_disp, n3subm_disp
        integer(kind=4), allocatable, dimension(:)      :: ddtype_data_write_recv, request_recv

        integer(kind=4) :: ddtype_temp, ddtype_temp2, ddtype_data_write_send, ddtype_write_gatherv, ddtype_file, request_send
        integer(kind=4) :: ddtype_temp3, ddtype_write_gatherv_halo
        integer(kind=4) :: sizes(0:2), subsizes(0:2), starts(0:2)
        integer(kind=4) :: i, j, k, ierr, r8size
        integer(kind=4) :: status(MPI_STATUS_SIZE), filep

        integer(kind=MPI_OFFSET_KIND) :: disp
        integer(kind=MPI_ADDRESS_KIND):: extent, lb
        integer(kind=4), allocatable, dimension(:)   :: cnt_gatherv
        integer(kind=4), allocatable, dimension(:)   :: disps_gatherv
        ! ascii file write

        allocate(n1subm_cnt(0:comm_1d_x%nprocs-1), n1subm_disp(0:comm_1d_x%nprocs-1))
        allocate(n2subm_cnt(0:comm_1d_y%nprocs-1), n2subm_disp(0:comm_1d_y%nprocs-1))
        allocate(n3subm_cnt(0:comm_1d_z%nprocs-1), n3subm_disp(0:comm_1d_z%nprocs-1))
        allocate(ddtype_data_write_recv(0:nprocs-1))
        allocate(request_recv(0:nprocs-1))
        allocate(disps_gatherv(0:nprocs-1), cnt_gatherv(0:nprocs-1))
        allocate(sol_all(1:g_domain%nx,1:g_domain%ny,1:g_domain%nz))
        allocate(global_x(0:g_domain%nx+1,0:g_domain%ny+1,0:g_domain%nz+1))

        sol_all(:,:,:) = 0.0d0
        allocate( cart_coord(0:2,0:nprocs-1) )
        do i = 0, nprocs-1
            call MPI_Cart_coords(mpi_world_cart, i, 3, cart_coord(:,i), ierr )
        enddo
    
        call MPI_Allgather(s_domain%nx, 1, MPI_INTEGER, n1subm_cnt,1, MPI_INTEGER, comm_1d_x%mpi_comm, ierr)
        call MPI_Allgather(s_domain%ny, 1, MPI_INTEGER, n2subm_cnt,1, MPI_INTEGER, comm_1d_y%mpi_comm, ierr)
        call MPI_Allgather(s_domain%nz, 1, MPI_INTEGER, n3subm_cnt,1, MPI_INTEGER, comm_1d_z%mpi_comm, ierr)
    
        n1subm_disp(0) = 0
        do i = 1, comm_1d_x%nprocs-1
            n1subm_disp(i)=sum(n1subm_cnt(0:i-1))
        enddo
    
        n2subm_disp(0) = 0
        do i = 1, comm_1d_y%nprocs-1
            n2subm_disp(i)=sum(n2subm_cnt(0:i-1))
        enddo
    
        n3subm_disp(0) = 0
        do i = 1, comm_1d_z%nprocs-1
            n3subm_disp(i)=sum(n3subm_cnt(0:i-1))
        enddo
    
        ! Make data types for MPI_Send/Recv
        ! ddtype_send is common
        sizes    = (/s_domain%nx+2,s_domain%ny+2,s_domain%nz+2/)
        subsizes = (/s_domain%nx,s_domain%ny,s_domain%nz/)
        starts   = (/1,      1,      1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, ddtype_temp, ierr)
        CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
        lb = 0
        extent = r8size * s_domain%nx
        call MPI_Type_create_resized(ddtype_temp, lb, extent, ddtype_data_write_send, ierr)
        call MPI_Type_commit( ddtype_data_write_send, ierr)

        ! ddtype_recv for MPI_Send/Recv
        do i = 0, nprocs-1
            sizes    = (/g_domain%nx, g_domain%ny, g_domain%nz/)
            subsizes = (/n1subm_cnt(cart_coord(0,i)), n2subm_cnt(cart_coord(1,i)), n3subm_cnt(cart_coord(2,i))/)
            starts   = (/n1subm_disp(cart_coord(0,i)), n2subm_disp(cart_coord(1,i)), n3subm_disp(cart_coord(2,i))/)
    
            call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_REAL8, ddtype_data_write_recv(i), ierr)
            call MPI_Type_commit( ddtype_data_write_recv(i), ierr)
        enddo

        ! Communication using MPI_Send/Recv
        call MPI_Isend(s_domain%x(0,0,0), 1, ddtype_data_write_send, 0, 101, MPI_COMM_WORLD, request_send, ierr)
    
        if(myrank == 0 ) then
            do i = 0, nprocs-1
                call MPI_Irecv(sol_all(1,1,1), 1, ddtype_data_write_recv(i), i, 101, MPI_COMM_WORLD, request_recv(i), ierr)
            enddo
        endif
    
        call MPI_Wait(request_send, status, ierr)
        if(myrank == 0 ) then
            call MPI_Waitall(nprocs, request_recv, MPI_STATUSES_IGNORE, ierr)
        endif

        ! File write after MPI_Send/Recv
        if(myrank == 0) then
            open(unit=myrank,form="formatted",file="solution_all.dat")
            do k = 1, g_domain%nz
                do j = 1, g_domain%ny
                    do i = 1, g_domain%nx
                        write(myrank,'(3D20.10,1E30.20)') g_domain%xg(i),g_domain%yg(j),g_domain%zg(k),sol_all(i,j,k)
                    enddo
                enddo
            enddo
            close(myrank)

            open(unit=myrank+1,form="unformatted",file="solution_all.bin")
            write(myrank+1) sol_all
            close(myrank+1)
        endif

        ! ddtype_recv for MPI_Gatherv

        sol_all = 0.0d0

        sizes    = (/g_domain%nx, g_domain%ny, g_domain%nz/)
        subsizes = (/s_domain%nx, s_domain%ny, s_domain%nz/)
        starts   = (/0, 0, 0/)

        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, ddtype_temp2, ierr)
        CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
        lb = 0
        extent = r8size
        call MPI_Type_create_resized(ddtype_temp2, lb, extent, ddtype_write_gatherv, ierr)
        call MPI_Type_commit(ddtype_write_gatherv, ierr)

        do i = 0, nprocs-1
            cnt_gatherv(i) = 1
            disps_gatherv(i) = n1subm_disp(cart_coord(0,i))  &
                             + n2subm_disp(cart_coord(1,i)) * g_domain%nx  &
                             + n3subm_disp(cart_coord(2,i)) * g_domain%nx * g_domain%ny 
        enddo

        call MPI_Gatherv(s_domain%x, 1, s_domain%ddt_inner_domain, sol_all(1,1,1), cnt_gatherv, disps_gatherv, ddtype_write_gatherv, 0, MPI_COMM_WORLD, ierr)

        ! File write after MPI_Gatherv
        if(myrank == 0) then
            open(unit=myrank,form="formatted",file="solution_all_gatherv.dat")
            do k = 1, g_domain%nz
                do j = 1, g_domain%ny
                    do i = 1, g_domain%nx
                        write(myrank,'(3D20.10,1E30.20)') g_domain%xg(i),g_domain%yg(j),g_domain%zg(k),sol_all(i,j,k)
                    enddo
                enddo
            enddo
            close(myrank)

            open(unit=myrank+1,form="unformatted",file="solution_all_gatherv.bin")
            write(myrank+1) sol_all
            close(myrank+1)
        endif

        ! ddtype_recv for MPI_Gatherv with halo cells

        global_x = 0.0d0

        sizes    = (/g_domain%nx+2, g_domain%ny+2, g_domain%nz+2/)
        subsizes = (/s_domain%nx, s_domain%ny, s_domain%nz/)
        starts   = (/1, 1, 1/)

        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, ddtype_temp3, ierr)
        CALL MPI_Type_size(MPI_REAL8, r8size, ierr)
        lb = 0
        extent = r8size
        call MPI_Type_create_resized(ddtype_temp3, lb, extent, ddtype_write_gatherv_halo, ierr)
        call MPI_Type_commit(ddtype_write_gatherv_halo, ierr)

        do i = 0, nprocs-1
            cnt_gatherv(i) = 1
            disps_gatherv(i) = n1subm_disp(cart_coord(0,i))  &
                             + n2subm_disp(cart_coord(1,i)) * (g_domain%nx + 2)  &
                             + n3subm_disp(cart_coord(2,i)) * (g_domain%nx + 2) * (g_domain%ny + 2)
        enddo

        call MPI_Gatherv(s_domain%x, 1, s_domain%ddt_inner_domain, global_x, cnt_gatherv, disps_gatherv, ddtype_write_gatherv_halo, 0, MPI_COMM_WORLD, ierr)

        ! File write after MPI_Gatherv with halo cells
        if(myrank == 0) then
            open(unit=myrank,form="formatted",file="solution_all_gatherv_halo.dat")
            do k = 1, g_domain%nz
                do j = 1, g_domain%ny
                    do i = 1, g_domain%nx
                        write(myrank,'(3D20.10,1E30.20)') g_domain%xg(i),g_domain%yg(j),g_domain%zg(k),global_x(i,j,k)
                    enddo
                enddo
            enddo
            close(myrank)

            open(unit=myrank+1,form="unformatted",file="solution_all_gatherv_halo.bin")
            write(myrank+1) sol_all
            close(myrank+1)
        endif

        ! ddtype_recv for MPI_Send/Recv
        do i = 0, nprocs-1
            sizes    = (/g_domain%nx, g_domain%ny, g_domain%nz/)
            subsizes = (/n1subm_cnt(cart_coord(0,i)), n2subm_cnt(cart_coord(1,i)), n3subm_cnt(cart_coord(2,i))/)
            starts   = (/n1subm_disp(cart_coord(0,i)), n2subm_disp(cart_coord(1,i)), n3subm_disp(cart_coord(2,i))/)
    
            call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_REAL8, ddtype_data_write_recv(i), ierr)
            call MPI_Type_commit( ddtype_data_write_recv(i), ierr)
        enddo

        ! For binary file write
        ! ddtype for MPI file write
        sizes    = (/g_domain%nx, g_domain%ny, g_domain%nz/)
        subsizes = (/s_domain%nx, s_domain%ny, s_domain%nz/)
        starts   = (/s_domain%ista-1, s_domain%jsta-1, s_domain%ksta-1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_REAL8, ddtype_file, ierr)
        call MPI_Type_commit( ddtype_file, ierr)

        ! Set view and file write
        disp = 0
        call MPI_File_open(MPI_COMM_WORLD, 'solution_all_MPI_IO.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_REAL8, ddtype_file, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_write(filep, s_domain%x(0,0,0), 1, ddtype_data_write_send, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)
        
        call MPI_Type_free(ddtype_data_write_send, ierr)
        do i = 0, nprocs-1
            call MPI_Type_free(ddtype_data_write_recv(i), ierr)
        enddo
        call MPI_Type_free(ddtype_write_gatherv_halo, ierr)
        call MPI_Type_free(ddtype_write_gatherv, ierr)
        call MPI_Type_free(ddtype_file, ierr)
        deallocate(sol_all)
        deallocate(global_x)
        deallocate(n1subm_cnt, n1subm_disp)
        deallocate(n2subm_cnt, n2subm_disp)
        deallocate(n3subm_cnt, n3subm_disp)
        deallocate(ddtype_data_write_recv, request_recv)
        deallocate(cnt_gatherv, disps_gatherv)
        deallocate( cart_coord )

    end subroutine file_write

    subroutine binary_file_check

        use mpi
        use mpi_topology
        implicit none

        integer(kind=4) :: filep, ierr
        integer(kind=4) :: status(MPI_STATUS_SIZE)
        real(kind=8), allocatable :: temp1(:,:,:), temp2(:,:,:)

        if(myrank.eq.0) then

            allocate(temp1(1:g_domain%nx,1:g_domain%ny,1:g_domain%nz))
            open(unit=1,form='unformatted',file="solution_all.bin")
            read(1) temp1
            close(1)

            allocate(temp2(1:g_domain%nx,1:g_domain%ny,1:g_domain%nz))
            call mpi_file_open(MPI_COMM_SELF, 'solution_all_MPI_IO.bin', MPI_MODE_RDONLY, &
            MPI_INFO_NULL, filep, ierr)
            call mpi_file_read(filep, temp2, g_domain%nx*g_domain%ny*g_domain%nz, MPI_REAL8, status, ierr)
            call mpi_file_close(filep,ierr)

            do k = 1, g_domain%nz
                do j = 1, g_domain%ny
                    do i = 1, g_domain%nx
                        if(temp1(i,j,k).ne.temp2(i,j,k)) then
                            print *, i, j, k, temp1(i,j,k), temp2(i,j,k)
                        endif
                    enddo
                enddo
            enddo

            deallocate(temp1, temp2)
        endif
        
    end subroutine binary_file_check

end program main