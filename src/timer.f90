module timer

    use mpi

    implicit none

    real(kind=8), private       :: t_zero, t_curr, t_comm_zero, t_comm_curr
    real(kind=8), public        :: t_array(64), t_array_reduce(64)
    character(len=64)  :: t_str(64)

    contains

    subroutine  timer_init

        t_array(:) = 0.0d0
        t_array_reduce(:) = 0.0d0

        t_str(:)    = 'null'

        t_str(1)    = 'Total : V-cycle_solver'
        t_str(2)    = 'Comp : LV0_Iteration'
        t_str(3)    = 'Comp : LV0_Residual'
        t_str(4)    = 'Comp : LV0_Res_dot_product'
        t_str(5)    = 'Comp : LV0_Restriction'
        t_str(6)    = 'Comp : LV0_Prolongation'
        t_str(7)    = 'Comp : LVs_Iterator'
        t_str(8)    = 'Comp : LVs_Residual'
        t_str(9)    = 'Comp : LVs_Restriction'
        t_str(10)   = 'Comp : LVs_Prolongation'
        t_str(11)   = 'Comp : Coarset_LV_Solve'
        t_str(12)   = 'Comp : LV_aggr1_Iteration'
        t_str(13)   = 'Comp : LV_aggr1_Residual'
        t_str(14)   = 'Comp : LV_aggr1_Restriction'
        t_str(15)   = 'Comp : LV_aggr1_Prolongation'
        t_str(16)   = 'Comp : LV_aggr2_Iteration'
        t_str(17)   = 'Comp : LV_aggr2_Residual'
        t_str(18)   = 'Comp : LV_aggr2_Restriction'
        t_str(19)   = 'Comp : LV_aggr2_Prolongation'
        t_str(20)   = 'Comp : LV_aggr3_Iteration'
        t_str(21)   = 'Comp : LV_aggr3_Residual'
        t_str(22)   = 'Comp : LV_aggr3_Restriction'
        t_str(23)   = 'Comp : LV_aggr3_Prolongation'

        t_str(24)   = 'Comm : LVs_PL_Update_halocells'
        t_str(25)   = 'none : None'
        t_str(26)   = 'Comm : MPI_Gather'
        t_str(27)   = 'Comm : MPI_Scatter'
        t_str(28)   = 'Comm : MPI_comm_setup'
        t_str(29)   = 'Others : LV0_Residual_mem_alloc'
        t_str(30)   = 'Others : etc'

        t_str(31)   = 'Comm : LVs_MPI_Allreduce'
        t_str(32)   = 'Comm : LVc_MPI_Allreduce'
        t_str(33)   = 'Comm : LVs_MPI_SendRecv'
        t_str(34)   = 'Comm : LVc_MPI_SendRecv'

    end subroutine timer_init

    subroutine timer_stamp0

        t_zero = MPI_Wtime()

    end subroutine timer_stamp0

    subroutine timer_stamp(timer_id)

        integer(kind=4), intent(in) :: timer_id

        t_curr = MPI_Wtime()
        t_array(timer_id) = t_array(timer_id) + t_curr - t_zero
        t_zero = t_curr

    end subroutine timer_stamp

    subroutine timer_comm_stamp0

        t_comm_zero = MPI_Wtime()

    end subroutine timer_comm_stamp0

    subroutine timer_comm_stamp(timer_id)

        integer(kind=4), intent(in) :: timer_id

        t_comm_curr = MPI_Wtime()
        t_array(timer_id) = t_array(timer_id) + t_comm_curr - t_comm_zero
        t_comm_zero = t_comm_curr

    end subroutine timer_comm_stamp

    subroutine timer_start(timer_id)

        integer(kind=4), intent(in) :: timer_id

        t_array(timer_id) = MPI_Wtime()

    end subroutine timer_start

    subroutine timer_end(timer_id)

        integer(kind=4), intent(in) :: timer_id

        t_array(timer_id) = MPI_Wtime() - t_array(timer_id)

    end subroutine timer_end

    function timer_elapsed(timer_id) result(t_elapsed)

        integer(kind=4), intent(in) :: timer_id
        real(kind=8)    :: t_elapsed

        t_elapsed = MPI_Wtime() - t_array(timer_id)

        return

    end function timer_elapsed

    subroutine timer_reduction

        integer(kind=4) :: ierr
        call MPI_Reduce(t_array, t_array_reduce, 40, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    end subroutine timer_reduction

    subroutine timer_output

        use mpi_topology, only  : myrank, nprocs

        integer(kind=4)     :: i
        character(len=64)   :: filename


        if(myrank.eq.0) then
            do i = 1, 40
                if(trim(t_str(i)).ne.'null') then
                    print '(a,a35,a,i3,a, f16.9)','[Timer] ', adjustl(t_str(i)),' : (',i,') : ',t_array_reduce(i) / nprocs
                endif
            enddo
        endif

        write(filename,'(a,i0.5)') 'timer_info.',myrank

        open(11, file=filename,action='write',form='formatted')
        do i = 1, 40
            if(trim(t_str(i)).ne.'null') then
                write(11, '(a,a35,a,i3,a, f16.9)') '[Timer] ', adjustl(t_str(i)),' : (',i,') : ',t_array(i)
            endif
        enddo
        close(11)

    end subroutine timer_output
    
end module timer
