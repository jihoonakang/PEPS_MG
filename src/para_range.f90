subroutine para_range(n1, n2, nprocs, myrank, ista, iend)

        implicit none

        integer(kind=4), intent(in)     :: n1, n2, nprocs, myrank
        integer(kind=4), intent(out)    :: ista, iend
        integer(kind=4) :: iwork1, iwork2

        iwork1 = (n2 - n1 + 1) / nprocs
        iwork2 = mod(n2 - n1 + 1, nprocs)
        ista = myrank * iwork1 + n1 + min(myrank, iwork2)
        iend = ista + iwork1 - 1
        if (iwork2 > myrank) iend = iend + 1

end subroutine para_range


