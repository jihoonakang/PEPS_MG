!======================================================================================================================
!> @file        para_range.f90
!> @brief       The para_range function assigns the computing range to each MPI process.
!> @author      
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>
!> @date        May 2021
!> @version     1.0
!> @par         Copyright
!>              Copyright (c) 2021 Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE in )
!======================================================================================================================

!> @brief       Compute the indices of the assigned range for each MPI process .
!> @param       n1      First index of total range
!> @param       n2      Last index of total range
!> @param       nprocs  Number of MPI process
!> @param       myrank  Rank ID of my MPI process
!> @param       ista    First index of assigned range for myrank
!> @param       iend    Last index of assigned range for myrank

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


