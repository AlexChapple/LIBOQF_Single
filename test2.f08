program main 


    implicit none 

    real, dimension(100000000) :: a, b, c 
    integer :: beginning, end, rate, i 

    call system_clock(beginning, rate)

    call random_number(a)
    call random_number(b)

    c = a * b 
    ! do i = 1, size(a)
    !     c(i) = a(i) * b(i)
    ! end do 

    call system_clock(end)

    print *, "Execution time: ", real(end - beginning) / real(rate), " seconds."

end program main 