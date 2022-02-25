! Simulation of the quantum system using quantum trajectories and the SDW model. 
! This file is for the case where only one photon is in the system 
! 
! Author: Alex Chapple 
!
!
! NOTES: not working, can't get the decay to show 

program main 

    implicit none
    
    ! Declare parameters  
    real (kind=8) :: start_time = 0.0d0  
    integer (kind=8), parameter :: N = 100d0 
    integer (kind=8), parameter :: end_time = 5d0 
    integer (kind=8), parameter :: time_steps = 200000d0  
    integer (kind=8), parameter :: num_of_simulations = 30d0
    real (kind=8), parameter :: pi = 3.14159265358979323846d0
    real (kind=8), parameter :: phase = pi !0.0d0 
    real (kind=8), parameter :: gammaL = 0.5d0 
    real (kind=8), parameter :: gammaR = 0.5d0
    real (kind=8), parameter :: Omega = 10.0d0 * pi 
    real (kind=8), parameter :: dt = 0.000025d0 
    integer (kind=8), parameter :: period = 80d0 
    real (kind=8), parameter :: tau = 0.2d0 
    complex (kind=8), parameter :: i = cmplx(0.0d0,1.0d0,kind=8) 

    ! Declare flexible variables 
    real (kind=8) :: total
    integer (kind=8) :: sim, index, q, j, k, beginning, end, rate, log_line, index1, index2
    real (kind=8), dimension(time_steps) :: time_list, rand_list
    complex (kind=8) :: lambdaL, lambdaR
    real (kind=8) :: psi_0, psi_1, prob, rand_num, spin_up_prob, spin_down_prob, spin_total

    ! Declare state expansion coefficients 
    complex (kind=8) :: g_0, g_0_new, e_0, e_0_new
    complex (kind=8), dimension(N) :: g_1, g_1_new, e_1, e_1_new 

    ! Coefficient tracking 
    real (kind=8), dimension(time_steps) :: eta_0, xi_0
    real (kind=8), dimension(N, time_steps) :: eta_1, xi_1 

    ! Spin parameters 
    real (kind=8), dimension(time_steps) :: spin_up_list, spin_down_list

    ! Photon counting parameters 
    integer (kind=8), parameter :: bin_width = 200d0 
    integer (kind=8) :: photon_number
    integer (kind=8), dimension(bin_width) :: photon_list

    ! Emission tracking variables 
    integer (kind=8) :: emission_index
    integer (kind=8), parameter :: tracking_bin_width = 200
    real (kind=8), dimension(num_of_simulations, tracking_bin_width) :: emission_tracking_list  

    ! -------------------------------------------------------------------
    !
    ! START PROGRAM
    ! 
    ! -------------------------------------------------------------------

    ! Print and store parameters 
    call print_info(tau, period, dt, num_of_simulations, Omega, end_time, &
                            time_steps, phase, gammaL, gammaR, N)

    ! Start the program duration stopwatch  
    call system_clock(beginning, rate)

    ! Initialise spin up and down lists 
    spin_up_list = 0.0d0; spin_down_list = 0.0d0 

    ! Initialise photon list 
    photon_list = 0d0 

    ! Construct time list 
    call linspace(start=start_time, end=end_time, list=time_list)

    ! Initialise lambdaL and lambdaR 
    lambdaL = exp(cmplx(0.0d0, phase / 2, kind=8)) * sqrt(gammaL) * sqrt(N/tau)
    lambdaR = exp(cmplx(0.0d0, -phase / 2, kind=8)) * sqrt(gammaR) * sqrt(N/tau)

    
    ! Iterate over each simulation 
    do sim = 1, num_of_simulations

        photon_number = 0d0 ! Initialises photon number 
    
        emission_index = 1d0 ! Initialises emission tracking line 

        ! Initialise arrays 
        call initialise_arrays(N, g_0, g_0_new, e_0, e_0_new, g_1, g_1_new,&
                                 e_1 ,e_1_new)

        ! Construct random number list 
        call random_number(rand_list)

        ! Iterate over each time step 
        do index = 1, size(time_list)

            ! Initialise some variables 
            psi_0 = 0.0d0; psi_1 = 0.0d0; prob = 0.0d0; rand_num = 0.0d0; spin_up_prob = 0.0d0; 
            spin_down_prob = 0.0d0; spin_total = 0.0d0; total = 0.0d0

            ! Update the expansion coefficients 
            g_0_new = 0.0d0; e_0_new = 0.0d0; g_1_new = 0.0d0; e_1_new = 0.0d0

            g_0_new = g_0_new + ((-i) * (Omega/2) * e_0)

            e_0_new = e_0_new + ((-i) * (Omega/2) * g_0) + ((-i) * (lambdaL * g_1(1))) + ((-i) * (lambdaR * g_1(N)))

            g_1_new(1) = g_1_new(1) + ((-i) * lambdaL * e_0)

            g_1_new(N) = g_1_new(N) + ((-i) * lambdaR * e_0)

            do j = 1, N 
                g_1_new(j) = g_1_new(j) + ((-i) * (Omega/2) * e_1(j))
                e_1_new(j) = e_1_new(j) + ((-i) * (Omega/2) * g_1(j))
            end do 

            g_0_new = g_0 + (dt*g_0_new)
            e_0_new = e_0 + (dt*e_0_new)
            g_1_new = g_1 + (dt*g_1_new)
            e_1_new = e_1 + (dt*e_1_new)



            ! Check if photon is in the Nth box  
            if (mod(index, period) == 0) then 

                ! Calculate probabilities 
                psi_0 = (modulo_func(g_0_new)**2) + (modulo_func(e_0_new)**2)

                do j = 1,(N-1) 
                    psi_0 = psi_0 + (modulo_func(g_1_new(j))**2) + (modulo_func(e_1_new(j))**2)
                end do 

                psi_1 = (modulo_func(g_1_new(N))**2) + (modulo_func(e_1_new(N))**2)

                prob = psi_1 / (psi_1 + psi_0)
                
                ! Grab a random number from rand_list 
                rand_num = rand_list(index)

                if (rand_num <= prob) then ! photon is in the Nth box 

                    g_0 = 1.00; e_0 = 0.0d0 

                    g_1 = 0.0d0; e_1 = 0.0d0 

                    ! Store the information in the local spin_up/down_list 
                    spin_up_list(index) = spin_up_list(index) + spin_up_prob 
                    spin_down_list(index) = spin_down_list(index) + spin_down_prob

                    ! Photon counting done here 
                    photon_number = photon_number + 1

                    ! Photon emission tracking done here
                    emission_tracking_list(sim, emission_index) = time_list(index) ! Saves the emission time to the end of the list 
                    emission_index = emission_index + 1 ! Increases the emission index for the next emission 
                

                else ! photon is not in the Nth box 

                    g_0 = g_0_new; e_0 = e_0_new 

                    g_1(1) = 0.0d0; e_1(1) = 0.0d0 

                    do j = 2, N
                        g_1(j) = g_1_new(j-1)
                        e_1(j) = e_1_new(j-1)
                    end do

                    !!! Normalise coefficients 
                    total = modulo_func(g_0)**2 + modulo_func(e_0)**2 

                    do j = 1,N 
                        total = total + modulo_func(g_1(j))**2 + modulo_func(e_1(j))**2 
                    end do 
                    
                    g_0 = g_0 / sqrt(total); e_0 = e_0 / sqrt(total) 
                    g_1 = g_1 / sqrt(total); e_1 = e_1 / sqrt(total) 

                    ! Calculate spin up and down 
                    spin_down_prob = modulo_func(g_0)**2 
                    spin_up_prob = modulo_func(e_0)**2 

                    do j = 1,N 
                        spin_down_prob = spin_down_prob + (modulo_func(g_1(j))**2)
                        spin_up_prob = spin_up_prob + (modulo_func(e_1(j))**2)
                    end do 

                    ! Store the information in the local spin_up/down_list 
                    spin_up_list(index) = spin_up_list(index) + spin_up_prob 
                    spin_down_list(index) = spin_down_list(index) + spin_down_prob


                end if 

            else ! Not a multiple of a period

                g_0 = g_0_new 
                e_0 = e_0_new 
                g_1 = g_1_new
                e_1 = e_1_new

                !!! Normalise coefficients 
                total = modulo_func(g_0)**2 + modulo_func(e_0)**2 

                do j = 1,N 
                    total = total + modulo_func(g_1(j))**2 + modulo_func(e_1(j))**2 
                end do 
                
                g_0 = g_0 / sqrt(total); e_0 = e_0 / sqrt(total) 
                g_1 = g_1 / sqrt(total); e_1 = e_1 / sqrt(total) 

                ! Calculate spin up and down 
                spin_down_prob = modulo_func(g_0)**2 
                spin_up_prob = modulo_func(e_0)**2 

                do j = 1,N 
                    spin_down_prob = spin_down_prob + (modulo_func(g_1(j))**2)
                    spin_up_prob = spin_up_prob + (modulo_func(e_1(j))**2)
                end do 

                ! Store the information in the local spin_up/down_list 
                spin_up_list(index) = spin_up_list(index) + spin_up_prob 
                spin_down_list(index) = spin_down_list(index) + spin_down_prob

            end if 

            ! Store coefficients 
            eta_0(index) = eta_0(index) + modulo_func(g_0)**2 
            xi_0(index) = xi_0(index) + modulo_func(e_0)**2 

            do j = 1, N 
                eta_1(j, index) = eta_1(j, index) + modulo_func(g_1(j))**2 
                xi_1(j, index) = xi_1(j, index) + modulo_func(e_1(j))**2 
            end do 

        end do ! All time steps completed 

        ! Store photon number into the photon list 
        if (photon_number < bin_width) then 

            photon_list(photon_number + 1) = photon_list(photon_number + 1) + 1

        end if 

        ! Prints simulation completions to console 
        ! if (end_time >= 15) then 
        !     print *, sim ,' simulations completed.'
        ! else
        !     if (mod(sim, 10) == 0) then 
        !         print *, sim ,' simulations completed.'
        !     end if     
        ! end if 

        print *, sim ,' simulations completed.'

    end do 

    call system_clock(end)

    print *, "All simulations completed. Execution time: ", real(end - beginning) / real(rate), " seconds."



    ! -------------------------------------------------------------------------------------------
    ! 
    !   Functions and Subroutines 
    !
    !-------------------------------------------------------------------------------------------
    

    contains 

    subroutine linspace(start, end, list) 

        real (kind=8), intent(in) :: start
        integer (kind=8), intent(in) :: end 
        real (kind=8), intent(out) :: list(:)
        real (kind=8) :: range
        integer (kind=8) :: p, q

        p = size(list)
        range = end - start
        
        do q = 1,p
            list(q) = start + (range * (q - 1) / (p - 1))
        end do

    end subroutine

    function modulo_func(z) result(c)

        ! Takes in a complex number z and returns its modulus

        implicit none

        ! Declare var types
        real (kind=8) :: a, b, c
        complex (kind=8), intent(in) :: z

        a = real(z)
        b = aimag(z)

        c = sqrt(a**2 + b**2)

    end function 

    subroutine initialise_arrays(N, g_0, g_0_new, e_0 , e_0_new, g_1, &
                                 g_1_new, e_1 ,e_1_new)

        ! This stores the initial condition of the simulations

        ! Declare types 
        integer (kind=8) :: N
        complex (kind=8) :: g_0, g_0_new, e_0 ,e_0_new
        complex (kind=8), dimension(N) :: g_1, g_1_new, e_1 ,e_1_new

        g_0 = 1.0d0; g_0_new = 0.0d0
        e_0 = 0.0d0; e_0_new = 0.0d0

        g_1 = 0.0d0; g_1_new = 0.0d0
        e_1 = 0.0d0; e_1_new = 0.0d0

    end subroutine

    subroutine print_info(tau, period, dt, num_of_simulations, Omega, end_time, &
                            time_steps, phase, gammaL, gammaR, N)

        integer (kind=8) :: N, end_time, time_steps, num_of_simulations, period
        real (kind=8) :: phase, gammaL, gammaR, dt, tau, Omega

        print *, "Number of Boxes: ", N 
        print *, "simulation duration: ", end_time 
        print *, "time steps: ", time_steps
        print *, "number of simulations: ", num_of_simulations
        print *, "phase: ", phase 
        print *, "gammaL: ", gammaL 
        print *, "gammaR: ", gammaR 
        print *, "Omega: ", real(Omega)
        print *, "dt: ", dt 
        print *, "tau: ", tau 
        print *, "period: ", period 
        print *, "Delta t: ", tau / N 

        ! Also will write this to an input.txt file 
        open(5, file="results/input.txt", status="replace")

        write (5,*) "Number of Boxes: ", N 
        write (5,*) "simulation duration: ", end_time 
        write (5,*) "time steps: ", time_steps
        write (5,*) "number of simulations: ", num_of_simulations
        write (5,*) "phase: ", phase 
        write (5,*) "gammaL: ", gammaL 
        write (5,*) "gammaR: ", gammaR 
        write (5,*) "Omega: ", real(Omega)
        write (5,*) "dt: ", dt 
        write (5,*) "tau: ", tau 
        write (5,*) "period: ", period 
        write (5,*) "Delta t: ", tau / N 

        close(5)

    end subroutine


end program main 