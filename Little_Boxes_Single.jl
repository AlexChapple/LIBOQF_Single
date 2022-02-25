# A julia version of little_boxes_single.f08 to see performance differences 

using SparseArrays
using Profile 
using DelimitedFiles
using ProfileView 
using StaticArrays 

# Declare constant variables 
const start_time = 0 
const N = 100 
const end_time = 5 
const time_steps = 200000 
const num_of_simulations = 30 
const phase = pi 
const gammaL = 0.5 
const gammaR = 0.5 
const Omega = 10 * pi 
const dt = 0.000025 
const period = 80 
const tau = 0.2 


function main()

    ### ----- Main function ----- 

    # Initialise arrays 
    g_0 :: ComplexF64 = 0 
    g_0_new :: ComplexF64 = 0 
    e_0 :: ComplexF64 = 0 
    e_0_new :: ComplexF64 = 0 

    g_1 = zeros(ComplexF64, N)
    g_1_new = zeros(ComplexF64, N)
    e_1 = zeros(ComplexF64, N)
    e_1_new = zeros(ComplexF64, N)

    eta_0 = zeros(ComplexF64, time_steps)
    xi_0 = zeros(ComplexF64, time_steps)

    eta_1 = zeros(ComplexF64, N, time_steps)
    xi_1 = zeros(ComplexF64, N, time_steps)

    # Spin parameters 
    spin_up_list = zeros(Float64, time_steps)
    spin_down_list = zeros(Float64, time_steps)

    # Photon counting parameters 
    bin_width :: Int64 = 200 
    photon_list = zeros(Int64, bin_width)

    # Emission tracking parameters 
    tracking_bin_width :: Int64 = 200 
    emission_tracking_list = zeros(num_of_simulations, tracking_bin_width)

    # Construct time list 
    time_list = range(start_time, end_time, time_steps)

    # Initialise lambdaL and lambdaR 
    lambdaL :: ComplexF64 = exp(im * phase / 2) * sqrt(gammaL) * sqrt(N/tau)
    lambdaR :: ComplexF64 = exp(im * -phase / 2) * sqrt(gammaR) * sqrt(N/tau)



    ### ----- Calls function to perform one simulation at a time -----
    for sim in 1:num_of_simulations

        run_simulation(sim, g_0, g_0_new, e_0, e_0_new, g_1, g_1_new, e_1, e_1_new, lambdaL, lambdaR, bin_width, time_list, spin_down_list, spin_up_list, eta_0, eta_1, xi_0, xi_1, photon_list, emission_tracking_list)
    
    end 

    ### ----- Write to files ----- 
    # write_results_to_files(time_list, spin_up_list, spin_down_list, eta_0, eta_1, xi_0, xi_1, photon_list, emission_tracking_list, num_of_simulations)


end 

function run_simulation(sim, g_0, g_0_new, e_0, e_0_new, g_1, g_1_new, e_1, e_1_new, lambdaL, lambdaR, bin_width, time_list, spin_down_list, spin_up_list, eta_0, eta_1, xi_0, xi_1, photon_list, emission_tracking_list) 

    ### Runs one simulation and stores the outcomes 

    photon_number :: Int64 = 0 
    emission_index :: Int64 = 1 

    # Call initialise arrays 
    g_0, g_0_new, e_0, e_0_new, g_1, g_1_new, e_1, e_1_new = initialise_arrays(g_0, g_0_new, e_0, e_0_new, g_1, g_1_new, e_1, e_1_new)

    # Create random number list 
    rand_list = rand(Float64, time_steps)

    for index in 1:time_steps

        # Initialise some variables 
        g_0_new :: ComplexF64 = 0 
        e_0_new :: ComplexF64 = 0 
        g_1_new = zeros(ComplexF64, N)
        e_1_new = zeros(ComplexF64, N)

        g_0_new = (-im * Omega / 2) * e_0 

        e_0_new = ((-im * Omega / 2) * g_0) + (-im * lambdaL * g_1[1]) + (-im * lambdaR * g_1[N])

        g_1_new[1] += (-im * lambdaL * e_0)

        e_1_new[N] += (-im * lambdaR * e_0)

        g_1_new .+= (-im .* Omega ./ 2 .* e_1)
        e_1_new .+= (-im .* Omega ./ 2 .* g_1)

        # Update values 
        g_0_new = g_0 + (dt * g_0_new)
        e_0_new = e_0 + (dt * e_0_new)

        g_1_new .= g_1 .+ (dt .* g_1_new)
        e_1_new .= e_1 .+ (dt .* e_1_new)

        # Check if photon is in the Nth box 
        if mod(index, period) === 0

            # Calculate probability 
            psi_0 :: Float64 = modulo_func(g_0_new)^2 + modulo_func(e_0_new)^2

            psi_0 += sum((modulo_func_vector(g_1_new[1:N-1]).^2) .+ (modulo_func_vector(e_1_new[1:N-1]).^2)) 

            psi_1 :: Float64 = modulo_func(g_1_new[N])^2 + modulo_func(e_1_new[N])^2 

            prob = psi_1 / (psi_0 + psi_1)

            # Grab a random number from rand_list 
            rand_num :: Float64 = rand_list[index]

            ### Check if photon is in the Nth box 
            if rand_num <= prob 

                g_0 :: ComplexF64 = 1 
                e_0 :: ComplexF64 = 0 
                g_1 = zeros(ComplexF64, N)
                e_1 = zeros(ComplexF64, N)

                # Photon counting 
                photon_number += 1 

                emission_tracking_list[sim, emission_index] = time_list[index]
                emission_index += 1 

            else 

                g_0 = g_0_new 
                e_0 = e_0_new 

                g_1[1] :: ComplexF64 = 0 
                e_1[1] :: ComplexF64 = 0 

                g_1[2:N] = g_1_new[1:N-1]
                e_1[2:N] = e_1_new[1:N-1]

                g_0, e_0, g_1, e_1 = normalise_coefficients(g_0, e_0, g_1, e_1)
     

            end 


        else 

            g_0 = g_0_new 
            e_0 = e_0_new 
            g_1 = g_1_new
            e_1 = e_1_new

            g_0, e_0, g_1, e_1 = normalise_coefficients(g_0, e_0, g_1, e_1)


        end 

        # Update spin 
        update_spin(index, g_0, e_0, g_1, e_1, spin_down_list, spin_up_list)

        # Store coefficients 
        eta_0[index] += modulo_func(g_0)^2 
        xi_0[index] += modulo_func(e_0)^2 

        eta_1[:, index] .+= modulo_func_vector(g_1).^2 
        xi_1[:, index] .+= modulo_func_vector(e_1).^2 

    end # All time steps completed 

    # Store photon number into the photon list 
    if photon_number < bin_width 
        photon_list[photon_number + 1] += 1 
    end 

    # if end_time >= 15 
    #     print(sim, " simulations completed.")
    # elseif mod(sim, 10) === 0
    #     print(sim, " simulations completed.")
    # end 
    
    print(sim, " simulations completed.\n")

end 

function initialise_arrays(g_0, g_0_new, e_0, e_0_new, g_1, g_1_new, e_1, e_1_new)

    g_0 :: ComplexF64 = 1 
    g_0_new :: ComplexF64 = 0 
    e_0 :: ComplexF64 = 0 
    e_0_new :: ComplexF64 = 0 

    g_1 = zeros(ComplexF64, N)
    g_1_new = zeros(ComplexF64, N)
    e_1 = zeros(ComplexF64, N)
    e_1_new = zeros(ComplexF64, N)

    return g_0, g_0_new, e_0, e_0_new, g_1, g_1_new, e_1, e_1_new

end

function modulo_func_vector(z::Vector{ComplexF64})

    return sqrt.(real.(z).^2 + imag.(z).^2)

end 

function modulo_func(z::ComplexF64)

    return sqrt(real(z)^2 + imag(z)^2)

end 

function normalise_coefficients(g_0, e_0, g_1, e_1)

    ### Normalise coefficients 
    total :: Float64 = modulo_func(g_0)^2 + modulo_func(e_0)^2 
    total += sum((modulo_func_vector(g_1).^2) .+ (modulo_func_vector(e_1).^2))

    g_0 /= sqrt(total) 
    e_0 /= sqrt(total)
    g_1 ./= sqrt(total)
    e_1 ./= sqrt(total)

    return g_0, e_0, g_1, e_1 

end 

function update_spin(index, g_0, e_0, g_1, e_1, spin_down_list, spin_up_list)

    spin_down_prob = modulo_func(g_0)^2 
    spin_up_prob = modulo_func(e_0)^2 

    spin_down_prob += sum(modulo_func_vector(g_1).^2) 
    spin_up_prob += sum(modulo_func_vector(e_1).^2) 

    spin_down_list[index] += spin_down_prob
    spin_up_list[index] += spin_up_prob

end 

function write_results_to_files(time_list, spin_up_list, spin_down_list, eta_0, eta_1, xi_0, xi_1, photon_list, emission_tracking_list, num_of_simulations)

    eta_0 /= num_of_simulations
    xi_0 /= num_of_simulations
    eta_1 ./= num_of_simulations
    xi_1 ./= num_of_simulations

    spin_up_list ./= num_of_simulations
    spin_down_list ./= num_of_simulations

    # Write out spin data 
    open("results_julia/spin_up.dat", "w") do io
        writedlm(io, [time_list spin_up_list])
    end

    open("results_julia/spin_down.dat", "w") do io
        writedlm(io, [time_list spin_down_list])
    end

    ##### WRITE EMISSION TRACKING LATER 

    # Write photon counting data 
    open("results_julia/photon_list", "w") do io 
        writedlm(io, photon_list)
    end 

    # Write coefficient data 
    open("results_julia/eta_0.dat", "w") do io 
        writedlm(io, [time_list eta_0])
    end 

    open("results_julia/xi_0.dat", "w") do io 
        writedlm(io, [time_list xi_0])
    end 

    open("results_julia/eta_1.dat", "w") do io 
        writedlm(io, eta_1)
    end 

    open("results_julia/xi_1.dat", "w") do io 
        writedlm(io, xi_1)
    end 



end 


# Profile.clear()
# @profile main()
# @profview main()
@time main()
# Profile.print()