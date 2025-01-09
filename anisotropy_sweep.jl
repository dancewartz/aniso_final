include("../HFModel.jl")
using .HFModel
using Base.Threads

println("Number of threads:", nthreads())

# Define fixed variables
Df = 1.0
Dh = 1.0
v0 = 1.0
λ = 20.0
s0 = 1.0 / 4
f0 = 0.1
L = 2000.0
dx = 0.1
T = 420.0 
dt = 0.002
thin_by = 4000
wave_type=:pushed

# Define parameter ranges for beta and alpha
β_values = 0.0:3.0:15.0  # Example range for β, adjust as needed
α_values = -0.1:0.01:0.3
front_speeds_matrix = zeros(length(α_values), length(β_values))
discontinuity_speeds_matrix = zeros(length(α_values), length(β_values))

for (j,β) in enumerate(β_values)
    front_speeds = zeros(length(α_values))
    discontinuity_speeds = zeros(length(α_values))

    @threads for i in eachindex(α_values)
        α = α_values[i]
        # Instantiate Parameters struct with the current set of variables
        par = Parameters(Df, Dh, v0, α, 0.0, λ, s0, f0, β, L, dx, T, dt, thin_by, wave_type)

        # Initialize frequencies and heights arrays
        frequencies = similar(par.xs)
        frequencies .= 0
        frequencies[par.xs .< -950] .= 1

        heights = similar(par.xs)
        heights .= 0

        # Run simulation with current parameters
        @time result, record_times = run_simulation!(frequencies, heights, par)

        # Calculate front locations
        front_locations = [sum(result[:, time_index, 1]) * par.dx for time_index in 1:length(record_times)]
        
        # Linear regression to find the invasion speed
        reg =  record_times .> par.T * 0.9
        X = hcat(record_times[reg], ones(length(record_times[reg])))
        invasion_speed = X \ front_locations[reg]
        
        # Store the calculated speed
        front_speeds[i] = invasion_speed[1]

        if α > 0.00
            # now to find the location of the discontinuity...
            # discontinuity_locations = [par.xs[findlast(result[:, time_index, 2] .> 0.0 * result[end, time_index, 2])] for time_index in 1:length(record_times)]
            discontinuity_locations = fill(NaN, size(record_times))
            for time_index = 1:length(record_times)
                if !isnothing(findlast(result[:, time_index, 2] .> 1.01 * result[end, time_index, 2]))
                    discontinuity_locations[time_index] = par.xs[findlast(result[:, time_index, 2] .> 1.01 * result[end, time_index, 2])]
                end
            end
            # X = hcat(record_times[reg], ones(length(record_times[reg])))
            invasion_speed = X \ discontinuity_locations[reg]
            discontinuity_speeds[i] = invasion_speed[1]
            
        end
        # if α > 0
        #     # now to find the location of the discontinuity...
        #     discontinuity_locations = Vector{Float64}(undef, length(record_times))
            
        #     for time_index in 1:length(record_times)
        #         discontinuity_index = findlast(result[:, time_index, 2] .> 0.0 * result[end, time_index, 2])
        #         if !isnothing(discontinuity_index)
        #             discontinuity_locations[time_index] = par.xs[discontinuity_index]
        #         else
        #             discontinuity_locations[time_index] = par.xs[end]  # corrected this line; the assignment to time_index was removed
        #         end
        #     end
    
        #     if any(!isnothing, discontinuity_locations)  # Ensure there's at least one valid location
        #         # Use only valid locations for regression
        #         valid_times = record_times[.!isnothing.(discontinuity_locations)]
        #         valid_locations = discontinuity_locations[.!isnothing.(discontinuity_locations)]
                
        #         # Perform linear regression if we have enough data
        #         if length(valid_times) > 1
        #             X = hcat(valid_times, ones(length(valid_times)))
        #             invasion_speed = X \ valid_locations
        #             discontinuity_speeds[i] = invasion_speed[1]
        #         end
        #     end
        # end
        
    end
    front_speeds_matrix[:, j] = front_speeds # Store speeds for this gamma
    discontinuity_speeds_matrix[:, j] = discontinuity_speeds
end

using CSV
using DataFrames

# Folder to save the files
folder_path = wave_type == :pulled ? joinpath(dirname(@__FILE__), "SimulationOutputs_pulled_new") : joinpath(dirname(@__FILE__), "SimulationOutputs_pushed_new")
isdir(folder_path) || mkdir(folder_path)

# Correcting the way to save front speeds matrix
# Automatically generate column names for the DataFrame created from the matrix
front_speeds_df = DataFrame(front_speeds_matrix, :auto)
front_speeds_file_path = joinpath(folder_path, "front_speeds_matrix.csv")
CSV.write(front_speeds_file_path, front_speeds_df, append=false)

discontinuity_speeds_df = DataFrame(discontinuity_speeds_matrix, :auto)
discontinuity_speeds_file_path = joinpath(folder_path, "discontinuity_speeds_matrix.csv")
CSV.write(discontinuity_speeds_file_path, discontinuity_speeds_df, append=false)


# Save α values
α_values_df = DataFrame(α=α_values)
α_file_path = joinpath(folder_path, "alpha_values.csv")
CSV.write(α_file_path, α_values_df, append=false)

# Save β values
β_values_df = DataFrame(β=β_values)
β_file_path = joinpath(folder_path, "beta_values.csv")
CSV.write(β_file_path, β_values_df, append=false)


# Save other parameter values
param_values = DataFrame(
    Df=[Df], 
    Dh=[Dh], 
    v0=[v0], 
    λ=[λ], 
    s0=[s0], 
    f0=[f0], 
    L=[L], 
    dx=[dx], 
    T=[T], 
    dt=[dt], 
    thin_by=[thin_by], 
    wave_type=[string(wave_type)]
)
param_file_path = joinpath(folder_path, "other_param_values.csv")
CSV.write(param_file_path, param_values, append=false)

println("Data saved in folder: $(folder_path)")
