module HFModel
    export Parameters, run_simulation!

    struct Parameters
        Df        :: Float64
        Dh        :: Float64
        v0        :: Float64
        α         :: Float64
        γ         :: Float64
        λ         :: Float64
        s0        :: Float64
        f0        :: Float64
        β         :: Float64
        L         :: Float64
        dx        :: Float64
        T         :: Float64
        dt        :: Float64
        thin_by   :: Int64
        xs        :: Vector{Float64}
        wave_type :: Symbol
        λ1        :: Float64
        β1        :: Float64
    
        # Inner constructor
        function Parameters(
            Df::Float64, Dh::Float64, v0::Float64, α::Float64, γ::Float64, λ::Float64, s0::Float64, f0::Float64,
            β::Float64, L::Float64, dx::Float64, T::Float64, dt::Float64,
            thin_by::Int64=1, wave_type::Symbol=:pushed, λ1::Float64=0.0, β1::Float64=0.0
        )
            xs = collect(-L/2:dx:L/2)  # Compute xs based on L and dx
            new(Df, Dh, v0, α, γ, λ, s0, f0, β, L, dx, T, dt, thin_by, xs, wave_type, λ1, β1)
        end
    end
    
    

    function run_simulation!(frequencies, heights, par::Parameters)
        @assert length(frequencies) == length(par.xs) && length(heights) == length(par.xs) "Array sizes mismatch"

        temp_data_1 = similar(frequencies)
        temp_data_2 = similar(heights)

        record_times = 0:par.dt*par.thin_by:par.T # store the times at which we will record the result
        result = fill(NaN,length(par.xs),length(record_times),2) # we will store the result here

        result[:,1,1] .= frequencies; result[:,1,2] .= heights;

        for time_idx = 2:length(record_times)
            update_fields!(frequencies, heights, temp_data_1, temp_data_2, par)
            result[:,time_idx,1] .= frequencies; result[:,time_idx,2] .= heights;
        end

        return result, record_times
    end

    function update_fields!(frequencies, heights, temp_data_1, temp_data_2, par::Parameters)
        for i = 1:par.thin_by
            apply_diffusion!(frequencies, heights, temp_data_1, temp_data_2, par)
            par.wave_type == :pushed ? apply_growth_and_advection!(frequencies, heights, temp_data_1, temp_data_2, par) : apply_growth_and_advection_pulled!(frequencies, heights, temp_data_1, temp_data_2, par)
        end
    end

    function apply_diffusion!(frequencies, heights, temp_data_1, temp_data_2, par::Parameters)
        for i = 2:length(frequencies)-1
            temp_data_1[i] = par.Df * par.dt / par.dx^2 * (frequencies[i-1] + frequencies[i+1] - 2 * frequencies[i])
            temp_data_2[i] = par.Dh * par.dt / par.dx^2 * (heights[i-1] + heights[i+1] - 2 * heights[i])
        end

        # boundary conditions are no flux 
        temp_data_1[1] = par.Df * par.dt / par.dx^2 * (frequencies[1] + frequencies[2] - 2 * frequencies[1]) # left end
        temp_data_1[end] = par.Df * par.dt / par.dx^2 * (frequencies[end-1] + frequencies[end] - 2 * frequencies[end]) # right end

        temp_data_2[1] = par.Dh * par.dt / par.dx^2 * (heights[1] + heights[2] - 2 * heights[1]) # left end 
        temp_data_2[end] = par.Dh * par.dt / par.dx^2 * (heights[end-1] + heights[end] - 2 * heights[end]) # right end

        # update fields in place
        frequencies .+= temp_data_1
        heights .+= temp_data_2
    end

    function apply_growth_and_advection!(frequencies, heights, temp_data_1, temp_data_2, par::Parameters)
        for i = 2:length(frequencies)-1
            temp_data_1[i] = par.dt * (par.s0 * (frequencies[i] - par.f0) * frequencies[i] * (1 - frequencies[i]) + 
                             (par.β + par.β1 * frequencies[i]) * (frequencies[i+1] - frequencies[i-1]) / (2 * par.dx) * (heights[i+1] - heights[i-1]) / (2 * par.dx))
            temp_data_2[i] = par.dt * (par.v0 + par.α * frequencies[i] + par.γ * frequencies[i] * (1 - frequencies[i]) +
                             (par.λ + par.λ1 * frequencies[i])/2 * (heights[i+1] - heights[i-1])^2 / (2 * par.dx)^2 )
        end
        
        temp_data_1[1] = par.dt * (par.s0 * (frequencies[1] - par.f0) * frequencies[1] * (1 - frequencies[1]) + 
                         (par.β + par.β1 * frequencies[1]) * (frequencies[2] - frequencies[1]) / (2 * par.dx) * (heights[2] - heights[1]) / (2 * par.dx))

        temp_data_1[end] = par.dt * (par.s0 * (frequencies[end] - par.f0) * frequencies[end] * (1 - frequencies[end]) + 
                           (par.β + par.β1 * frequencies[end]) * (frequencies[end] - frequencies[end-1]) / (2 * par.dx) * (heights[end] - heights[end-1]) / (2 * par.dx))

        temp_data_2[1] = par.dt * (par.v0 + par.α * frequencies[1] + par.γ * frequencies[1] * (1 - frequencies[1]) +
                         (par.λ + par.λ1 * frequencies[1])/2 * (heights[2] - heights[1])^2 / (2 * par.dx)^2 )

        temp_data_2[end] = par.dt * (par.v0 + par.α * frequencies[end] + par.γ * frequencies[end] * (1 - frequencies[end]) +
                         (par.λ + par.λ1 * frequencies[end])/2 * (heights[end] - heights[end-1])^2 / (2 * par.dx)^2 )

        frequencies .+= temp_data_1
        heights .+= temp_data_2
    end

    function apply_growth_and_advection_pulled!(frequencies, heights, temp_data_1, temp_data_2, par::Parameters)
        for i = 2:length(frequencies)-1
            temp_data_1[i] = par.dt * (par.s0 * frequencies[i] * (1 - frequencies[i]) + 
                             (par.β + par.β1 * frequencies[i]) * (frequencies[i+1] - frequencies[i-1]) / (2 * par.dx) * (heights[i+1] - heights[i-1]) / (2 * par.dx))
            temp_data_2[i] = par.dt * (par.v0 + par.α * frequencies[i] + par.γ * frequencies[i] * (1 - frequencies[i]) +
                             (par.λ + par.λ1 * frequencies[i])/2 * (heights[i+1] - heights[i-1])^2 / (2 * par.dx)^2 )
        end
        
        temp_data_1[1] = par.dt * (par.s0 * frequencies[1] * (1 - frequencies[1]) + 
                         (par.β + par.β1 * frequencies[1]) * (frequencies[2] - frequencies[1]) / (2 * par.dx) * (heights[2] - heights[1]) / (2 * par.dx))

        temp_data_1[end] = par.dt * (par.s0 * frequencies[end] * (1 - frequencies[end]) + 
                           (par.β + par.β1 * frequencies[end]) * (frequencies[end] - frequencies[end-1]) / (2 * par.dx) * (heights[end] - heights[end-1]) / (2 * par.dx))

        temp_data_2[1] = par.dt * (par.v0 + par.α * frequencies[1] + par.γ * frequencies[1] * (1 - frequencies[1]) +
                         (par.λ + par.λ1 * frequencies[1])/2 * (heights[2] - heights[1])^2 / (2 * par.dx)^2 )

        temp_data_2[end] = par.dt * (par.v0 + par.α * frequencies[end] + par.γ * frequencies[end] * (1 - frequencies[end]) +
                         (par.λ + par.λ1 * frequencies[end])/2 * (heights[end] - heights[end-1])^2 / (2 * par.dx)^2 )

        frequencies .+= temp_data_1
        heights .+= temp_data_2
    end
end
