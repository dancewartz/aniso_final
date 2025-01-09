using CSV
using DataFrames
using Plots
using Colors
using ColorSchemes
using LaTeXStrings
pgfplotsx()  # Use PGFPlotsX backend

# Define the folder path where your files are located
folder_path = joinpath(dirname(@__FILE__), "SimulationOutputs_pulled_new")

# Load the data
front_speeds_df = CSV.read(joinpath(folder_path, "front_speeds_matrix.csv"), DataFrame)
discontinuity_speeds_df = CSV.read(joinpath(folder_path, "discontinuity_speeds_matrix.csv"), DataFrame)
α_values_df = CSV.read(joinpath(folder_path, "alpha_values.csv"), DataFrame)
β_values_df = CSV.read(joinpath(folder_path, "beta_values.csv"), DataFrame)
param_values_df = CSV.read(joinpath(folder_path, "other_param_values.csv"), DataFrame)

# Convert the front speeds DataFrame back to a matrix for easier indexing
front_speeds_matrix = Matrix(front_speeds_df)
discontinuity_speeds_matrix = Matrix(discontinuity_speeds_df)

# Extract alpha and beta values
α_values = α_values_df.α
β_values = β_values_df.β

# Extracting parameters from the DataFrame
v0 = param_values_df.v0[1]
λ = param_values_df.λ[1]
s0 = param_values_df.s0[1]
Df = param_values_df.Df[1]
f0 = param_values_df.f0[1]
u0 = sqrt(s0 * Df ) * 2

# Kap constant
kap = 0.0
gam = 1.0
# u0 = sqrt(s0 * Df / 2) * (1 - 2 * f0)
colors = get(ColorSchemes.viridis, range(0, stop=1, length=length(β_values)+2))

# Initialize the plot
p = plot(
    xlabel=L"expansion speed difference, $\alpha$",
    ylabel=L"invasion speed, $u/u_0$",
    legend=:topleft,
    xguidefontsize=20,
    yguidefontsize=20,
    legendfontsize=8,
    tickfontsize=15
)

# Define marker types
markers = [:circle, :diamond, :square, :utriangle, :cross, :dtriangle]

# Define RGB colors for regions
color_negative = RGB(0.4, 0.2, 0.6)  # Red for alpha < 0
color_positive_above = RGB(0.9, 0.6, 0.3)  # Green for alpha > 0 above sqrt curve
color_positive_below = RGB(0.2, 0.6, 0.8)  # Blue for alpha > 0 below sqrt curve

# Fill regions
α_values_fill = LinRange(-maximum(abs.(α_values)), maximum(abs.(α_values)), 1000)
α_values_fill = LinRange(α_values[1],α_values[end], 3000 )
velocities_sqrt = sqrt.(2 * λ * abs.(α_values_fill)) / u0

# Fill alpha < 0
plot!(p, α_values_fill[α_values_fill .< 0], zeros(size(velocities_sqrt[α_values_fill .< 0])), ribbon=1.1*maximum(velocities_sqrt), color=color_negative, fillalpha=0.5, label=:none)

# Fill alpha > 0 above sqrt curve
plot!(p, α_values_fill[α_values_fill .> 0], velocities_sqrt[α_values_fill .> 0], fillrange=1.1*maximum(velocities_sqrt) .- 0*velocities_sqrt[α_values_fill .> 0], color=color_positive_above, fillalpha=0.5, label=:none)

# Fill alpha > 0 below sqrt curve
plot!(p, α_values_fill[α_values_fill .> 0], zeros(size(α_values_fill[α_values_fill .> 0])), ribbon=velocities_sqrt[α_values_fill .> 0], color=color_positive_below, fillalpha=0.5, label=:none)

# Plot each curve
for i in length(β_values):-1:1
    β = β_values[i]

    if β < 0
        continue
    end

    kap_eff = kap * β / λ

    # Theoretical expression calculation
    reg = α_values .<= u0^2 / 2 / λ / (1 - kap_eff)^2
    u_theo = (u0 * (-1 + kap_eff) .+ kap_eff .* sqrt.(u0^2 .+ 2 .* α_values[reg] .* (-1 + 2 * kap_eff) * λ)) ./ (-1 + 2 * kap_eff)

    # Scatter plots (black, different markers)
    marker = markers[mod1(i, length(markers))]
    scatter!(p, α_values[1:2:end], front_speeds_matrix[1:2:end, i] / u0, label="β/λ = $(β/λ)", linewidth=2, markercolor=colors[i], marker=marker, markeralpha=0.6)

    # Dashed theoretical line (black)
    plot!(p, α_values[reg], u_theo / u0, color=colors[i], lw=2, linestyle=:dash, label="β/λ = $(β/λ)")

    # Analytical line (dash-dotted, black)
    reg = α_values .> 0
    analytical_line = (u0 .+ gam .* β ./ param_values_df.λ .* sqrt.(2 * param_values_df.λ .* α_values[reg])) ./ (1 - β / param_values_df.λ[1])
    analytical_line[analytical_line .> sqrt.(2 * param_values_df.λ .* α_values[reg])] .= NaN
    plot!(p, α_values[reg], analytical_line / u0, color=colors[i], lw=2, linestyle=:dashdot, label=:none)
end

# Circular arc curve (black, dashed)
α_values_theoretical = LinRange(0, α_values[end], 1000)
plot!(p, α_values_theoretical, sqrt.(2 * λ * α_values_theoretical) / u0, lw=5, linestyle=:solid, color=RGB(0.6, 0.9, 0.6), label="Circular Arc", alpha=0.9, fillalpha=0.8)
plot!(p, α_values_theoretical, sqrt.(2 * λ * α_values_theoretical) / u0, lw=2, linestyle=:dash, color=:black, label=false)

# Adjust the y-limits if necessary
ylims!(p, 0, 1.1 * maximum(sqrt.(2 * λ * α_values_theoretical)) / u0)
xlims!(p, α_values[1], α_values[end])
# Save the plot
plot_file_path = joinpath(folder_path, "pulled_speed_vs_alpha_theoretical.pdf")
savefig(p, plot_file_path)

println("Plot saved to: $(plot_file_path)")

#########################################################
# Initialize the plot
p =plot(
    xlabel=L"expansion speed difference, $\alpha$",
    ylabel=L"slope discontinuity speed, $c$",
    legend=:topleft,
    xguidefontsize=20,
    yguidefontsize=20,
    legendfontsize=8,
    tickfontsize=15
)

# Plot each curve
for i in length(β_values):-1:1
    β = β_values[i]

    if β < 0
        continue
    end

    color = colors[i]
    marker = markers[mod1(i, length(markers))]
        
    reg = α_values .> 0.0
    analytical_line = (2 * sqrt.(param_values_df.s0 .* param_values_df.Df) .+ β ./ param_values_df.λ .* sqrt.(2 * param_values_df.λ .* α_values[reg])) ./ (1 - β / param_values_df.λ[1])
    # plot!(p, α_values[reg], 1/2 * (analytical_line .+ sqrt.(2 * λ * α_values[reg])), color=color, lw=2, linestyle=:dot, label=:none)
    omg = 1/2 * (analytical_line .+ sqrt.(2 * λ * α_values[reg]))
    omg[omg .> sqrt.(2 * λ * α_values[reg])] .= NaN
    plot!(p, α_values[reg],α_values[reg], color=color, lw=2, linestyle=:dot, label=:none, alpha=0)
    plot!(p, α_values[reg],analytical_line, color=color, lw=2, linestyle=:dot, label=:none, alpha=0)
    plot!(p, α_values[reg],α_values[reg] .+ analytical_line, color=color, lw=2, linestyle=:dot, label=:none, alpha=0)
    plot!(p, α_values[reg], omg, color=color, lw=2, linestyle=:dot, label="β/λ = $(β/λ)")
    println(size(α_values[reg]))
    println(size(analytical_line))
    println(any(isnan.(analytical_line)))
    println(any(isnan.(α_values[reg])))
    
    
    α_values_positive = α_values[α_values .> 0.0]
    disc_speeds_positive = discontinuity_speeds_matrix[α_values .> 0.0, i]
    reg = analytical_line .< sqrt.(2 * λ * α_values_positive)
    # disc_speeds_positive[disc_speeds_positive .> sqrt.(2 * λ * α_values_positive)] .= NaN

    # Plot simulation results
    scatter!(p, α_values_positive[1:2:end], disc_speeds_positive[1:2:end], label="β/λ = $(β/λ)", linewidth=2, markercolor=color, marker=marker, markeralpha=0.6)
end

plot!(p, α_values_theoretical, sqrt.(2 * λ * α_values_theoretical), lw = 2, label="Circular Arc", linestyle=:dash, color=:black)

ylims!(p, 0, 1.1 * maximum(sqrt.(2 * λ * α_values_theoretical)))
# Save the plot
plot_file_path = joinpath(folder_path, "pulled_discontinuity_vs_alpha_theoretical.pdf")
savefig(p, plot_file_path)

println("Plot saved to: $(plot_file_path)")
###########################################################
# Initialize the plot
p = plot(
    xlabel="Beta (β)",
    ylabel="Measured Speed",
    title="Measured Speed vs Beta for different Alpha values",
    legend=:best
)

# Define colors for the plot, ensuring there are enough colors for each α
colors = distinguishable_colors(length(α_values))

# Plot each curve
for i in 1:length(α_values)
    α = α_values[i]
    speeds_for_alpha = [front_speeds_matrix[i, j] for j in 1:length(β_values)]  # Extract speeds for the current α across all β

    # Plot simulation results for this α
    scatter!(p, β_values, speeds_for_alpha, label="α = $α", linewidth=2, markercolor=colors[i], marker=:circle, markeralpha=0.5)
    
end

# Adjust the y-limits if necessary
ylims!(p, min(0, minimum(front_speeds_matrix)), 2 * maximum(front_speeds_matrix))

# Save the plot
plot_file_path = joinpath(folder_path, "velocity_vs_beta_different_alpha.pdf")
savefig(p, plot_file_path)

println("Plot saved to: $(plot_file_path)")