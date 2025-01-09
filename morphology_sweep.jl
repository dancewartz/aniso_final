include("HFModel.jl")
using .HFModel
using Plots
using Plots.PlotMeasures
using ColorSchemes
using Colors
using LaTeXStrings
using Base.Threads

println(nthreads())

# Define constants and parameters
Df = 1.0; Dh = 1.0; v0 = 0.15; γ = 0.0; λ = 20.0; s0 = 0.25; f0 = 0.1;
L = 2400.0; dx = 0.2; T = 520.0; dt = 0.01; thin_by = 100;

# Compute critical alpha and define alpha values
critical_alpha = (2 * sqrt(s0 * Df))^2 / (2 * λ)
alphas = [-0.1, 0.7 * critical_alpha, critical_alpha, 1.5 * critical_alpha, 4 * critical_alpha]
# alphas = collect(LinRange(-0.1, 5 * critical_alpha, 5))

# Define β values for sweep
betas = [0.0, 2.0, 10.0, 20.0, 25.0]  # Adjust as needed
# betas = collect(LinRange(0.0, 10.0, 5))

# Define colors and gradient
light_sky_blue = RGB(0/255, 174/255, 239/255)
yellow = RGB(249/255, 237/255, 50/255)
custom_gradient = cgrad([yellow, light_sky_blue])

# Storage for results
results_storage = Vector{Tuple{Array, Parameters}}(undef, length(alphas) * length(betas))

# Simulation function
function run_sweep(alpha, beta)
    par = Parameters(Df, Dh, v0, alpha, γ, λ, s0, f0, beta, L, dx, T, dt, thin_by, :pulled)
    frequencies = similar(par.xs); frequencies .= 0; frequencies[abs.(par.xs) .< 1] .= 1;
    heights = similar(par.xs); heights .= 0;
    result, record_times = run_simulation!(frequencies, heights, par)
    return result, par
end

# Run simulations in parallel and store results
Threads.@threads for idx in 1:length(alphas) * length(betas)
    i = div(idx - 1, length(betas)) + 1
    j = rem(idx - 1, length(betas)) + 1
    if j == 0
        j = length(betas)
        i -= 1
    end
    println(i, j)
    α = alphas[i]
    β = betas[j]
    result, par = run_sweep(α, β)  # Get result and parameters
    results_storage[idx] = (result, par)  # Store results using a flat index
end

# Plot and save results after the threaded loop
plot_directory = "anisotropy/individual_plots"
if !isdir(plot_directory)
    mkdir(plot_directory)
end

# for idx in 1:length(results_storage)
#     i = div(idx - 1, length(betas)) + 1
#     j = rem(idx - 1, length(betas)) + 1
#     if j == 0
#         j = length(betas)
#         i -= 1
#     end
#     println(i, j)
#     result, par = results_storage[idx]
#     # Plot results
#     p = plot(par.xs, result[:, 1:1:end, 2], line_z=result[:, 1:1:end, 1],
#              color=custom_gradient, label=:none, linewidth=4, xaxis=false, yaxis=false,
#              colorbar=false, grid=false)
    
#     # Save each plot individually using indices
#     println("saving", "plot_alpha_" * string(i) * "_beta_" * string(j))
#     savefig(p, joinpath(plot_directory, "plot_alpha_" * string(i) * "_beta_" * string(j) * ".pdf"))
#     savefig(p, joinpath(plot_directory, "plot_alpha_" * string(i) * "_beta_" * string(j) * ".png"))
#     savefig(p, joinpath(plot_directory, "plot_alpha_" * string(i) * "_beta_" * string(j) * ".svg"))
# end


using Printf

# Create a grid of all individual plots with titles for each alpha and beta
grid_plots = [plot() for _ in 1:length(results_storage)]
for idx in 1:length(results_storage)
    i = div(idx - 1, length(betas)) + 1
    j = rem(idx - 1, length(betas)) + 1
    if j == 0
        j = length(betas)
        i -= 1
    end
    result, par = results_storage[idx]
    
    # Format alpha and beta to two decimal places for the title
    plot_title = @sprintf("α = %.3f, β = %.3f", alphas[i], betas[j])

    # Generate the individual plot with specified title and adjustments for size
    grid_plots[idx] = plot(
        par.xs, result[:, 1:1:end, 2], line_z=result[:, 1:1:end, 1],
        color=custom_gradient, label=:none, linewidth=1.5,
        xaxis=false, yaxis=false, colorbar=false, grid=false,
        title=plot_title, titlefontsize=10  # Set title and font size
    )

    println(i, j)
end

# Arrange and save the grid of all plots with modified layout
final_grid_plot = plot(
    grid_plots..., layout=(length(alphas), length(betas)),
    size=(1200, 800),  # Adjust the grid size to make individual plots larger
    padding=5mm        # Optional: add padding between plots
)
savefig(final_grid_plot, joinpath(plot_directory, "morphology_sweep_grid.png"))
# savefig(final_grid_plot, joinpath(plot_directory, "morphology_sweep_grid.pdf"))
# savefig(final_grid_plot, joinpath(plot_directory, "morphology_sweep_grid.svg"))