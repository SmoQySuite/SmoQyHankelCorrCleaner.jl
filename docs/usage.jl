# # Usage
#
# Here we demonstrate how to use the [SmoQyHankelCorrCleaner.jl](https://github.com/SmoQySuite/SmoQyHankelCorrCleaner.jl) package.

using SmoQyHankelCorrCleaner
using Statistics
using CairoMakie
CairoMakie.activate!(type = "svg")

# As an initial example, we consider the atomic Hubbard model
# ```math
# H = U (\hat{n}_\uparrow - \tfrac{1}{2})(\hat{n}_\downarrow - \tfrac{1}{2}) - \mu (\hat{n}_\uparrow + \hat{n}_\downarrow)
# ```
# with the correspoding imaginary-time Green's function given by
# ```math
# G_\sigma(\tau) = \langle \hat{c}_\sigma^{\dagger}(\tau) \hat{c}_\sigma^{\phantom\dagger}(0) \rangle
#     = \frac{e^{\tau U/2 + \tau \mu - \beta U/4} + e^{-\tau U/2 + \tau \mu + \beta U/4 + \beta \mu}}{e^{-\beta U/4} + 2 e^{\beta (U/4+\mu)} + e^{-\beta(U/4-2\mu)}},
# ```
# where ``\tau \in [0,\beta)`` and ``\beta = 1/T`` is the inverse temperature.

## Evaluate G(τ) for the atomic Hubbard model
function atomic_hubbard_greens(τ, β, U, μ)

    Z   = exp(-β*U/4) + 2*exp(β*(U/4+μ)) + exp(-β*(U/4-2*μ))
    ZGτ = exp(τ*U/2 + τ*μ - β*U/4) + exp(-τ*U/2 + τ*μ + β*U/4 + β*μ)
    Gτ  = ZGτ / Z
    return Gτ
end

# Now we assign parameter values in our model, including the inverse temperature.

## Hubbard repulsion
U = 5.0

## Chemical potential
μ = 2.50

## Inverse temperature
β = 5.0;

# Next we evaulate the Green's funciton on a regular imaginary-time $\tau$ grid.

## Discretization interval
Δτ = 0.05

## Imaginary-time values
τ = collect(range(start = 0.0, stop = β, step = Δτ))

## Evaluate G(τ)
Gτ_exact = atomic_hubbard_greens.(τ, β, U, μ)

# Now let us define some synthetic noisy data.

## Add normal random noise to Green's function
Gτ_noisy = Gτ_exact + 0.007 * randn(length(τ))

## Reduce noise for G(τ = 0) point
Gτ_noisy[1] = Gτ_exact[1] +  0.002 * randn()

## Ensure that noisy data still satisty the sum rule G(τ=0) + G(τ=β) = 1
Gτ_noisy[end] = 1 - Gτ_noisy[1]

Gτ_noisy

# Now let us denoise the noisy Green's function using the hankel projection method.

## Denoise imaginary-time Green's function
Gτ_clean, iter, err =  hankel_correlation_cleaner(
    Gτ_noisy,
    maxiter = 1000,
    tol = 1e-3,
    positive_curvature = true,
    fixed_endpoints = true,
    symmetric = false,
    verbose = false
)

(iter, err)

# Now let us plot the result to see how this all looks.

fig = Figure(
    size = (700, 700),
    fonts = (; regular= "CMU Serif"),
    figure_padding = 10
)

ax = Axis(fig[1, 1],
    aspect = 1,
    xlabel = L"\tau",
    ylabel = L"G_\sigma(\tau)",
    xlabelsize = 36,
    ylabelsize = 36,
    xticklabelsize = 30,
    yticklabelsize = 30,
)

xlims!(ax, 0.0, β)

scatter!(
    τ, Gτ_noisy,
    color = :blue, markersize = 9, label = "Noisy"
)

lines!(
    τ, Gτ_exact,
    linewidth = 7, alpha = 1.0, color = :black, linestyle = :dash, label = "Exact"
)

lines!(
    τ, Gτ_clean,
    linewidth = 4, color = :red, linestyle = :solid, label = "Denoised"
)
axislegend(ax, halign = :left, valign = :top, labelsize = 34)

fig

# Often noisy imaginary-time correlation data is generated by a quantum Monte Carlo simulation,
# or some other stochastic process, in which you have a set of samples or bins.
# In this case it is desirable to propagate errors through the denoising process to get a final
# error for the denoised correlation function. Below we will demonstrate two ways of doing this,
# one using the jackknife method and the other via bootstrap resampling.
#
# To start let us generated binned noisy synthetic ``G(\tau)`` data.

## Number of samples
N_bins = 32

## Generate synthetic binned data
Gτ_binned = hcat((Gτ_exact for b in 1:N_bins)...) + 0.08 * randn(length(τ), N_bins)

## Reduce noise for G(τ = 0)
@. Gτ_binned[1,:] = Gτ_exact[1] +  0.02 * randn()

## Ensure that noisy data still satisty the sum rule G(τ=0) + G(τ=β) = 1
@. Gτ_binned[end,:] = 1.0 - Gτ_binned[1,:]

Gτ_binned

# Before we demonstrate how to propagate errors through the denoising process,
# let us first calculate mean and error of our binned ``G(\tau)`` data first.
# Having calculated the mean and error of our synthetic binned ``G(\tau)`` data,
# let us also plot it really quickly as a reference.

## Calculate mean of binned data
Gτ_avg = vec(mean(Gτ_binned, dims=2))

## Calculate standard deviation of the mean of binned data
Gτ_err = vec(std(Gτ_binned, dims=2)) / sqrt(N_bins)

## Get the average plus or minus the error
Gτ_avg_lower = Gτ_avg - Gτ_err
Gτ_avg_upper = Gτ_avg + Gτ_err;

## Now we plot the average with the one standard deviation confidence interval
fig = Figure(
    size = (700, 700),
    fonts = (; regular= "CMU Serif"),
    figure_padding = 10
)

ax = Axis(fig[1, 1],
    aspect = 1,
    xlabel = L"\tau",
    ylabel = L"G_\sigma(\tau)",
    xlabelsize = 36,
    ylabelsize = 36,
    xticklabelsize = 30,
    yticklabelsize = 30,
)

b_err = band!(
    τ, Gτ_avg_lower, Gτ_avg_upper,
    color = (:red, 0.2)
)
translate!(b_err, 0, 0, 0.0)

l_err = lines!(
    τ, Gτ_exact,
    linewidth = 3, alpha = 1.0, color = :black, linestyle = :solid, label = "Exact"
)
translate!(l_err, 0, 0, 2.0)

l_avg = lines!(
    τ, Gτ_avg,
    linewidth = 3, color = :red, linestyle = :solid, label = "Mean"
)
translate!(l_avg, 0, 0, 1.0)

l_lower = lines!(
    τ, Gτ_avg_lower,
    linewidth = 2, alpha = 0.6, color = :red, linestyle = :solid
)
translate!(l_lower, 0, 0, 1.0)

l_upper = lines!(
    τ, Gτ_avg_upper,
    linewidth = 2, alpha = 0.6, color = :red, linestyle = :solid
)
translate!(l_upper, 0, 0, 1.0)

xlims!(ax, 0.0, β)

axislegend(ax, halign = :left, valign = :top, labelsize = 34)

fig

# Now let us denoise the binned ``G(\tau)`` while propagating the errors
# using the jackknife method.

## Denoise and propagate errors with jackknife
Gτ_jackknife_avg, Gτ_jackknife_err, Gτ_jackknife_cov = jackknife_hankel_correlation_cleaner(
    correlation_bins = Gτ_binned,
    sign_bins = ones(length(τ)),
    maxiter = 1000,
    tol= 1e-3,
    positive_curvature = true,
    fixed_endpoints = true,
    symmetric = false,
    covariance = true
)

## Get the average plus or minus the error
Gτ_jackknife_avg_lower = Gτ_jackknife_avg - Gτ_jackknife_err
Gτ_jackknife_avg_upper = Gτ_jackknife_avg + Gτ_jackknife_err;

## Now we plot the average with the one standard deviation confidence interval
fig = Figure(
    size = (700, 700),
    fonts = (; regular= "CMU Serif"),
    figure_padding = 10
)

ax = Axis(fig[1, 1],
    aspect = 1,
    xlabel = L"\tau",
    ylabel = L"G_\sigma(\tau)",
    xlabelsize = 36,
    ylabelsize = 36,
    xticklabelsize = 30,
    yticklabelsize = 30,
)

b_err = band!(
    τ, Gτ_jackknife_avg_lower, Gτ_jackknife_avg_upper,
    color = (:red, 0.2)
)
translate!(b_err, 0, 0, 0.0)

l_err = lines!(
    τ, Gτ_exact,
    linewidth = 3, alpha = 1.0, color = :black, linestyle = :solid, label = "Exact"
)
translate!(l_err, 0, 0, 2.0)

l_avg = lines!(
    τ, Gτ_jackknife_avg,
    linewidth = 3, color = :red, linestyle = :solid, label = "Denoised Jackknife Mean"
)
translate!(l_avg, 0, 0, 1.0)

l_lower = lines!(
    τ, Gτ_jackknife_avg_lower,
    linewidth = 2, alpha = 0.6, color = :red, linestyle = :solid
)
translate!(l_lower, 0, 0, 1.0)

l_upper = lines!(
    τ, Gτ_jackknife_avg_upper,
    linewidth = 2, alpha = 0.6, color = :red, linestyle = :solid
)
translate!(l_upper, 0, 0, 1.0)

xlims!(ax, 0.0, β)

axislegend(ax, halign = :left, valign = :top, labelsize = 34)

fig

# Lastly, let us denoise the binned ``G(\tau)`` but instead propagate error
# using bootstrap resampling.

## Denoise and propagate errors with bootstrap resampling
Gτ_bootstrap_avg, Gτ_bootstrap_err, Gτ_bootstrap_cov = bootstrap_hankel_correlation_cleaner(
    correlation_bins = Gτ_binned,
    sign_bins = ones(length(τ)),
    N_bootstrap = 100,
    maxiter = 1000,
    tol= 1e-3,
    positive_curvature = true,
    fixed_endpoints = true,
    symmetric = false,
    covariance = true
)

## Get the average plus or minus the error
Gτ_bootstrap_avg_lower = Gτ_bootstrap_avg - Gτ_bootstrap_err
Gτ_bootstrap_avg_upper = Gτ_bootstrap_avg + Gτ_bootstrap_err;

## Now we plot the average with the one standard deviation confidence interval
fig = Figure(
    size = (700, 700),
    fonts = (; regular= "CMU Serif"),
    figure_padding = 10
)

ax = Axis(fig[1, 1],
    aspect = 1,
    xlabel = L"\tau",
    ylabel = L"G_\sigma(\tau)",
    xlabelsize = 36,
    ylabelsize = 36,
    xticklabelsize = 30,
    yticklabelsize = 30,
)

b_err = band!(
    τ, Gτ_bootstrap_avg_lower, Gτ_bootstrap_avg_upper,
    color = (:red, 0.2)
)
translate!(b_err, 0, 0, 0.0)

l_err = lines!(
    τ, Gτ_exact,
    linewidth = 3, alpha = 1.0, color = :black, linestyle = :solid, label = "Exact"
)
translate!(l_err, 0, 0, 2.0)

l_avg = lines!(
    τ, Gτ_bootstrap_avg,
    linewidth = 3, color = :red, linestyle = :solid, label = "Denoised Bootstrap Mean"
)
translate!(l_avg, 0, 0, 1.0)

l_lower = lines!(
    τ, Gτ_bootstrap_avg_lower,
    linewidth = 2, alpha = 0.6, color = :red, linestyle = :solid
)
translate!(l_lower, 0, 0, 1.0)

l_upper = lines!(
    τ, Gτ_bootstrap_avg_upper,
    linewidth = 2, alpha = 0.6, color = :red, linestyle = :solid
)
translate!(l_upper, 0, 0, 1.0)

xlims!(ax, 0.0, β)

axislegend(ax, halign = :left, valign = :top, labelsize = 34)

fig