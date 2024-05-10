@testitem "Test bootstrap_hankel_correlation_cleaner()" begin

    using Statistics

    # single-site hubbard model green's function
    function eval_Ghub(τ, β, U, μ)

        Zhub = exp(-β*U/4) + 2*exp(β*(U/4+μ)) + exp(-β*(U/4-2*μ))
        Gτ = exp(τ*U/2 + τ*μ - β*U/4) + exp(-τ*U/2 + τ*μ + β*U/4 + β*μ)
        Gτ = Gτ / Zhub
        return Gτ
    end
    
    # calculate single-site hubbard model greens function
    U  = 5.00
    μ  = 2.50
    β  = 5.00
    Δτ = 0.05
    τ  = collect(range(start = 0.0, stop = β, step = Δτ))
    Lτ = length(τ)
    G_exact = eval_Ghub.(τ, β, U, μ)

    # generate synthetic noisy binned greens function data
    N_bins = 32
    G_binned = hcat((G_exact for b in 1:N_bins)...) + 0.01 * randn(Lτ, N_bins)
    @views @. G_binned[1,:] = G_binned[1,:] + 0.002 * randn()
    @views @. G_binned[end,:] = 1 - G_binned[1,:] # enforce sum rule G(τ=0) + G(τ=β) = 1

    # get mean of noisy binned greens function data
    G_mean = mean(G_binned, dims=2)

    # perform jackknife correlation cleaning
    G_clean_mean, G_clean_err, G_cov = bootstrap_hankel_correlation_cleaner(
        correlation_bins = G_binned,
        sign_bins = ones(N_bins),
        N_bootstrap = N_bins,
        maxiter = 100,
        tol = 1e-3,
        positive_curvature = true,
        fixed_endpoints = true,
        symmetric = false,
        covariance = true
    )

    # make sure sum rule still statisfied
    @test G_clean_mean[1] + G_clean_mean[end] ≈ 1.0

    # check that the curvature is strictly positive
    @test all(g -> g > -1e-9, diff(diff(G_clean_mean)))

    # check that the error is reduced
    G_clean_err = sum((G_clean_mean[i]-G_exact[i])^2 for i in eachindex(G_clean_mean))
    G_noisy_err = sum((G_mean[i]-G_exact[i])^2 for i in eachindex(G_mean))
    @test G_noisy_err > G_clean_err
end