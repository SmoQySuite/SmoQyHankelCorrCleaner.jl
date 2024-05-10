@testitem "Test hankel_correlation_cleaner()" begin

    # single-site hubbard model green's function
    function eval_Ghub(τ, β, U, μ)

        Zhub = exp(-β*U/4) + 2*exp(β*(U/4+μ)) + exp(-β*(U/4-2*μ))
        Gτ = exp(τ*U/2 + τ*μ - β*U/4) + exp(-τ*U/2 + τ*μ + β*U/4 + β*μ)
        Gτ = Gτ / Zhub
        return Gτ
    end
    
    U  = 5.00
    μ  = 2.50
    β  = 5.00
    Δτ = 0.05
    τ  = collect(range(start = 0.0, stop = β, step = Δτ))
    Lτ = length(τ)
    G_exact = eval_Ghub.(τ, β, U, μ)
    G_noisy = G_exact + 0.01 * randn(Lτ)
    G_noisy[1] = G_exact[1] + 0.002 * randn()
    G_noisy[end] = 1.0 - G_noisy[1] # enforce sum rule G(τ=0) + G(τ=β) = 1.0
    G_clean, iter, err = hankel_correlation_cleaner(
        G_noisy;
        maxiter = 100,
        tol = 1e-3,
        fixed_endpoints = true,
        positive_curvature = true,
        symmetric = false,
        verbose = false
    )

    # test that the endpoint did not change
    @test G_clean[1]   == G_noisy[1]
    @test G_clean[end] == G_noisy[end]

    # check that the curvature is strictly positive
    @test all(g -> g > -1e-9, diff(diff(G_clean)))

    # check that the error is reduced
    G_clean_err = sum((G_clean[i]-G_exact[i])^2 for i in eachindex(G_clean))
    G_noisy_err = sum((G_noisy[i]-G_exact[i])^2 for i in eachindex(G_noisy))
    @test G_noisy_err > G_clean_err
end