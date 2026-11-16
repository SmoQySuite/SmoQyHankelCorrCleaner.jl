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

@testitem "Test hankel_correlation_cleaner() PSD constraints" begin

    using LinearAlgebra
    using Random
    import SmoQyHankelCorrCleaner: init_hankel_matrix

    # single-site hubbard model green's function
    function eval_Ghub(τ, β, U, μ)

        Zhub = exp(-β*U/4) + 2*exp(β*(U/4+μ)) + exp(-β*(U/4-2*μ))
        Gτ = exp(τ*U/2 + τ*μ - β*U/4) + exp(-τ*U/2 + τ*μ + β*U/4 + β*μ)
        Gτ = Gτ / Zhub
        return Gτ
    end

    Random.seed!(1234)
    U  = 5.00
    μ  = 2.50
    β  = 5.00
    Δτ = 0.10
    τ  = collect(range(start = 0.0, stop = β, step = Δτ))
    G_exact = eval_Ghub.(τ, β, U, μ)
    G_noisy = G_exact + 0.01 * randn(length(τ))
    G_clean, iter, err = hankel_correlation_cleaner(
        G_noisy;
        maxiter = 5000,
        tol = 1e-8,
        fixed_endpoints = true,
        positive_curvature = false,
        symmetric = false,
        verbose = false
    )
    @test err < 1e-8

    # both the full Hankel matrix H[i,j] = G(τ[i+j-1]) and the shifted Hankel matrix
    # H′[i,j] = G(τ[i+j]), obtained by dropping the first column and last row of H,
    # must be positive semidefinite
    H = init_hankel_matrix(G_clean)
    H′ = H[1:end-1, 2:end]
    @test minimum(eigvals(Symmetric(H)))  > -1e-8
    @test minimum(eigvals(Symmetric(H′))) > -1e-8

    # for reference, the noisy input violates both constraints
    H_noisy = init_hankel_matrix(G_noisy)
    @test minimum(eigvals(Symmetric(H_noisy[1:end-1, 2:end]))) < -1e-4
end


@testitem "Test project_positive_curvature!() is an exact Frobenius projection" begin

    using Random
    import SmoQyHankelCorrCleaner: project_positive_curvature!, antidiagonal_average, JuMP, HiGHS

    # brute-force projection: optimize over all N² entries of X directly
    function brute_force_projection(H, endpoints)
        N = size(H, 1)
        L = 2N - 1
        model = JuMP.Model(HiGHS.Optimizer)
        JuMP.set_silent(model)
        JuMP.@variable(model, X[1:N, 1:N])
        avg(l) = sum(X[i, l-i+1] for i in max(1, l-N+1):min(N, l)) / min(l, L-l+1)
        JuMP.@constraint(model, [l in 2:L-1], avg(l+1) + avg(l-1) - 2avg(l) ≥ 0)
        JuMP.@constraint(model, X[1, 1] == endpoints[1])
        JuMP.@constraint(model, X[N, N] == endpoints[2])
        JuMP.@objective(model, Min, sum((X[i, j] - H[i, j])^2 for i in 1:N, j in 1:N))
        JuMP.optimize!(model)
        return JuMP.value.(X)
    end

    Random.seed!(3)
    N = 6
    H = randn(N, N) # deliberately not Hankel, as happens during Dykstra iterations
    endpoints = (0.3, 0.1)
    X = copy(H)
    project_positive_curvature!(X, endpoints)
    @test maximum(abs, X - brute_force_projection(H, endpoints)) < 1e-5

    # the result satisfies the constraints
    a = [antidiagonal_average(X, l) for l in 1:2N-1]
    @test all(≥(-1e-9), diff(diff(a)))
    @test X[1, 1] ≈ endpoints[1]
    @test X[N, N] ≈ endpoints[2]
end


@testitem "Test hankel_correlation_cleaner() with symmetric = true" begin

    using Random
    import SmoQyHankelCorrCleaner: project_symmetrized!, antidiagonal_average

    # the symmetrization projection is idempotent and symmetrizes anti-diagonal averages
    Random.seed!(5)
    N = 7
    H = randn(N, N)
    project_symmetrized!(H)
    a = [antidiagonal_average(H, l) for l in 1:2N-1]
    @test a ≈ reverse(a)
    H2 = copy(H)
    project_symmetrized!(H2)
    @test H2 ≈ H

    # particle-hole symmetric single-site Hubbard model: G(τ) = G(β-τ)
    function eval_Ghub(τ, β, U, μ)
        Zhub = exp(-β*U/4) + 2*exp(β*(U/4+μ)) + exp(-β*(U/4-2*μ))
        Gτ = exp(τ*U/2 + τ*μ - β*U/4) + exp(-τ*U/2 + τ*μ + β*U/4 + β*μ)
        return Gτ / Zhub
    end
    β  = 5.00
    τ  = collect(range(start = 0.0, stop = β, step = 0.05))
    G_exact = eval_Ghub.(τ, β, 5.0, 0.0)
    G_noisy = G_exact + 0.01 * randn(length(τ))
    G_noisy[1] = G_exact[1] + 0.002 * randn()
    G_noisy[end] = G_noisy[1]
    G_clean, iter, err = hankel_correlation_cleaner(
        G_noisy;
        maxiter = 100,
        tol = 1e-3,
        fixed_endpoints = true,
        positive_curvature = true,
        symmetric = true
    )
    @test G_clean ≈ reverse(G_clean) atol = 1e-3
    @test sum(abs2, G_clean - G_exact) < sum(abs2, G_noisy - G_exact)
end


@testitem "Test even-length input is rejected" begin
    @test_throws ArgumentError hankel_correlation_cleaner(rand(10); maxiter = 10, tol = 1e-3)
end
