@doc raw"""
    hankel_correlation_cleaner(
        noisy_correlations::AbstractVector{T};
        # KEYWORD ARGUMENTS
        maxiter::Int,
        tol::T,
        positive_curvature::Bool = true,
        fixed_endpoints::Bool = true,
        symmetric::Bool = false,
        verbose::Bool = false
    ) where {T<:AbstractFloat}

Denoise imaginary-time correlation data using the Hankel projection method.
Returns the tuple `(clean_correlations iter, err)` where
`clean_correlation::Vector{T}` is a vector of the cleaned correlations,
`iter::Int` is the number of iterations and `err::T` is the final error.

# Arguments

- `noisy_correlations::AbstractVector{T}`: Vector of imaginary time correlation data to be cleaned.

# Keyword Arguments

- `maxiter::T`: Maximum number of iteration used in Dykstra's algorithm.
- `tol::T`: Tolerance threshold used in Dykstra's algorithm.
- `positive_curvature::Bool = false`: Whether to project onto a matrix with average anti-diagonals that have strictly positive curvature.
- `fixed_endpoints::Bool = true`: Whether to fix the correlation values at ``\tau = 0`` and ``\tau = \beta``.
- `symmetric::Bool = false`: Whether the imaginary-time correlation data is symmetric about ``\tau = \beta/2``. 
- `verbose::Bool = false`: Whether to print statements tracking the convergence of the denoising process.
"""
function hankel_correlation_cleaner(
    noisy_correlations::AbstractVector{T};
    # KEYWORD ARGUMENTS
    maxiter::Int,
    tol::T,
    positive_curvature::Bool = true,
    fixed_endpoints::Bool = true,
    symmetric::Bool = false,
    verbose::Bool = false
) where {T<:AbstractFloat}

    clean_correlations = zero(noisy_correlations)
    iter, err = hankel_correlation_cleaner!(
        clean_correlations, noisy_correlations,
        maxiter = maxiter,
        tol = tol,
        positive_curvature = positive_curvature,
        fixed_endpoints = fixed_endpoints,
        symmetric = symmetric,
        verbose = verbose
    )

    return clean_correlations, iter, err
end


@doc raw"""
    hankel_correlation_cleaner!(
        correlations::AbstractVector{T};
        # KEYWORD ARGUMENTS
        maxiter::Int,
        tol::T,
        positive_curvature::Bool = true,
        fixed_endpoints::Bool = true,
        symmetric::Bool = false,
        verbose::Bool = false
    ) where {T<:AbstractFloat}

Denoise imaginary-time correlation data using the Hankel projection method.
Also return the tuple `(iter, err)` where `iter::Int` is the number of iterations
and `err::T` is the final error.

# Arguments

- `correlations::AbstractVector{T}`: Vector of imaginary time correlation data to be cleaned.

# Keyword Arguments

- `maxiter::T`: Maximum number of iteration used in Dykstra's algorithm.
- `tol::T`: Tolerance threshold used in Dykstra's algorithm.
- `positive_curvature::Bool = false`: Whether to project onto a matrix with average anti-diagonals that have strictly positive curvature.
- `fixed_endpoints::Bool = true`: Whether to fix the correlation values at ``\tau = 0`` and ``\tau = \beta``.
- `symmetric::Bool = false`: Whether the imaginary-time correlation data is symmetric about ``\tau = \beta/2``. 
- `verbose::Bool = false`: Whether to print statements tracking the convergence of the denoising process.
"""
function hankel_correlation_cleaner!(
    correlations::AbstractVector{T};
    # KEYWORD ARGUMENTS
    maxiter::Int,
    tol::T,
    positive_curvature::Bool = true,
    fixed_endpoints::Bool = true,
    symmetric::Bool = false,
    verbose::Bool = false
) where {T<:AbstractFloat}

    return hankel_correlation_cleaner!(
        correlations, correlations,
        maxiter = maxiter,
        tol = tol,
        positive_curvature = positive_curvature,
        fixed_endpoints = fixed_endpoints,
        symmetric = symmetric,
        verbose = verbose
    )
end


@doc raw"""
    hankel_correlation_cleaner!(
        clean_correlations::AbstractVector{T},
        noisy_correlations::AbstractVector{T};
        # KEYWORD ARGUMENTS
        maxiter::Int,
        tol::T,
        positive_curvature::Bool = true,
        fixed_endpoints::Bool = true,
        symmetric::Bool = false,
        verbose::Bool = false
    ) where {T<:AbstractFloat}

Denoise imaginary-time correlation data using the Hankel projection method.
Also return the tuple `(iter, err)` where `iter::Int` is the number of iterations
and `err::T` is the final error.

# Arguments

- `clean_correlations::AbstractVector{T}`: Vector to contain cleaned/denoised imaginary-time correlation data.
- `noisy_correlations::AbstractVector{T}`: Noisy imaginary-time correlation data.

# Keyword Arguments

- `maxiter::T`: Maximum number of iteration used in Dykstra's algorithm.
- `tol::T`: Tolerance threshold used in Dykstra's algorithm.
- `positive_curvature::Bool = false`: Whether to project onto a matrix with average anti-diagonals that have strictly positive curvature.
- `fixed_endpoints::Bool = true`: Whether to fix the correlation values at ``\tau = 0`` and ``\tau = \beta``.
- `symmetric::Bool = false`: Whether the imaginary-time correlation data is symmetric about ``\tau = \beta/2``. 
- `verbose::Bool = false`: Whether to print statements tracking the convergence of the denoising process.
"""
function hankel_correlation_cleaner!(
    clean_correlations::AbstractVector{T},
    noisy_correlations::AbstractVector{T};
    # KEYWORD ARGUMENTS
    maxiter::Int,
    tol::T,
    positive_curvature::Bool = true,
    fixed_endpoints::Bool = true,
    symmetric::Bool = false,
    verbose::Bool = false
) where {T<:AbstractFloat}


    # construct initial Hankel matrix based on correlationd data
    H = init_hankel_matrix(noisy_correlations)

    # initialize storage matrix
    A  = zero(H)
    A′ = @view A[2:end, 2:end]

    # initialize vector to contain projection functions
    projections! = Function[]

    # initialize vector to contain projection names
    projection_names = String[]

    # initilaize partial PSD projection function
    proj_psd_partial!(X) = project_psd!(view(X,2:size(X,1),2:size(X,2)), A′)
    push!(projections!, proj_psd_partial!)
    if verbose
        push!(projection_names, "Partial PSD Projection")
    end

    # initialize full PSD projection function
    proj_psd_full!(X) = project_psd!(X, A)
    push!(projections!, proj_psd_full!)
    if verbose
        push!(projection_names, "Full PSD Projection")
    end

    # if fixing the values of C(τ=0) and C(τ=β) to initial values
    if fixed_endpoints
        endpoints = (noisy_correlations[1], noisy_correlations[end])
    else
        endpoints = (NaN, NaN)
    end

    # initialize hankel projection function
    proj_hankel!(X) = project_hankel!(X, endpoints)
    push!(projections!, proj_hankel!)
    if verbose
        push!(projection_names, "Hankel Projection")
    end

    # initialize positive curvature projection
    if positive_curvature
        # model = define_project_positive_curvature(H, endpoints)
        # proj_pos_curv!(X) = project_positive_curvature!(X, model)
        proj_pos_curv!(X) = project_positive_curvature!(X, endpoints)
        push!(projections!, proj_pos_curv!)
        if verbose
            push!(projection_names, "Positive Curvature Projection")
        end
    end

    # if symmetric correlation function about τ = β/2
    if symmetric
        push!(projections!, project_symmetrized!)
        if verbose
            push!(projection_names, "Symmetrization Projection")
        end
    end

    # perform Dykstra's algorithm with given projections
    iter, err = dykstra!(
        H = H,
        projections! = projections!,
        tol = tol,
        maxiter = maxiter,
        verbose = verbose,
        projection_names = projection_names
    )

    # record the cleaned correlations by averaging anti-diagonals
    Lτo2p1 = size(H,1)
    fill!(clean_correlations, 0)
    @simd for j in 1:Lτo2p1
        for i in 1:Lτo2p1
            clean_correlations[i+j-1] += H[i, j]
        end
    end

    # normalize the clean correlations
    @simd for i in 1:Lτo2p1
        clean_correlations[i] /= i
    end
    @simd for j in 2:Lτo2p1
        clean_correlations[Lτo2p1+j-1] /= (Lτo2p1-j+1)
    end

    return iter, err
end