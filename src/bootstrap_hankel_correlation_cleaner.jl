@doc raw"""
    bootstrap_hankel_correlation_cleaner(;
        # KEYWORD ARGUMENTS
        correlation_bins::AbstractMatrix{T},
        sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2)),
        N_bootstrap::Int = 100,
        maxiter::Int = 1000,
        tol::T = 1e-4,
        positive_curvature::Bool = false,
        fixed_endpoints::Bool = true,
        symmetric::Bool = false,
        covariance::Bool = false
    ) where {T<:AbstractFloat}

Denoise binned imaginary-time correlation data, propagating the error via bootstrapping.

# Arguments

- `correlation_bins::AbstractMatrix{T}`: Binned imaginary time correlations where the columns correspond to the bins and the rows to the imaginary time slices.

# Keyword Arguments

- `sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2))`: For data generated via a quantum Monte Carlo simulation with a sign problem, the average sign associated with each bin.
- `N_bootstrap::Int = 100`: Number of bootstrap samples used to calculate the error.
- `maxiter::T`: Maximum number of iteration used in Dykstra's algorithm.
- `tol::T`: Tolerance threshold used in Dykstra's algorithm.
- `positive_curvature::Bool = false`: Whether to project onto a matrix with average anti-diagonals that have strictly positive curvature.
- `fixed_endpoints::Bool = true`: Whether to fix the correlation values at ``\tau = 0`` and ``\tau = \beta``.
- `symmetric::Bool = false`: Whether the imaginary-time correlation data is symmetric about ``\tau = \beta/2``. 
- `covariance::Bool = false`: If true, this method also returns the covariance matrix.
"""
function bootstrap_hankel_correlation_cleaner(;
    # KEYWORD ARGUMENTS
    correlation_bins::AbstractMatrix{T},
    sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2)),
    N_bootstrap::Int = 100,
    maxiter::Int = 1000,
    tol::T = 1e-4,
    positive_curvature::Bool = false,
    fixed_endpoints::Bool = true,
    symmetric::Bool = false,
    covariance::Bool = false
) where {T<:AbstractFloat}

    clean_correlations_mean = zeros(T, size(correlation_bins,1))
    clean_correlations_error = zeros(T, size(correlation_bins,1))

    covariance_matrix = bootstrap_hankel_correlation_cleaner!(
        clean_correlations_mean = clean_correlations_mean,
        clean_correlations_error = clean_correlations_error,
        correlation_bins = correlation_bins,
        sign_bins = sign_bins,
        N_bootstrap = N_bootstrap,
        maxiter = maxiter,
        tol = tol,
        positive_curvature = positive_curvature,
        fixed_endpoints = fixed_endpoints,
        symmetric = symmetric,
        covariance = covariance
    )

    if covariance
        vals = (clean_correlations_mean, clean_correlations_error, covariance_matrix)
    else
        vals = (clean_correlations_mean, clean_correlations_error)
    end

    return vals
end

@doc raw"""
    bootstrap_hankel_correlation_cleaner!(;
        # KEYWORD ARGUMENTS
        clean_correlations_mean::AbstractVector{T},
        clean_correlations_error::AbstractVector{T},
        correlation_bins::AbstractMatrix{T},
        sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2)),
        N_bootstrap::Int = 100,
        maxiter::Int = 1000,
        tol::T = 1e-4,
        positive_curvature::Bool = false,
        fixed_endpoints::Bool = true,
        symmetric::Bool = false,
        covariance::Bool = false
    ) where {T<:AbstractFloat}

Denoise binned imaginary-time correlation data, propagating the error via bootstrapping.

# Arguments

- `clean_correlations_mean::AbstractVector{T}`: Vector to contain denoised imaginary time correlations.
- `clean_correlations_error::AbstractVector{T}`: Vector to contain error associated with denoised imaginary time correlations.
- `correlation_bins::AbstractMatrix{T}`: Binned imaginary time correlations where the columns correspond to the bins and the rows to the imaginary time slices.

# Keyword Arguments

- `sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2))`: For data generated via a quantum Monte Carlo simulation with a sign problem, the average sign associated with each bin.
- `N_bootstrap::Int = 100`: Number of bootstrap samples used to calculate the error.
- `maxiter::T`: Maximum number of iteration used in Dykstra's algorithm.
- `tol::T`: Tolerance threshold used in Dykstra's algorithm.
- `positive_curvature::Bool = false`: Whether to project onto a matrix with average anti-diagonals that have strictly positive curvature.
- `fixed_endpoints::Bool = true`: Whether to fix the correlation values at ``\tau = 0`` and ``\tau = \beta``.
- `symmetric::Bool = false`: Whether the imaginary-time correlation data is symmetric about ``\tau = \beta/2``. 
- `covariance::Bool = false`: If true, this method also returns the covariance matrix.
"""
function bootstrap_hankel_correlation_cleaner!(;
    # KEYWORD ARGUMENTS
    clean_correlations_mean::AbstractVector{T},
    clean_correlations_error::AbstractVector{T},
    correlation_bins::AbstractMatrix{T},
    sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2)),
    N_bootstrap::Int = 100,
    maxiter::Int = 1000,
    tol::T = 1e-4,
    positive_curvature::Bool = false,
    fixed_endpoints::Bool = true,
    symmetric::Bool = false,
    covariance::Bool = false
) where {T<:AbstractFloat}

    # get length of imaginary time axis plus one and number of samples
    (Lτp1, N_bin) = size(correlation_bins)

    # vector to contain single bootstrap sample mean
    bootstrap_sample = zeros(T, Lτp1)

    # to track bootstrap variance and full covariance
    bootstrap_cov_tracker = OnlineStats.CovMatrix(Lτp1)

    # iterate over number of bootstrap samples
    for n in 1:N_bootstrap

        # initialize bootstrap sample to zero
        fill!(bootstrap_sample, 0)

        # initialize boostrap sign to zero
        bootstrap_sign = zero(T)

        # iterate over bins
        for i in 1:N_bin

            # randomly sample bin
            bin = rand(1:N_bin)

            # record correlation bin
            @views @. bootstrap_sample += correlation_bins[:,bin]

            # record sign bin
            bootstrap_sign += sign_bins[bin]
        end

        # calculate final bootstrap sample mean, accounting for sign
        @. bootstrap_sample = bootstrap_sample / bootstrap_sign

        # denoise the bootstrap sample
        hankel_correlation_cleaner!(
            bootstrap_sample,
            bootstrap_sample,
            maxiter = maxiter,
            tol = tol,
            positive_curvature = positive_curvature,
            fixed_endpoints = fixed_endpoints,
            symmetric = symmetric,
            verbose = false
        )

        # update bootstrap mean
        @. clean_correlations_mean = clean_correlations_mean + (bootstrap_sample - clean_correlations_mean)/n

        # update variances and covariance
        OnlineStats.fit!(bootstrap_cov_tracker, bootstrap_sample)
    end

    # calculate jackknife error
    bootstrap_var = var(bootstrap_cov_tracker, corrected = false)
    @. clean_correlations_error = sqrt( N_bin/(N_bin - 1) * bootstrap_var )

    # calculate bootstrap covariance matrix
    bootstrap_cov = nothing
    if covariance
        bootstrap_cov = N_bin/(N_bin - 1) * cov(bootstrap_cov_tracker, corrected = false)
    end

    return bootstrap_cov
end