@doc raw"""
    jackknife_hankel_correlation_cleaner(;
        # KEYWORD ARGUMENTS
        correlation_bins::AbstractMatrix{T},
        sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2)),
        maxiter::Int = 1000,
        tol::T = 1e-4,
        positive_curvature::Bool = false,
        fixed_endpoints::Bool = true,
        symmetric::Bool = false,
        covariance::Bool = false
    ) where {T<:AbstractFloat}

Denoise binned imaginary-time correlation data, propagating the error using the
jackknife algorithm. Returns the mean and error of the denoised imaginary-time correlation data.

# Arguments

- `clean_correlations_mean::AbstractVector{T}`: Vector to contain denoised imaginary-time correlations.
- `clean_correlations_error::AbstractVector{T}`: Vector to contain error associated with denoised imaginary-time correlations.
- `correlation_bins::AbstractMatrix{T}`: Binned imaginary-time correlations where the columns correspond to the bins and the rows to the imaginary time slices.

# Keyword Arguments

- `sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2))`: For data generated via a quantum Monte Carlo simulation with a sign problem, the average sign associated with each bin.
- `maxiter::T`: Maximum number of iteration used in Dykstra's algorithm.
- `tol::T`: Tolerance threshold used in Dykstra's algorithm.
- `positive_curvature::Bool = false`: Whether to project onto a matrix with average anti-diagonals that have strictly positive curvature.
- `fixed_endpoints::Bool = true`: Whether to fix the correlation values at ``\tau = 0`` and ``\tau = \beta``.
- `symmetric::Bool = false`: Whether the imaginary-time correlation data is symmetric about ``\tau = \beta/2``. 
- `covariance::Bool = false`: If true, this method also returns the covariance matrix.
"""
function jackknife_hankel_correlation_cleaner(;
    # KEYWORD ARGUMENTS
    correlation_bins::AbstractMatrix{T},
    sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2)),
    maxiter::Int = 1000,
    tol::T = 1e-4,
    positive_curvature::Bool = false,
    fixed_endpoints::Bool = true,
    symmetric::Bool = false,
    covariance::Bool = false
) where {T<:AbstractFloat}

    clean_correlations_mean = zeros(T, size(correlation_bins,1))
    clean_correlations_error = zeros(T, size(correlation_bins,1))

    covariance_matrix = jackknife_hankel_correlation_cleaner!(
        clean_correlations_mean = clean_correlations_mean,
        clean_correlations_error = clean_correlations_error,
        correlation_bins = correlation_bins,
        sign_bins = sign_bins,
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
    jackknife_hankel_correlation_cleaner!(
        clean_correlations_mean::AbstractVector{T},
        clean_correlations_error::AbstractVector{T},
        correlation_bins::AbstractMatrix{T};
        # KEYWORD ARGUMENTS
        sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2)),
        maxiter::Int = 1000,
        tol::T = 1e-4,
        positive_curvature::Bool = false,
        fixed_endpoints::Bool = true,
        symmetric::Bool = false,
        covariance::Bool = false
    ) where {T<:AbstractFloat}

Denoise binned imaginary-time correlation data, propagating the error using the
jackknife algorithm.

# Arguments

- `clean_correlations_mean::AbstractVector{T}`: Vector to contain denoised imaginary time correlations.
- `clean_correlations_error::AbstractVector{T}`: Vector to contain error associated with denoised imaginary time correlations.
- `correlation_bins::AbstractMatrix{T}`: Binned imaginary time correlations where the columns correspond to the bins and the rows to the imaginary time slices.

# Keyword Arguments

- `sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2))`: For data generated via a quantum Monte Carlo simulation with a sign problem, the average sign associated with each bin.
- `maxiter::T`: Maximum number of iteration used in Dykstra's algorithm.
- `tol::T`: Tolerance threshold used in Dykstra's algorithm.
- `positive_curvature::Bool = false`: Whether to project onto a matrix with average anti-diagonals that have strictly positive curvature.
- `fixed_endpoints::Bool = true`: Whether to fix the correlation values at ``\tau = 0`` and ``\tau = \beta``.
- `symmetric::Bool = false`: Whether the imaginary-time correlation data is symmetric about ``\tau = \beta/2``. 
- `covariance::Bool = false`: If true, this method also returns the covariance matrix.
"""
function jackknife_hankel_correlation_cleaner!(;
    # KEYWORD ARGUMENTS
    clean_correlations_mean::AbstractVector{T},
    clean_correlations_error::AbstractVector{T},
    correlation_bins::AbstractMatrix{T},
    sign_bins::AbstractVector{T} = ones(eltype(correlation_bins), size(correlation_bins, 2)),
    maxiter::Int = 1000,
    tol::T = 1e-4,
    positive_curvature::Bool = false,
    fixed_endpoints::Bool = true,
    symmetric::Bool = false,
    covariance::Bool = false
) where {T<:AbstractFloat}

    # get length of imaginary time axis plus one and number of samples
    (Lτp1, N) = size(correlation_bins)

    # calculate mean correlations
    mean_correlations = zeros(T, Lτp1)
    mean!(mean_correlations, correlation_bins)

    # calculate the average sign
    mean_sign = mean(sign_bins)

    # to contain jackknife samples
    jackknife_samples = zero(correlation_bins)

    # construct jackknife samples
    for i in axes(correlation_bins, 2)
        # calculate the jackknife sign by updating the mean sign to reflect removing the
        # average sign associated with the i'th bin
        jackknife_sign = mean_sign + (mean_sign - sign_bins[i])/(N-1)
        # calculate the i'th jackknife sample by updating the mean with respect to removing the i'th bin
        correlation_sample = @view correlation_bins[:,i]
        jackknife_sample = @view jackknife_samples[:,i]
        @. jackknife_sample = (mean_correlations + (mean_correlations - correlation_sample)/(N-1)) / jackknife_sign
    end

    # go through and denoise each jackknife sample
    for i in axes(jackknife_samples, 2)
        jackknife_sample = @view jackknife_samples[:, i]
        hankel_correlation_cleaner!(
            jackknife_sample,
            jackknife_sample;
            maxiter = maxiter,
            tol = tol,
            positive_curvature = positive_curvature,
            fixed_endpoints = fixed_endpoints,
            symmetric = symmetric,
            verbose = false
        )
    end

    # calculate the jackknife mean
    mean!(clean_correlations_mean, jackknife_samples)

    # calculate the jackknife error
    for l in eachindex(clean_correlations_mean)
        jackknife_samples_l = @view jackknife_samples[l,:]
        clean_correlations_error[l] = sqrt((N-1)*varm(jackknife_samples_l, clean_correlations_mean[l], corrected = false))
    end

    # calculate the full covariance matrix
    jackknife_cov = nothing
    if covariance
        jackknife_cov = (N-1) * cov(jackknife_samples, dims = 2, corrected = false)
    end

    return jackknife_cov
end