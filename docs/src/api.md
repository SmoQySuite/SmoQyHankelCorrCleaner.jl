# API

- [`hankel_correlation_cleaner`](@ref)
- [`hankel_correlation_cleaner!`](@ref)
- [`jackknife_hankel_correlation_cleaner`](@ref)
- [`jackknife_hankel_correlation_cleaner!`](@ref)
- [`bootstrap_hankel_correlation_cleaner`](@ref)
- [`bootstrap_hankel_correlation_cleaner!`](@ref)

```@docs
hankel_correlation_cleaner
hankel_correlation_cleaner!(::AbstractVector{T}) where {T<:AbstractFloat}
hankel_correlation_cleaner!(::AbstractVector{T}, ::AbstractVector{T}) where {T<:AbstractFloat}
jackknife_hankel_correlation_cleaner
jackknife_hankel_correlation_cleaner!
bootstrap_hankel_correlation_cleaner
bootstrap_hankel_correlation_cleaner!
```

## Developer API

```@docs
SmoQyHankelCorrCleaner.dykstra!
SmoQyHankelCorrCleaner.project_psd!
SmoQyHankelCorrCleaner.project_hankel!
SmoQyHankelCorrCleaner.project_positive_curvature!
SmoQyHankelCorrCleaner.init_hankel_matrix
SmoQyHankelCorrCleaner.init_hankel_matrix!
```