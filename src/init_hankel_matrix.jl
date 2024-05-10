@doc raw"""
    init_hankel_matrix(
        correlations::AbstractVector{T}
    ) where {T<:AbstractFloat}

Allocate, initialize and return the Hankel matrix based on the vector of `correlations`.
"""
function init_hankel_matrix(
    correlations::AbstractVector{T}
) where {T<:AbstractFloat}

    # length of imaginary time axis plus one
    Lτp1 = length(correlations)

    # half the length of the imaginary time axis
    Lτo2 = (Lτp1-1)÷2
    Lτo2p1 = Lτo2 + 1

    # allocate hankel matrix
    H = zeros(T, Lτo2p1, Lτo2p1)

    # initialize hankel matrix
    init_hankel_matrix!(H, correlations)

    return H
end

@doc raw""""
    init_hankel_matrix!(
        H::AbstractMatrix{T},
        correaltions::AbstractVector{T}
    ) where {T<:AbstractFloat}

Initialize the Hankel matrix `H` in-place based on the vector of `correlations`.
"""
function init_hankel_matrix!(
    H::AbstractMatrix{T},
    correaltions::AbstractVector{T}
) where {T<:AbstractFloat}

    # iterate over columns of Hankel matrix
    @simd for j in axes(H,2)
        # iterate over rows of Hankel matrix
        for i in axes(H,1)
            H[i,j] = correaltions[i+j-1]
        end
    end

    return nothing
end