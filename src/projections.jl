@doc raw"""
    project_positive_curvature!(
        H::AbstractMatrix{T},
        endpoints::Tuple{T,T} = (NaN, NaN)
    ) where {T<:AbstractFloat}

Project the matrix `H` onto the nearest matrix, as measured by the Frobenius norm, that
has anti-diagonal average values that have strictly positive curvature i.e. if you take the second
finite difference of anti-diagonal average values the result values are strictly positive.
If the `endpoints` are real numbers, then the `H[1,1]` and `H[end:end]` set set equal to the
passed values.
"""
function project_positive_curvature!(
    H::AbstractMatrix{T},
    endpoints::Tuple{T,T} = (NaN, NaN)
) where {T<:AbstractFloat}

    # get size of Hankel matrix
    N = size(H,1)
    # get number of imaginary time slices
    L = 2*N-1
    # calculate antidiagonal averates
    target = [antidiagonal_average(H,l) for l in 1:L]
    # initialize model
    model = JuMP.Model(HiGHS.Optimizer)
    # setting strings off and silen on for efficiency
    JuMP.set_string_names_on_creation(model, false)
    JuMP.set_silent(model)
    # define variables to represent correlations that get fitted
    JuMP.@variable(model, fit[1:L])
    # enforce positive curvature
    JuMP.@constraint(model, [i in 2:L-1], fit[i+1] + fit[i-1] - 2.0 * fit[i] ≥ 0.0)
    # fix endpoints of correlation function
    if isfinite(endpoints[1]) && isfinite(endpoints[2])
        JuMP.@constraint(model, fit[1] == endpoints[1])
        JuMP.@constraint(model, fit[L] == endpoints[2])
    end
    # minimize Frobenius norm between fitted Hankel matrix and target matrix
    JuMP.@objective(model, Min, sum((fit[i] - target[i])^2 for i in 1:L))
    # find the solution
    JuMP.optimize!(model)
    # get fit
    result = JuMP.value.(fit)
    # record solution
    for j in axes(H,2)
        for i in axes(H,1)
            l = i+j-1
            H[i,j] += result[l] - target[l]
        end
    end

    return nothing
end

# calculate the average of the l'th anti-diagonal of the matrix H
function antidiagonal_average(
    H::AbstractMatrix{T},
    l::Int
) where {T<:AbstractFloat}

    N = size(H, 1)
    L = 2 * N - 1
    n = min(l, L-l+1)
    h = zero(T)
    for k in 1:n
        i = min(l, N) - k + 1
        j = max(0, l-N) + k
        h += H[i,j] / n
    end

    return h
end


@doc raw"""
    project_hankel!(
        H::AbstractMattrix{T},
        endpoints::Tuple{T,T} = (NaN, NaN)
    ) where {T<:AbstractFloat}

Modifying the matrix `H` in-place, project it onto the closest Hankel matrix as defined
by the Frobenius norm, by replacing each anti-diagonal by its average value.
If finite `endpoints` are passed, then the values of `H[1,1]` and `H[end,end]` are fixed
to those values respectively.
"""
function project_hankel!(
    H::AbstractMatrix{T},
    endpoints::Tuple{T,T} = (NaN, NaN)
) where {T<:AbstractFloat}

    Lτo2p1 = size(H, 1)
    Lτo2   = Lτo2p1 - 1

    # sum anti-diagonals, writing the sums to the first column and bottom row
    @simd for j in 2:Lτo2p1
        for i in 1:Lτo2p1-1
            (m, n) = (i+j-1 <= Lτo2) ? (i+j-1, 1) : (Lτo2p1, i+j-1-Lτo2)
            H[m, n] += H[i, j]
        end
    end
    
    # normalize sums in the first column
    @simd for i in 1:Lτo2p1
        H[i,1] /= i
    end
    
    # normalize the sums in the bottom column
    @simd for j in 2:Lτo2p1
        H[Lτo2p1,j] /= (Lτo2p1-j+1)
    end
    
    # apply the anti-diagonal averages
    @simd for j in 2:Lτo2p1
        for i in 1:Lτo2p1-1
            (m, n) = (i+j-1 <= Lτo2) ? (i+j-1, 1) : (Lτo2p1, i+j-1-Lτo2)
            H[i, j] = H[m, n]
        end
    end

    # apply fixed endpoint
    if isfinite(endpoints[1]) && isfinite(endpoints[2])
        H[1,1] = endpoints[1]
        H[end,end] = endpoints[2]
    end

    return nothing
end


function project_symmetrized!(
    H::AbstractMatrix{T}
) where {T<:AbstractFloat}

    # get Lτ/2 + 1
    Lτo2p1 = size(H, 1)

    # symmetrize around τ = β/2
    if symmetrize
        @simd for i in 1:Lτo2p1
            j = Lτo2p1 - i + 1
            Hsym = (H[i,1] + H[Lτo2p1, j])/2
            H[i,1] = Hsym
            H[Lτo2p1, j] = Hsym
        end
    end

    return nothing
end


@doc raw"""
    project_psd!(
        H::AbstractMatrix{T},
        A::AbstractMatrix{T} = similar(H)
    ) where {T<:AbstractFloat}

Modifying `H` in-place, project it onto the nears positive semidefinite matrix
by setting all the negative eigenvalues to zero.
"""
function project_psd!(
    H::AbstractMatrix{T},
    A::AbstractMatrix{T} = similar(H)
) where {T<:AbstractFloat}

    ϵ, U = eigen(Hermitian(H))
    @. ϵ = max(ϵ, 0.0)
    Ut = adjoint(U)
    D = Diagonal(ϵ)
    mul!(A, D, Ut)
    mul!(H, U, A)

    return nothing
end