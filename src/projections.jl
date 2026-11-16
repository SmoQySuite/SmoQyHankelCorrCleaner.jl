@doc raw"""
    project_positive_curvature!(
        H::AbstractMatrix{T},
        endpoints::Tuple{T,T} = (NaN, NaN)
    ) where {T<:AbstractFloat}

Modifying `H` in-place, project it onto the nearest matrix (in Frobenius norm) whose
anti-diagonal averages ``a_l`` have non-negative second finite difference,
``a_{l+1} - 2a_l + a_{l-1} \ge 0``. If the `endpoints` are finite, the first and last
anti-diagonal averages, i.e. `H[1,1]` and `H[end,end]`, are also fixed to those values.

Only the anti-diagonal averages are constrained, so the projection solves a weighted QP
for the new averages ``f_l``, minimizing ``\sum_l n_l (f_l - a_l)^2`` where
``n_l = \min(l, L-l+1)`` is the length of the ``l``-th anti-diagonal, and then shifts every
entry of anti-diagonal ``l`` by ``f_l - a_l``. See the Algorithm page of the docs.
"""
function project_positive_curvature!(
    H::AbstractMatrix{T},
    endpoints::Tuple{T,T} = (NaN, NaN)
) where {T<:AbstractFloat}

    # get size of Hankel matrix
    N = size(H,1)
    # get number of imaginary time slices
    L = 2*N-1
    # calculate antidiagonal averages
    target = [antidiagonal_average(H,l) for l in 1:L]
    # initialize model
    model = JuMP.Model(HiGHS.Optimizer)
    # setting strings off and silent on for efficiency
    JuMP.set_string_names_on_creation(model, false)
    JuMP.set_silent(model)
    # Cap the number of QP iterations. On some degenerate inputs the HiGHS QP solver cycles
    # without ever certifying optimality and would otherwise run forever. Well-behaved solves
    # need O(L) iterations, so the cap is never reached for them.
    JuMP.set_attribute(model, "qp_iteration_limit", 100 * L)
    # define variables to represent correlations that get fitted
    JuMP.@variable(model, fit[1:L])
    # enforce positive curvature
    JuMP.@constraint(model, [i in 2:L-1], fit[i+1] + fit[i-1] - 2.0 * fit[i] ≥ 0.0)
    # fix endpoints of correlation function
    if isfinite(endpoints[1]) && isfinite(endpoints[2])
        JuMP.@constraint(model, fit[1] == endpoints[1])
        JuMP.@constraint(model, fit[L] == endpoints[2])
    end
    # Minimize the Frobenius distance. Shifting the l'th anti-diagonal, of length n_l, by Δ
    # changes the squared Frobenius norm by n_l Δ², so each term is weighted by n_l.
    JuMP.@objective(model, Min, sum(min(l, L-l+1) * (fit[l] - target[l])^2 for l in 1:L))
    # find the solution
    JuMP.optimize!(model)
    # When the iteration cap is hit while cycling, the solver is already at a feasible point
    # that is (near-)optimal, so accept it. Any other failure is a genuine error.
    status = JuMP.termination_status(model)
    if !(status == JuMP.MOI.OPTIMAL ||
         (status == JuMP.MOI.ITERATION_LIMIT && JuMP.primal_status(model) == JuMP.MOI.FEASIBLE_POINT))
        error("Positive curvature projection failed: termination status $status, primal status $(JuMP.primal_status(model))")
    end
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
        H::AbstractMatrix{T},
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


@doc raw"""
    project_symmetrized!(H::AbstractMatrix{T}) where {T<:AbstractFloat}

Modifying `H` in-place, project it onto the nearest matrix (in Frobenius norm) whose
anti-diagonal averages are symmetric about ``\tau = \beta/2``, i.e. ``a_l = a_{L+1-l}``.
Each anti-diagonal ``l`` is shifted uniformly by ``(a_l + a_{L+1-l})/2 - a_l``.
"""
function project_symmetrized!(
    H::AbstractMatrix{T}
) where {T<:AbstractFloat}

    N = size(H, 1)
    L = 2 * N - 1

    # anti-diagonals l and L+1-l have equal length, so the weighted projection onto
    # a_l = a_{L+1-l} reduces to the plain pairwise average
    a = [antidiagonal_average(H, l) for l in 1:L]
    for j in axes(H, 2), i in axes(H, 1)
        l = i + j - 1
        H[i, j] += (a[L+1-l] - a[l]) / 2
    end

    return nothing
end


@doc raw"""
    project_psd!(
        H::AbstractMatrix{T},
        A::AbstractMatrix{T} = similar(H)
    ) where {T<:AbstractFloat}

Modifying `H` in-place, project it onto the nearest (in Frobenius norm) symmetric positive
semidefinite matrix by symmetrizing `H` and then setting all negative eigenvalues to zero.
The matrix `A` is used as workspace and must not alias `H`.
"""
function project_psd!(
    H::AbstractMatrix{T},
    A::AbstractMatrix{T} = similar(H)
) where {T<:AbstractFloat}

    # symmetrize, since H is not necessarily symmetric during Dykstra iterations
    ϵ, U = eigen(hermitianpart!(H))
    @. ϵ = max(ϵ, 0)
    mul!(A, Diagonal(ϵ), adjoint(U))
    mul!(H, U, A)

    return nothing
end