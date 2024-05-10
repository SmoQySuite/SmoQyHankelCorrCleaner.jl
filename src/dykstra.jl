@doc raw"""
    dykstra!(;
        # KEYWORD ARGUMENTS
        H::AbstractArray{T},
        projections!::Vector{Function},
        tol::T,
        maxiter::Int,
        H0::AbstractArray{T} = zero(H),
        P0::AbstractArray{T} = zero(H),
        Ps::AbstractArray{T} = zeros(eltype(H), size(H)..., length(projections!)),
        verbose::Bool = false,
        projection_names::Vectors{String} = String[]
    ) where {T<:AbstractFloat}

Perform Dykstra's algorithm on an input array `H` given a vector `projections!` for
projection functions that modify `H` in place.
This method is based on the implementation found in the
[Dykstra](https://github.com/mjhough/Dykstra.git) python package.
"""
function dykstra!(;
    # KEYWORD ARGUMENTS
    H::AbstractArray{T},
    projections!::Vector{Function},
    tol::T,
    maxiter::Int,
    H0::AbstractArray{T} = zero(H),
    P0::AbstractArray{T} = zero(H),
    Ps::AbstractArray{T} = zeros(eltype(H), size(H)..., length(projections!)),
    verbose::Bool = false,
    projection_names::Vector{String} = String[]
) where {T<:AbstractFloat}

    # generate default projection names
    if verbose && length(projection_names) < length(projections!)
        projection_names = ["Projection $i" for i in eachindex(projections!)]
    end

    # set previous projections to zero
    fill!(Ps, 0)

    # iteration record
    iter = 0

    # error record
    err = zero(T)

    if verbose
        println("Dykstra's Algorithm Start")
        println()
    end

    # iterate over updates
    for i in 1:maxiter

        if verbose
            println("Iteration = ", i)
            println()
        end

        # initialize squared error to zero
        sqrd_err = zero(T)

        # recrod the iteration
        iter = i

        # iterate over projections
        for p in eachindex(projections!)

            # get the current projection
            projection! = projections![p]

            # get relevant projections matrix
            P = selectdim(Ps, ndims(Ps), p)

            # record H matrix
            copyto!(H0, H)

            # update H by subtracting off previous projection result
            @. H = H - P

            # apply projection
            projection!(H)

            # record P matrix
            copyto!(P0, P)

            # update projection matrix
            @. P = H - (H0 - P0)

            # update current squared error
            @. P0 = P0 - P
            P0_norm = norm(P0)
            sqrd_err += P0_norm^2

            if verbose
                println("\t", projection_names[p], " Error = ", P0_norm)
                println("\tTotal Error = ", sqrt(sqrd_err))
                println()
            end
        end

        # calculate current total error
        err = sqrt(sqrd_err)

        # terminate if error below threshold
        if err < tol 
            break
        end
    end

    if verbose
        println("Dykstra's Algorithm End")
        println()
        println("Final Total Error = ", err)
        println("Number of Iterations = ", iter)
    end

    return (iter, err)
end