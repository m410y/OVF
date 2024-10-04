function identity_matrix(dim::Integer; type = Float64)
    return dim > 1 ? one(Matrix{type}(undef, dim, dim)) : 1
end

function thomas_algorithm(
    bot::AbstractVector{T},
    diag::AbstractVector{T},
    top::AbstractVector{T},
    vals::AbstractVector{T},
) where {T<:Number}
    N = length(diag)
    length(bot) == N - 1 || throw(ArgumentError("wrong bot size"))
    length(top) == N - 1 || throw(ArgumentError("wrong top size"))
    length(vals) == N || throw(ArgumentError("wrong vals size"))
    vals_top = [vals [top; 0.0]]
    vals_top[begin, :] /= diag[begin]
    for i = 1:(N-1)
        vals_top[i+1, 1] -= bot[i] * vals_top[i, 1]
        vals_top[i+1, :] /= diag[i+1] - bot[i] * vals_top[i, 2]
    end
    for i = (N-1):-1:1
        vals_top[i, 1] -= vals_top[i, 2] * vals_top[i+1, 1]
    end
    return vals_top[:, 1]
end

function thomas_algorithm!(
    bot::AbstractVector{T},
    diag::AbstractVector{T},
    top::AbstractVector{T},
    vals::AbstractVector{T},
) where {T<:Number}
    N = length(diag)
    length(bot) == N - 1 || throw(ArgumentError("wrong bot size"))
    length(top) == N - 1 || throw(ArgumentError("wrong top size"))
    length(vals) == N || throw(ArgumentError("wrong vals size"))
    for i = 1:(N-1)
        diag[i+1] -= bot[i] / diag[i] * top[i]
        vals[i+1] -= bot[i] / diag[i] * vals[i]
    end
    vals[end] /= diag[end]
    for i = (N-1):-1:1
        vals[i] -= top[i] * vals[i+1]
        vals[i] /= diag[i]
    end
    return vals
end
