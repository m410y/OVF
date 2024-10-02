function identity_matrix(dim::Integer; type=Float64)
    return dim > 1 ? one(zeros(type, dim, dim)) : 1
end

function thomas_algorithm(bot::V, diag::V, top::V, vals::V) where V <: AbstractVector{T} where T <: Number
    (   length(vals) != length(diag) ||
        length(top) != length(diag) - 1 ||
        length(bot) != length(diag) - 1
    ) && throw(ArgumentError("wrong sizes"))
    vals_top = [vals [top; 0.0]]
    vals_top[begin, :] /= diag[begin]
    for i in 1:(length(diag) - 1)
        vals_top[i+1, 1] -= bot[i] * vals_top[i, 1]
        vals_top[i+1, :] /= diag[i+1] - bot[i] * vals_top[i, 2]
    end
    for i in (length(diag) - 1):-1:1
        vals_top[i, 1] -= vals_top[i, 2] * vals_top[i+1, 1]
    end
    return vals_top[:, 1]
end

function thomas_algorithm!(bot::V, diag::V, top::V, vals::V) where V <: AbstractVector{T} where T <: Number
    (   length(vals) != length(diag) ||
        length(top) != length(diag) - 1 ||
        length(bot) != length(diag) - 1
    ) && throw(ArgumentError("wrong sizes"))
    for i in 1:(length(diag) - 1)
        diag[i+1] -= bot[i] / diag[i] * top[i]
        vals[i+1] -= bot[i] / diag[i] * vals[i]
    end
    vals[end] /= diag[end]
    for i in (length(diag) - 1):-1:1
        vals[i] -= top[i] * vals[i+1]
        vals[i] /= diag[i]
    end
    return vals
end