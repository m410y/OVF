function identity_matrix(dim::Integer; type = Float64)
    return dim > 1 ? one(Matrix{type}(undef, dim, dim)) : 1
end

function norm(vector::AbstractVector{T}) where {T<:Number}
    return sqrt(vector'vector)
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

function power_method(
    A::AbstractMatrix{T};
    start = rand(T, first(size(A))),
    error = eps(T),
    iter_max = 9999,
) where {T<:Real}
    N = length(start)
    size(A) == (N, N) || throw(ArgumentError("wrong A size"))
    X = start / norm(start)
    λ_old, λ = Inf, 1
    iter = 0
    while abs(λ_old - λ) > error && iter < iter_max
        X = A * X
        λ_old, λ = λ, norm(X)
        X /= λ
    end
    return λ, X
end

function power_tridiag_method(
    bot::AbstractVector{T},
    diag::AbstractVector{T},
    top::AbstractVector{T};
    start = rand(T, first(size(diag))),
    error = eps(T),
    iter_max = 9999,
) where {T<:Real}
    N = length(start)
    length(bot) == N - 1 || throw(ArgumentError("wrong bot size"))
    length(top) == N - 1 || throw(ArgumentError("wrong top size"))
    length(diag) == N || throw(ArgumentError("wrong diag size"))
    X = start / norm(start)
    λ_old, λ = Inf, 1
    iter = 0
    while abs(λ_old - λ) > error && iter < iter_max
        X_new = diag .* X
        X_new[begin+1:end] += bot .* X[begin:end-1]
        X_new[begin:end-1] += top .* X[begin+1:end]
        λ_old, λ = λ, norm(X_new)
        X = X_new / λ
    end
    return λ, X
end

function inverse_tridiag_method(
    bot::AbstractVector{T},
    diag::AbstractVector{T},
    top::AbstractVector{T};
    start = rand(T, first(size(diag))),
    error = eps(T),
    iter_max = 9999,
) where {T<:Real}
    N = length(start)
    length(bot) == N - 1 || throw(ArgumentError("wrong bot size"))
    length(top) == N - 1 || throw(ArgumentError("wrong top size"))
    length(diag) == N || throw(ArgumentError("wrong diag size"))
    X = start / norm(start)
    λ_old, λ = Inf, 1
    iter = 0
    while abs(λ_old - λ) > error && iter < iter_max
        X = thomas_algorithm(bot, diag, top, X)
        λ_old, λ = λ, 1 / norm(X)
        X *= λ
    end
    return λ, X
end
