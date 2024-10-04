function dichotomy_root(
    f::Function,
    left::T,
    right::T;
    error = eps(T),
    show_iter = false,
    iter_max = 9999,
) where {T<:Real}
    iter = 0
    f_wrap(x) = begin
        iter += 1
        f(x)
    end
    retval(x) = show_iter ? (x, iter) : x
    left_val = f_wrap(left)
    iszero(left_val) && return retval(left)
    right_val = f_wrap(right)
    iszero(right_val) && return retval(right)
    if sign(left_val) == sign(right_val)
        throw(ArgumentError("function xs at borders have same sign"))
    end
    if left_val > right_val
        left, right = right, left
    end
    while abs(right - left) > error && iter < iter_max
        middle = (left + right) / 2
        middle_val = f_wrap(middle)
        if middle_val < 0
            left = middle
            left_val = middle_val
        elseif middle_val > 0
            right = middle
            right_val = middle_val
        else
            return retval(middle)
        end
    end
    return retval((left + right) / 2)
end

function iteration_root(
    f::Function,
    start::Union{T,AbstractArray{T}};
    error = eps(T),
    show_iter = false,
    λ = 1,
    iter_max = 9999,
) where {T<:Number}
    iter = 0
    f_wrap(x) = begin
        iter += 1
        f(x)
    end
    retval(x) = show_iter ? (x, iter) : x
    value = f_wrap(start)
    iszero(value) && return retval(start)
    old_x, x = start, start - λ * value
    while maximum(abs.(x - old_x)) > error && iter < iter_max
        value = f_wrap(x)
        old_x, x = x, x - λ * value
    end
    return retval(x)
end

function newton_root(
    f::Function,
    J::Function,
    start::Union{T,AbstractArray{T}};
    error = eps(T),
    show_iter = false,
    iter_max = 9999,
) where {T<:Number}
    iter = 0
    f_wrap(x) = begin
        iter += 1
        inv(J(x)) * f(x)
    end
    retval(x) = show_iter ? (x, iter) : x
    value = f_wrap(start)
    iszero(value) && return retval(start)
    old_x, x = start, start - value
    while maximum(abs.(x - old_x)) > error && iter < iter_max
        value = f_wrap(x)
        old_x, x = x, x - value
    end
    return retval(x)
end

function secant_root(
    f::Function,
    start::Union{T,AbstractArray{T}};
    error = eps(T),
    show_iter = false,
    λ = 1,
    iter_max = 9999,
) where {T<:Number}
    iter = 0
    f_wrap(x) = begin
        iter += 1
        f(x)
    end
    retval(x) = show_iter ? (x, iter) : x
    value = f_wrap(start)
    iszero(value) && return retval(start)
    old_x, x = start, start - λ * value
    while maximum(abs.(x - old_x)) > error && iter < iter_max
        old_value, value = value, f_wrap(x)
        old_x, x = x, x - (x - old_x) / (value - old_value) * value
    end
    return retval(x)
end
