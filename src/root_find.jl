function dichotomy_root(f::Function, left::T, right::T; error=eps(T), verbose=false, max_iter=9999) where T <: Real
    counter = 0
    f_wrap(x) = begin
        counter += 1
        f(x)
    end
    retval(val) = verbose ? (val, counter) : val
    left_value = f_wrap(left)
    iszero(left_value) && return retval(left)
    right_value = f_wrap(right)
    iszero(right_value) && return retval(right)
    if sign(left_value) == sign(right_value)
        throw(ArgumentError("function xs at borders have same sign"))
    end
    if left_value > right_value
        left, right = right, left
    end
    while abs(right - left) > error
        middle = (left + right) / 2
        middle_value = f_wrap(middle)
        if middle_value < 0
            left = middle
            left_value = middle_value
        elseif  middle_value > 0
            right = middle
            right_value = middle_value
        else
            return retval(middle)
        end
        if counter > max_iter
            throw(OverflowError("too much iterations (> $(max_iter))"))
        end
    end
    return retval((left + right) / 2)
end

function iteration_root(f::Function, start::T; error=eps(T), λ=1.0, verbose=false, max_iter=9999) where T <: Union{Number, AbstractArray}
    counter = 0
    f_wrap(x) = begin
        counter += 1
        f(x)
    end
    retval(x) = verbose ? (x, counter) : x
    value = f_wrap(start)
    iszero(value) && return retval(start)
    old_x, x = start, start - λ * value
    while abs(x - old_x) > error
        value = f_wrap(x)
        old_x, x = x, x - λ * value
        if counter > max_iter
            throw(OverflowError("too much iterations (> $(max_iter))"))
        end
    end
    return retval(x)
end

function newton_root(f::Function, J::Function, start::T; error=eps(T), verbose=false, max_iter=9999) where T <: Union{Number, AbstractArray}
    counter = 0
    f_wrap(x) = begin
        counter += 1
        inv(J(x)) * f(x)
    end
    retval(x) = verbose ? (x, counter) : x
    value = f_wrap(start)
    iszero(value) && return retval(start)
    old_x, x = start, start - value
    while abs(x - old_x) > error
        value = f_wrap(x)
        old_x, x = x, x - value
        if counter > max_iter
            throw(OverflowError("too much iterations (> $(max_iter))"))
        end
    end
    return retval(x)
end

function secant_root(f::Function, start::T; error=eps(T), verbose=false, λ=1.0, max_iter=9999) where T <: Union{Number, AbstractArray}
    counter = 0
    f_wrap(x) = begin
        counter += 1
        f(x)
    end
    retval(x) = verbose ? (x, counter) : x
    value = f_wrap(start)
    iszero(value) && return retval(start)
    old_x, x = start, start - λ * value
    while abs(x - old_x) > error
        old_value, value = value, f_wrap(x)
        old_x, x = x, x - (x - old_x) / (value - old_value) * value
        if counter > max_iter
            throw(OverflowError("too much iterations (> $(max_iter))"))
        end
    end
    return retval(x)
end