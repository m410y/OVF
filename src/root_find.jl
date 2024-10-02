function dichotomy_root(f::Function, left::Real, right::Real; error=1e-5, verbose=false, max_iter=9999)
    counter = 0
    f_wrap(val) = begin
        counter += 1
        f(val)
    end
    retval(val) = verbose ? (val, counter) : val
    left_value = f_wrap(left)
    left_value == 0 && return retval(left)
    right_value = f_wrap(right)
    right_value == 0 && return retval(right)
    if sign(left_value) == sign(right_value)
        throw(ArgumentError("function xs at borders have same sign"))
    end
    if left_value > right_value
        left, right = right, left
    end
    while abs(right_value - left_value) > error
        middle = (left + right) / 2
        middle_x = f_wrap(middle)
        if middle_x < 0
            left = middle
            left_value = middle_x
        elseif  middle_x > 0
            right = middle
            right_value = middle_x
        else
            return retval(middle)
        end
        if counter > max_iter
            throw(OverflowError("too much iterations (> $(max_iter))"))
        end
    end
    return retval((left + right) / 2)
end

function iteration_root(f::Function, start::Real; error=1e-5, λ=1.0, verbose=false, max_iter=9999)
    counter = 0
    f_wrap(val) = begin
        counter += 1
        val - λ*f(val)
    end
    retval(val) = verbose ? (val, counter) : val
    x = f_wrap(start)
    x == start && return retval(start)
    old_x, x = x, f_wrap(x)
    while abs(x - old_x) > error
        old_x, x = x, f_wrap(x)
        if counter > max_iter
            throw(OverflowError("too much iterations (> $(max_iter))"))
        end
    end
    return retval(x)
end

function newton_root(f::Function, J::Function, start::Real; error=1e-5, verbose=false, max_iter=9999)
    counter = 0
    f_wrap(val) = begin
        counter += 1
        val - inv(J(val))*f(val)
    end
    retval(val) = verbose ? (val, counter) : val
    x = f_wrap(start)
    x == start && return retval(start)
    old_x, x = x, f_wrap(x)
    while abs(x - old_x) > error
        old_x, x = x, f_wrap(x)
        if counter > max_iter
            throw(OverflowError("too much iterations (> $(max_iter))"))
        end
    end
    return retval(x)
end