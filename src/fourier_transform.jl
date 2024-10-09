function DFT(vals::AbstractArray{T}) where {T<:Number}
    N = length(vals)
    idxs = 0:N-1
    F_p = 
    return exp.(-1im * 2pi / N * idxs * idxs') * vals
end
