function DFT(vals::AbstractArray{T}) where {T<:Number}
    N = length(vals)
    idxs = 0:N-1
    F = exp.(-1im * 2pi / N * idxs * idxs')
    fourier = F * vals
    return fourier
end

function DFT_inverse(fourier::AbstractArray{T}) where {T<:Number}
    N = length(fourier)
    idxs = 0:N-1
    F = exp.(1im * 2pi / N * idxs * idxs') / N
    vals = F * fourier
    return vals
end

function FFT(vals::AbstractArray{T}) where {T<:Number}
    ispow2(length(vals)) || throw(ArgumentError("length is not power of 2"))

end
