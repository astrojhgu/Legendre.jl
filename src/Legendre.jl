module Legendre

import SpecialFunctions
import Polynomials

CONV_MAT=[]
MAX_ORD_CACHED=10
function __init__()
    global CONV_MAT
    CONV_MAT=leg2poly_matrix(Float64, MAX_ORD_CACHED)
    println("module initialized")
end


function log_factorial(x::T)::T where {T<:AbstractFloat}
    SpecialFunctions.lgamma(x+1)
end

function coefficient(n::T, k::T, R::DataType)::R where {T<:Integer}
    if k>n
        zero(R)
    else
        two=one(T)+one(T)
        m=div((n-k), two)
        if m*two!=(n-k)
            zero(R)
        else
            if (m%two==0)  one(R) else -one(R) end * exp(log_factorial(R(two*(n-m)))-log_factorial(R(m))-log_factorial(R(n-m))-log_factorial(R(n-two*m))-R(n)*log(two))
        end
    end
end

function leg2poly_coeff(p::AbstractArray{Float64, 1})::AbstractArray{Float64,1}
    n=length(p)
    order=n-1
    #if order
    if order<=MAX_ORD_CACHED
        return CONV_MAT[1:order+1, 1:order+1]*p
    end
    result=zeros(Float64, n)
    for i in 1:n
        for j in 0:i-1
            result[j+1]+=p[i]*coefficient(i-1, j, T)
        end
    end
    result
end


function leg2poly_coeff(p::AbstractArray{T, 1})::AbstractArray{T,1} where {T<:AbstractFloat}
    n=length(p)

    result=zeros(T, n)
    for i in 1:n
        for j in 0:i-1
            result[j+1]+=p[i]*coefficient(i-1, j, T)
        end
    end
    result
end

function leg2poly_matrix(T::DataType, order::Int)::Array{T,2}
    result=zeros(T, order+1, order+1)
    for i in 1:order+1
        for j in 0:order
            result[j+1,i]=coefficient(i-1, j, T)
        end
    end
    result
end

function leg2poly(p::AbstractArray{T, 1})::Polynomials.Poly{T} where {T<:AbstractFloat}
    Polynomials.Poly(leg2poly_coeff(p))
end

end # module
