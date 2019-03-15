module Legendre

import SpecialFunctions
import Polynomials

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

function leg2poly(p::AbstractArray{T, 1})::Polynomials.Poly{T} where {T<:AbstractFloat}
    Polynomials.Poly(leg2poly_coeff(p))
end

end # module
