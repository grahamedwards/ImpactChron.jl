## Interior summation loop in HeveySanders equation
using LoopVectorization

function ΣT(n, r, t, R, λ, κ)
    Σ = zero(typeof(λ))
    @inbounds for nᵢ ∈ n
        Σ += ( ((-1)^nᵢ) / (nᵢ*((nᵢ^2)-( λ*(R^2)/(κ*π^2) ) ) ) ) *
            sin(nᵢ*π*r/R) * exp(-κ*(nᵢ^2)*(π^2)*t/(R^2))
    end
    return Σ
end


function ΣTt(n, r, t, R, λ, κ)
    Σ = zero(typeof(λ))
    @tturbo for nᵢ ∈ n
        a = ifelse(isodd(nᵢ), -1.0, 1.0)
        b = nᵢ*((nᵢ^2)-( λ*(R^2)/(κ*π^2) ) )
        c = sin(nᵢ*π*r/R)
        d = exp(-κ*(nᵢ^2)*(π^2)*t/(R^2))
        Σ += (a / b ) * c * d
    end
    return Σ
end
