## Interior summation loop in HeveySanders equation
using LoopVectorization

function ΣT(n::UnitRange, r, t, R, λ, κ)
    Σ = zero(typeof(λ))
    @inbounds for nᵢ ∈ n
        Σ += ( ((-1)^nᵢ) / (nᵢ*((nᵢ^2)-( λ*(R^2)/(κ*π^2) ) ) ) ) *
            sin(nᵢ*π*r/R) * exp(-κ*(nᵢ^2)*(π^2)*t/(R^2))
    end
    return Σ
end


function ΣTt(n::UnitRange, r, t, R, λ, κ)
    Σ = zero(typeof(λ))
    @tturbo for nᵢ ∈ n
        a = ifelse(isodd(nᵢ), -1.0, 1.0)
        b = nᵢ*((nᵢ^2)-( λ*(R^2)/(κ*π^2) ) )
        c = sin(nᵢ*π*r/R)   #
        d = exp(-κ*(nᵢ^2)*(π^2)*t/(R^2))
        Σ += (a / b ) * c * d
    end
    return Σ
end


# Compare to original math to confirm same output.
#n=1:1000
sum( @. ( ((-1)^n) / (n*((n^2)-( λ*(R^2)/(κ*π^2) ) ) ) ) *
    sin(n*π*r/R) * exp(-κ*(n^2)*(π^2)*t/(R^2)) )
ΣTt(n,r,t,R,λ,κ)
# EQUIVALENT essentially to eps().


function ΣTtp(n::UnitRange, r, t, R, λ, κ)
    Σ = zero(typeof(λ))
    #κπ2 = κ * π * π
    #R2 = R * R
    @tturbo for nᵢ ∈ n
        a = ifelse(isodd(nᵢ), -1.0, 1.0)
        b = nᵢ*((nᵢ^2)-( λ*(R*R)/(κ*π*π) ) )
        c = sin(nᵢ*π*r/R)   #
        d = exp(-(nᵢ^2)*(κ*π*π)*t/(R*R))
        Σ += (a / b ) * c * d
    end
    return Σ
end

## Interior summation loop in HeveySanders equation
# Part II:
"""
How close to Inf does n have to be?
n=1:300 is plenty.
"""

# Fixed values for a given body...
R = 150e3
λ =  log(2)/7.17e5/(365*24*60*60)
κ =  4. / (3210. * 900.)

# Values that change
t = 0
r = 1e3


function Σntest(n::AbstractArray,r::AbstractArray,t::AbstractArray,;R=150e3, λ =  log(2)/7.17e5/(365*24*60*60), κ = 4. / (3210. * 900.))
    t *= (365*24*60*60*1e6)

    nₙ = length(n)
    nₜ = length(t)
    nᵣ = length(r)

    Σout = Array{Float64}(undef, nₙ , nₜ , nᵣ)
    for i ∈ 1:nₙ
        for j ∈ 1:nₜ
            for k ∈ 1:nᵣ
                Σout[i,j,k] = ΣTt(1:n[i],r[k],t[j],R,λ,κ)
            end
        end
    end

    dd = Tuple(findall(x -> x == 1, size(Σout)))
    return dropdims(Σout,dims=dd)
end


plot(n,Σntest(n,[0],[1e3]),xaxis=:log)


# Time dependence:
plot(0:.1:100,Σntest([100000],[1e3,100e3],0:.1:100), xlabel="Time (Ma after CAIs)",labels=["r = 1 km" "r = 100 km"])
    # Absolute magnitude decays with time

# radial distance dependence
plot(1:150,Σntest([100000],1e3:1e3:150e3,[0,1,5,10,100])', xlabel="radial distance (km)",labels=["0 Ma" "1 Ma" "5 Ma" "10 Ma" "100 Ma"])
    # Varies sinusoidally across the radius.
    # Dampens with time.

# Summation limit (n)
nplot1=plot(1:1000, Σntest(1:1000,[50e3],[0]),xaxis=:log,xlabel="n",labels="")
    nplot2=plot(10:10000, Σntest(10:10000,[50e3],[0]),xaxis=:log,xlabel="n",labels="")
    nplot3=plot(100:10000, Σntest(100:10000,[50e3],[0]),xaxis=:log,xlabel="n",labels="")
    plot(nplot1,nplot2,nplot3,layout=(3,1))
#  100 looks sufficient at the per mille level.

n100= ΣTt(1:100, 50e3, 0, 150e3, log(2)/7.17e5/(365*24*60*60), 4. / (3210. * 900.))
n101= ΣTt(1:101, 50e3, 0, 150e3, log(2)/7.17e5/(365*24*60*60), 4. / (3210. * 900.))

n300= ΣTt(1:300, 50e3, 0, 150e3, log(2)/7.17e5/(365*24*60*60), 4. / (3210. * 900.))
n500= ΣTt(1:500, 50e3, 0, 150e3, log(2)/7.17e5/(365*24*60*60), 4. / (3210. * 900.))
n1000= ΣTt(1:1000, 50e3, 0, 150e3, log(2)/7.17e5/(365*24*60*60), 4. / (3210. * 900.))
n1M= ΣTt(1:1000, 50e3, 0, 150e3, log(2)/7.17e5/(365*24*60*60), 4. / (3210. * 900.))

Δn100 = abs(n101-n100) / ((n101+n100)/2)
Δn1000 = abs(n1000-n1M) / n1M
Δn500 = abs(n500-n1M) / n1M
Δn300 = abs(n300-n1M) / n1M

## Testing speed of midpoint calculations
    # choice works for LinRange style vectors

dot_add(a,b,sh) = (0.5 * a / b) .+ sh[1:end-1] # > order of magnitude faster
diff_add(sh) = 0.5 .* diff(sh) .+ sh[1:end-1]
