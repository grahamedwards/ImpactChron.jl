## Functions used for statistical purposes & math support:
    # rangemidpoints 
    # rangemidbounds
    # downscale!
    # ImpactChron.histogramify[!]
    # log-likelihood calculators
        # ll_param
        # ll_params
        # ll_dist
        # ll_dist_params
        # weight_petro_types!


"""

```julia
rangemidpoints(x)
```

Calculate a `LinRange` of the midpoints of each step in `x`(`<:AbstractRange`).

"""
rangemidpoints(x::AbstractRange) = LinRange(first(x) + 0.5step(x), last(x) - 0.5step(x), length(x)-1)

"""

```julia
rangebinbounds(x)
```

Calculate a `LinRange` of the linear bounds for each "midpoint" step in `x` (`<:AbstractRange`).

"""
rangemidbounds(x::AbstractRange) = LinRange(first(x) - 0.5step(x), last(x) + 0.5step(x), length(x)+1)


"""

```julia
ImpactChron.downscale!(B, A)
```

Downscales elements of Array `A` into smaller Array `B` by summing.
Scales the downscaled values in `B` by the ratio of length(A)÷length(B) to preserve any normalizations.

Requires that length(A) % length(B) = 0.

"""
function downscale!(B::AbstractArray{T}, A) where T
    ratio = length(A)÷length(B)
    @assert ratio*length(B) == length(A)

    iA₀, iB₀ = firstindex(A), firstindex(B)
    @inbounds for i ∈ eachindex(B).-1
        Bᵢ = zero(T)
        for n = 0:ratio-1
            Bᵢ += A[iA₀+i*ratio+n]
        end
        B[iB₀+i] = Bᵢ / ratio
    end
    B
end


"""
```julia
ImpactChron.histogramify(domain::AbstractRange,x::AbstractVector,y::AbstractVector)
```

Constructs histogram over bins defined by `domain` from model outputs in `x` with corresponding
abundances in `y`. 
Does not require a constant step in `x`, so this calculates histograms of outputs from `planetesimal_cooling_dates`.

`histogramify` normalizes its output, such that for output `dist` ∑ dist * Δd = 1, so long as all x ∈ domain.
If any x ∉ domain, ∑ dist * Δd = 1- (∑yₒᵤₜ / ∑yₐₗₗ ) where the xₒᵤₜ corresponding to yₒᵤₜ are ∉ domain.

---

Returns only the histogram masses, centers of bins must be calculated separately with `rangemidpoints(domain)`

"""
function histogramify(domain::AbstractRange, x::AbstractVector, y::AbstractVector)
    dist = Vector{float(eltype(y))}(undef,length(domain)-1)
    histogramify!(dist,domain,x,y)
    return dist
end


"""

```julia
ImpactChron.histogramify!(dist::AbstractVector, domain::AbstractRange, x::AbstractVector, y::AbstractVector)
```

In-place `histogramify` that overwites a pre-allocated vector `dist`.
See  [`histogramify`](@ref) for details.

"""
function histogramify!(dist::AbstractVector,domain::AbstractRange,x::AbstractVector,y::AbstractVector)
# start with a fresh zero distribution
    fill!(dist,zero(eltype(dist)))
# Calculate Δd
    Δd = step(domain)
# Sort (x,y) values in order of ascending x (e.g. dates)
    i_sorted = sortperm(x)
    x_sort = x[i_sorted]
    y_sort = y[i_sorted]
# Remove any NaNs to ensure math works
    firstNaN = searchsortedfirst(x_sort,NaN)
    xₛ = view(x_sort,1:firstNaN-1)
    yₛ = view(y_sort,1:firstNaN-1)
# Ensure NaN removal did not delete all elements of x & y
    if length(xₛ)>0
# Calculate area under an unbounded binspace for normalizing relative to all of the (x,y) space
        ∑yΔd = vreduce(+,yₛ) * Δd
# Identify indices of domain that bound ALL values of x
        xmin = searchsortedfirst(domain,first(xₛ)) - 1
        xmax = searchsortedlast(domain,last(xₛ)) + 1
# Ensure that xmin & xmax are within the domain
        xmin = ifelse(iszero(xmin),1,xmin) # # If first(xₛ) < first(domain), xmin = 0
        xmax = ifelse(xmax>length(domain),length(domain)-1,xmax) # If last(xₛ) > last(domain), xmax = length(domain) + 1
            # this should never happen. Could remove, but takes <2ns and prevents a segfault if disaster strikes.

# Ensure any x ∈ domain. # If all x > domain xmin(=len)>xmax(=len-1). If all x < domain xmin=xmax=1. Either case iszero(dist) = true.
        if (xmax - xmin) > 0
            @batch for i ∈ xmin:xmax
                l = searchsortedfirst(xₛ , domain[i]) # lower index
                u = searchsortedlast(xₛ , domain[i+1]) # upper index
# If xₛ does not fall within domain[i:i+1], l > u and for-loop won't run.
                @inbounds for j = l:u
                    dist[i] += yₛ[j]
                end
                dist[i] /= ∑yΔd
            end
        end
    end
    return dist
end

## Log-likelihood calculators

"""

```julia
ll_param(x::Number,D<:PriorDistribution)
```

Calculate the log-likelihood that `x` is drawn from a distribution `D`, which may be a Normal ([`Nrm`](@ref)), Lognormal ([`lNrm`](@ref)), or Uniform ([`Unf`](@ref)) distribution.
For `D::Unf`, `ll_param` returns `0` and bounds checks are performed by [`ImpactChron.prior_bounds`](@ref) within [`ImpactChron.strict_priors`](@ref).

"""
ll_param(x::Number,D::Nrm) = -(x-D.μ)*(x-D.μ)/(2*D.σ*D.σ)
ll_param(x::Number,D::lNrm) = -(x-D.μ)*(x-D.μ)/(2*D.σ*D.σ)
    #lnx = log(x); return -lnx-(lnx-D.μ)*(lnx-D.μ) / (2*D.σ*D.σ)
ll_param(::T,D::Unf) where T<:Number = zero(T)


"""

```julia
ll_params(p::NamedTuple,d::NamedTuple)
```
Calculate log-likelihoods for the proposals in `p`
with corresponding distributions in `d`

Does not evaluate impact (_χ_) parameters, since these were only modeled on Unf distributions.
"""
function ll_params(p::NamedTuple,d::NamedTuple)
    ll = 0.
    ll += ll_param(p.Cp, d.Cp)
    ll += ll_param(p.R, d.R)
    ll += ll_param(p.Tc, d.Tc)
    ll += ll_param(p.Tm, d.Tm)
    ll += ll_param(p.cAl, d.cAl)
    ll += ll_param(p.k, d.k)
    ll += ll_param(p.rAlo, d.rAlo)
    ll += ll_param(p.ta, d.ta)
    ll += ll_param(p.tss, d.tss)
    ll += ll_param(p.ρ, d.ρ)
    return ll
end

"""

```julia
ll_dist(x, dist, mu, sigma)
```

Calculate loglikelihood that observations in Vectors `mu` and `sigma`
are drawn from modeled distribution described by Vectors `x` and `dist`,
where `x` contains the bincenters of a normalized histogram `dist` and
`mu` and `sigma` respectively contain the mean and 1σ of the observations.

"""
function ll_dist(   x::AbstractRange, dist::AbstractVector,
                    mu::AbstractVector,      # 1D Array/Vector of observed μ's (sorted)
                    sigma::AbstractVector)   # 1D Array/Vector of observed σ's (sorted)

# Define some frequently used variables
    ll = zero(float(eltype(dist)))
    nₘ = length(mu)            # n of "Measured" data
    nₚ = length(x)     # n of bincenters
    nbtwns = nₚ - 1   # n of spaces between bin centers.
    Δx = abs(step(x))
    xᵣ = abs(last(x) - first(x)) # range of x
    #tkr = zeros(nₘ)

# Cycle through each datum in (mu,sigma)
    @inbounds for j ∈ 1:nₘ 
# Find index of μ in the `dist` array
        iₓ = ceil(Int, (mu[j]-x[1]) / xᵣ * nbtwns + 1) # x[iₓ] ≥ mu[j], equivalent to searchsortedfirst(x,mu[j])
# If possible, prevent aliasing problems by interpolation
        if (1 < iₓ<=nₚ) && ( (2sigma[j]) < (x[iₓ]-mu[j]) ) && ( (2sigma[j])<(mu[j]-x[iₓ-1]) )
            # && (sigma[j] < (x[iₓ]-x[iₓ-1]) ) # original threshold, see notes.
# Interpolate corresponding distribution value, note: Δx cancels in second term
            likelihood = dist[iₓ] * Δx - (x[iₓ]-mu[j]) * (dist[iₓ]-dist[iₓ-1])
# Otherwise, sum contributions from Gaussians at each point in distribution
        else
            likelihood = zero(float(eltype(dist)))
            @turbo for i ∈ 1:nₚ  
# Likelihood curve follows a Gaussian PDF about mu[j]
                likelihood += ( dist[i] / (sigma[j] * sqrt(2*π)) ) *
                        exp( - (x[i]-mu[j])*(x[i]-mu[j]) / (2*sigma[j]*sigma[j]) )
            end
            likelihood*=Δx
        end
        ll += log(ifelse(iszero(likelihood),eps(),likelihood))
# Normalize by total area under curve for intercomparability among proposals.
    end
    return ll
end



"""

```julia
ll_dist_params(a, p, plims, mu,sig)
```

Calculate the combined log-likelihood that 

1. Observations with mean `mu` and corresponding 1σ uncertainties `sig` (both ::`Array`) are drawn from the downscaled age distribution and timesteps contained in `a` (::`AsteroidHistory`)

2. The proposal parameters `p` used to calculate this history are drawn from the prior distributions `plims` (both `::NamedTuple`)

If the age distribution is zero, quickly returns a log-liklihood of `-Inf`.

"""
function ll_dist_params(a::AsteroidHistory, p::NamedTuple, plims::NamedTuple, mu::AbstractArray,sig::AbstractArray)
# Calculate log-likelihood, skipping if already `-Inf`
    if iszero(a.agedist_downscaled)
        ll =-Inf
    else
        ll = ll_dist(a.t_downscaled, a.agedist_downscaled, mu, sig) + ll_params(p,plims)
    end
    ll
end


"""

```julia
ImpactChron.weight_petro_types!(v, T, petrotypes::PetroTypes)
```

Reweight volumetric fractions relative to the abundance of each petrologic type in the meteorite record.
Takes `Array`s of volumetric fraction (`v`) and peak temperature (`T`) as output by [`planetesimal_cooling_dates`](@ref) or contained within an [`AsteroidHistory`](@ref).

Requires all petrologic types to occupy at least one radial node, otherwise returns `zero` in all `v`.

see also: [`PetroTypes`](@ref)

"""
function weight_petro_types!(v::AbstractArray,T::AbstractArray,petrotypes::PetroTypes)
# Find first non-melted index (volume fraction > 0)
   imelt = findfirst(!iszero,v) 

# Find the indices where temperatures are < max temp of petrologic type (except type 6)
    i3 = findfirst(x -> x<=petrotypes.type3.T, T)
    i4 = findfirst(x -> x<= petrotypes.type4.T, T)
    i5 = findfirst(x -> x<=petrotypes.type5.T, T)

# Require all petrologic types to occur and occupy at least one layer. 
    if !isnothing(i3) && lastindex(T) >= i3 > i4 > i5 > imelt
# Calculate and apply the conversion between fractional volume of body and fractional abundance in the meteorite record.
      pv3 = petrotypes.type3.p / vsum(@view v[i3:lastindex(T)])
      pv4 = petrotypes.type4.p / vsum(@view v[i4:i3-1])
      pv5 = petrotypes.type5.p / vsum(@view v[i5:i4-1])
      pv6 = petrotypes.type6.p / vsum(@view v[imelt:i5-1])

      @inbounds for i = i3:lastindex(T)
         v[i] *= pv3
      end
      @inbounds for i = i4:i3-1
         v[i] *= pv4
      end
      @inbounds for i = i5:i4-1
         v[i] *= pv5
      end
      @inbounds for i = firstindex(T):i5-1
         v[i] *= pv6
      end
# If all petrologic types do not occupy ≥1 layer, reject all dates. This thing's weird.
   else
      fill!(v,zero(eltype(v)))
      printstyled("missing petrologic types: "; color=:light_magenta);flush(stdout)
   end
   v
end
