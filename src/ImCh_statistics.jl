## Functions used for statistical purposes & math support:
    # turbosum & tturbosum
    # rangemidpoints
    # histogramify ~ converts data into binned histogram
    # log-likelihood calculators
        # ll_param
        # ll_params
        # ll_dist

"""

```julia
turbosum(x::AbstractArray)
```
Fast summing of x with the power of LoopVectorization.jl's @turbo.
Since `vreduce` does not accept `view`s, we `turbosum`!

Faster than reduce for length(x)>≈200
"""

function turbosum(x::AbstractArray)
    ∑x = zero(eltype(x))
    @turbo for i in eachindex(x)
        ∑x += x[i]
    end
    return ∑x
end

"""

```julia
tturbosum(x::AbstractArray)
```
Fast summing of x with the power of LoopVectorization.jl's @tturbo (multithreaded turbo).
"""
function tturbosum(x::AbstractArray)
    ∑x = zero(eltype(x))
    @tturbo for i in eachindex(x)
        ∑x += x[i]
    end
    return ∑x
end

"""

```julia
rangemidpoints(x::AbstractRange)
```
Calculate a `LinRange` of the midpoints of each step in `x`(`<:AbstractRange`).
"""
rangemidpoints(x::AbstractRange) = LinRange(first(x) + 0.5step(x), last(x) - 0.5step(x), length(x)-1)

"""
```julia
histogramify(domain::AbstractRange,x::AbstractVector,y::AbstractVector)
```

Constructs histogram over (linear) midpoints of `domain` from model outputs in x with corresponding
abundances in y. Does not require a constant step in `x`.

Normalizes the output, such that for output `dist` ∑ dist[dᵢ] * Δd = 1
(for each dᵢ in the bincenters of domain with step-size Δd), so long as all x ∈ domain.
If any x ∉ domain, ∑ dist[dᵢ] * Δd = 1- (∑yₒᵤₜ / ∑yₐₗₗ ) where the corresponding xₒᵤₜ of each yₒᵤₜ is ∉ domain.

***
***
Returns only the histogram masses, centers of time bins must be calculated externally.
e.g.
```julia
Δd = step(domain)
bincenters= LinRange( first(domain)+Δd/2, last(domain)-Δd/2, length(domain)-1)
```
or see `rangemidpoints`

"""
function histogramify(domain::AbstractRange, x::AbstractVector, y::AbstractVector)
    dist = Vector{float(eltype(y))}(undef,length(domain)-1)
    histogramify!(dist,domain,x,y)
    return dist
end


"""

```julia
histogramify!(dist::AbstractVector, domain::AbstractRange, x::AbstractVector, y::AbstractVector)

```
In-place `histogramify` that overwites a pre-allocated vector `dist`.

see `histogramify` for details
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



"""

```julia
ll_dist(x::AbstractVector,dist::AbstractVector,mu::AbstractVector,sigma::AbstractVector)
```

Calculate loglikelihood that observations in `mu` and `sigma`
are drawn from modeled distribution described by `x` and `dist`
where
`x` contains the bincenters of a normalized histogram `dist` and
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
    Δx = step(x)
    xᵣ = abs(last(x) - first(x)) # range of x
    #tkr = zeros(nₘ)

# Cycle through each datum in (mu,sigma)
    @inbounds for j ∈ 1:nₘ #might batch speed this up?
# Find index of μ in the `dist` array
        iₓ = ceil(Int, (mu[j]-x[1]) / xᵣ * nbtwns + 1) # x[iₓ] ≥ mu[j], equivalent to searchsortedfirst(x,mu[j])
# If possible, prevent aliasing problems by interpolation
        if (iₓ>1) && (iₓ<=nₚ) && ( (2sigma[j]) < (x[iₓ]-mu[j]) ) && ( (2sigma[j])<(mu[j]-x[iₓ-1]) )
            # && (sigma[j] < (x[iₓ]-x[iₓ-1]) ) # original threshold, see notes.
# Interpolate corresponding distribution value, note: Δx cancels in second term
            likelihood = dist[iₓ] * Δx - (x[iₓ]-mu[j]) * (dist[iₓ]-dist[iₓ-1])
# Otherwise, sum contributions from Gaussians at each point in distribution
        else
            likelihood = zero(float(eltype(dist)))
            @turbo for i ∈ 1:nₚ     # @turbo faster than @tturbo
# Likelihood curve follows a Gaussian PDF about mu[j]
                likelihood += ( dist[i] / (sigma[j] * sqrt(2*π)) ) *
                        exp( - (x[i]-mu[j])*(x[i]-mu[j]) / (2*sigma[j]*sigma[j]) )
            end
            likelihood*=Δx
        end
        ll += log(likelihood)
# Normalize by total area under curve for intercomparability among proposals.
    end
    return ll
end
