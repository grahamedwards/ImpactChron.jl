## Functions used for statistical purposes:
    # histogramify ~ converts data into binned histogram
    # log-likelihood calculators
        # ll_param
        # ll_params
        # ll_dist


"""

```julia
histogramify(domain::AbstractRange,x::AbstractVector,y::AbstractVector)
```

Constructs histogram over (linear) midpoints of `domain` from model outputs in x with corresponding
abundances in y. Requires that each value of x falls within the bounds of `domain`.

Normalizes the output, such that for output `dist` ∑ dist[xᵢ] * Δx = 1
(for each xᵢ in the bincenters of x)

Returns only the histogram masses, bincenters must be calculated externally.
e.g.
```julia
Δd = step(domain)
bincenters= LinRange( first(domain)+Δd/2, last(domain)-Δd/2, length(domain)-1)
```

"""
function histogramify(domain::AbstractRange,x::AbstractVector,y::AbstractVector)
    dist = Vector{float(eltype(y))}(undef,length(domain)-1)
    histogramify!(dist,domain,x,y)
    return dist
end

"""
```julia
histogramify!(dist::AbstractVector,domain::AbstractRange,x::AbstractVector,y::AbstractVector)
```
In-place `histogramify` that ouputs to a pre-allocated vector `dist`

see `histogramify`
"""

function histogramify!(dist::AbstractVector,domain::AbstractRange,x::AbstractVector,y::AbstractVector)
# Declare distribution vector
    fill!(dist,zero(eltype(dist)))
# Calculate Δd
    Δd = step(domain)
# Sort (x,y) values in order of ascending x (e.g. dates) for search efficiency
    i_sorted = sortperm(x)
    x_sort = x[i_sorted]
    y_sort = y[i_sorted]
# Remove any NaNs to ensure
    firstNaN = searchsortedfirst(x_sort,NaN)
    xₛ = view(x_sort,1:firstNaN-1)
    yₛ = view(y_sort,1:firstNaN-1)
# Ensure NaN removal did not delete all elements of x & y
    if length(xₛ)>0
# Identify indices of domain that bound ALL values of x
        xmin = searchsortedfirst(domain,first(xₛ)) - 1
        xmax = searchsortedlast(domain,last(xₛ)) + 1
# Ensure that xmin & xmax are within the domain
    # If first(xₛ) < first(domain), xmin = 0
        xmin = ifelse(iszero(xmin),1,xmin)
    # If last(xₛ) > last(domain), xmax = length(domain) + 1
        xmax = ifelse(xmax>length(domain),length(domain)-1,xmax) # this should never happen. Could remove, but takes <2ns and prevents a segfault if disaster strikes.

        if (xmax - xmin) > 1 # if only 1 bin filled, xmax-xmin=1
            @inbounds for i ∈ xmin:xmax
                l = searchsortedfirst(xₛ , domain[i]) # lower index
                u = searchsortedlast(xₛ , domain[i+1]) # upper index
# Ensure values of (xₛ,yₛ) fall within bounds (if not, searchsortedfirst/searchsortedlast return l > u)
                u >= l && ( dist[i] = vreduce(+,yₛ[l:u]) )
            end
        elseif (xmax - xmin) == 1
            dist[xmin] = 1.0
        end # Note: if each x > OR < domain , xmin >= xmax, and iszero(dist)=true
# Normalize
        ∫distdx = vreduce(+,dist) * Δd
        vmap!(x -> x/∫distdx,dist,dist)
    end
    return dist
end

## Log-likelihood calculators

"""
```julia
ll_param(x::Number,D::T) -> T ∈ {Nrm,lNrm,Unf}

Calculate the log-likelihood that `x` is drawn from a distribution
`D`, where D may be...
Normal (`Nrm`), with mean D.μ and 1σ = D.σ

Lognormal (`lNrm`), with logspace mean D.μ and 1σ = D.σ. Assumes x is already in logspace.

Uniform (`Unf`) with lowerbound D.a and upperbound D.b
(D::Unf always returns loglikelihhood of zero, bounds test is done earlier to speed up code.)

see `ImCh_parameters.jl` for construction of T-structs

"""
ll_param(x::Number,D::Nrm) = -(x-D.μ)*(x-D.μ)/(2*D.σ*D.σ)
ll_param(x::Number,D::lNrm) = -(x-D.μ)*(x-D.μ)/(2*D.σ*D.σ)
    #lnx = log(x); return -lnx-(lnx-D.μ)*(lnx-D.μ) / (2*D.σ*D.σ)
ll_param(x::T,D::Unf) where T<:Number = zero(T)


"""
```julia
ll_params(p::NamedTuple,d::NamedTuple)
```
Calculate log-likelihood for a number of proposals in `p`
with corresponding distributions in `d`

Currently does not evaluate impact (_χ_) parameters.
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
ll_dist(x::AbstractVector,dist::AbstractVector,mu::AbstractVector,sigma::AbstractVector)
```

Calculate loglikelihood that observations in `mu` and `sigma`
are drawn from modeled distribution described by `x` and `dist`
where
`x` contains the bincenters of a normalized histogram `dist` and
`mu` and `sigma` respectively contain the mean and 1σ of the observations.

"""
function ll_dist(   x::AbstractRange, dist::AbstractVector,#p_dist::Tuple{AbstractVector,AbstractVector},    # Proposed distribution (ages,proportions)
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
        iₓ = ceil(Int, (mu[j]-x[1]) / xᵣ * nbtwns + 1)
                # formerly iₓ = searchsortedfirst(x,mu[j]) # x[iₓ] ≥ mu[j]
# If possible, prevent aliasing problems by interpolation
        if (iₓ>1) && (iₓ<=nₚ) && ( (2sigma[j]) < (x[iₓ]-mu[j]) ) && ( (2sigma[j])<(mu[j]-x[iₓ-1]) )
            # && (sigma[j] < (x[iₓ]-x[iₓ-1]) ) # original threshold, see notes.
# Interpolate corresponding distribution value, note: Δx cancels in second term
            likelihood = dist[iₓ] * Δx - (x[iₓ]-mu[j]) * (dist[iₓ]-dist[iₓ-1])
            #likelihood = 6 * sigma[j] * ( dist[iₓ] - (dist[iₓ]-dist[iₓ-1]) * (x[iₓ]-mu[j]) / Δx)
                    #alternate likelihood calculation that only integrates some "width"(e.g. 3σ) of mu distribution...
# Otherwise, sum contributions from Gaussians at each point in distribution
        else
            likelihood = zero(float(eltype(dist)))
            @turbo for i ∈ 1:nₚ     # @turbo faster than @tturbo
# Likelihood curve follows a Gaussian PDF.
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
