## Data Management: functions to aid visualization of ImpactChron outputs.

using Statistics, VectorizedStatistics
using LoopVectorization
using NaNStatistics
import ImpactChron: downscale!, rangemidbounds,rangemidpoints

"""
```julia
interleave(a)
```
Interleave the sequential values of `a` (`<:AbstractVector`) so that each value occurs twice in a `Vector` of `length(a)*2`

Examples
==========

```
julia> interleave([1,2,3])
6-element Vector{Int64}:
 1
 1
 2
 2
 3
 3
```

"""
function interleave(a::AbstractVector)
    isa(a,AbstractRange) && (a=collect(a))
    c = Vector{eltype(a)}(undef, 2*(length(a)))
    i = 0
    for x in 1:length(a)
        c[i += 1] = a[x]
        c[i += 1] = a[x]
    end
    return c
end


"""
```julia
binweave(a)
```
Interleave the sequential values of `a` (`<:AbstractVector`) using the first and last values only once. To be used with `interleave` to make plotting-ready histograms.

See also: `interleave`

Example
==========

```
julia> binweave([1,2,3])
4-element Vector{Int64}:
 1
 2
 2
 3
```
"""
function binweave(a::AbstractVector)
    isa(a,AbstractRange) && (a=collect(a))
    c = Vector{eltype(a)}(undef, 2*(length(a)-1))
    i = 0
    for x in 1:length(a)-1
        c[i += 1] = a[x]
        c[i += 1] = a[x+1]
    end
    return c
end


"""
```julia
distmeans
`
"""
function distmeans(d::Dict;stdev::Bool=false)
    if stdev
        dₒᵤₜ = Dict{Symbol,Tuple{Number,Number}}()
        for i ∈ keys(d)
            if eltype(d[i]) <: Number
                dₒᵤₜ[i]= (vmean(d[i]), vstd(d[i]))
            end
        end
    else
        dₒᵤₜ = Dict{Symbol,Number}()
        for i ∈ keys(d)
            if eltype(d[i]) <: Number
                dₒᵤₜ[i]= vmean(d[i])
            end
        end
    end
    return (; dₒᵤₜ...)
end


function distmedians(d::Dict;ci::Number=0.)
    if ci>0
        α=1-ci
        dₒᵤₜ = Dict{Symbol,Tuple}()
        for i ∈ keys(d)
            if eltype(d[i]) <: Number
                med = median(d[i])
                if length(d[i])>1
                    lower = med - quantile(d[i],α/2)
                    upper = quantile(d[i],1-α/2) - med
                else
                    lower = upper = zero(med)
                end
                dₒᵤₜ[i]= (med,lower,upper)
            end
        end
    else
        dₒᵤₜ = Dict{Symbol,Number}()
        for i ∈ keys(d)
            if eltype(d[i]) <: Number
                dₒᵤₜ[i]= median(d[i])
            end
        end
    end
    return (; dₒᵤₜ...)
end


function cleanhist(x::AbstractArray; nbins::Int=50, scooch_nbins::Int=4)
    x_scooch = (maximum(x)-minimum(x))/ (nbins-scooch_nbins)
    binedges = LinRange(minimum(x)-2*x_scooch,maximum(x)+2*x_scooch,nbins+1)
    bincenters = LinRange( first(binedges)+step(binedges)/2, last(binedges)-step(binedges)/2,nbins)
    y=histcounts(x,binedges)
    return (x=binweave(binedges), y=interleave(y),)
end


"""

"""
function summedpdfhist(x::AbstractRange,d::AbstractVector,ds::AbstractVector; bins::Int=50)
# Ensure downscaling will work correctly
    if !iszero(length(x)%bins)
        x = first(x):step(x): (last(x) - step(x)*length(x)%bins)
        @warn "input range adjusted for binning to x=$x"
    end

# Prepare output range of x
    downscale_factor = div(length(x),bins)
    Δbin = downscale_factor*step(x)
    bincenters = sum(x[1:downscale_factor])/downscale_factor : Δbin : last(x)

# Calculate pdf and downscale into bins
    summedpdf = sumpdfs(x,d,ds)
    binnedpdf = Vector{eltype(d)}(undef,length(bincenters))
    downscale!(binnedpdf,summedpdf)

    return (binweave(rangemidbounds(bincenters)),interleave(binnedpdf))
end

"""

```julia
normdens(x,m,s)
```

Calculate the probability of a normal distribution with mean `m` and standard deviation `s` at a value `x`

"""
function normdens(x::Number,m::Number,s::Number)
	(s*sqrt(2*π))^-1 * exp(-(x-m)*(x-m) / (2*s*s))
end

"""

```julia
sumpdfs(z,x,δx)
```

Sum pdfs of values in `x` with 1σ uncertainties `δx` over a domain defined by `z`, which can be either a `Vector` of constant spacing or a `Range`.

"""
function sumpdfs(z,x,δx)
	densities=Array{eltype(x)}(undef,length(z),length(x));
	
	@tturbo for j in eachindex(x)
		xj = x[j]
		δxj = δx[j]
		for i in eachindex(z)
			densities[i,j] = normdens(z[i],xj,δxj)
		end
	end
# Calculate the stepsize of z for normalization
	isa(z,AbstractRange) ? Δz = step(z) : Δz = z[2]-z[1]
# flatten and normalize z
	vec(vsum(densities,dims=2)) ./ (vsum(densities)*Δz)
	
end

"""

```julia
pdfsample(x::AbstractVector,p::AbstractVector;n::Integer=100)
```

Draw `n` samples from a pdf corresponding to values `x` and probabilities `p`. Calculates a CDF from `p` and linearly interpolates the values of `x` from a `rand()` on the CDF.

"""
function pdfsample(x::AbstractVector,p::AbstractVector;n::Integer=100)
	cp = cumsum(p)
	cp ./= last(cp)
	samples = Vector{eltype(x)}(undef,n)
	for i = 1:n
		r = rand() # draw a random value from the cdf
		ir = searchsortedfirst(cp,r)
# Since searchsortedfirst returns the index ≥ the cumulative probability, interpolate between the returned step and the preceding step where the sample likely lies.
		if ir>1 # If no preceding timestep, just set to age of first timestep
			samples[i] =  x[ir-1] +  (r - cp[ir-1]) * (x[ir]-x[ir-1]) / (cp[ir]-cp[ir-1])
		else			
			samples[i] = x[ir]
		end
	end
	return samples
end

println("""
Data Management Functions Loaded Successfully:

Summary statistics: `distmeans`, `distmedians`
Visualization prep: `interleave`, `binweave`, `cleanhist`, `summedpdfhist`, `normdens`, `sumpdfs`, `pdfsample`
""")
