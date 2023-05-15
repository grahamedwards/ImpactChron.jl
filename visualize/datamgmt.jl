## Data Management: functions to aid visualization of ImpactChron outputs.

using VectorizedStatistics # fast summary statistics, e.g. vmean, vmedian, vquantile, vstd
using LoopVectorization, Polyester # We need speed.
using ImpactChron # This is probably already loaded, but just in case.
using NaNStatistics



const ParamTitles = Dict(   
    :tss => "Solar age (CAIs, Ma)", 
    :rAlo=>"Initial ²⁶Al/²⁷Al (×10⁻⁵)",
    :Tm=> "Disk temperature (K, 2.5 AU)",
    :R => "Radius (km)",
    :ta=> "Accretion time (Myₛₛ)", 
    :cAl=> "Al abundance (wt%)",
    :ρ=> "Bulk density (kg/m³)",
    :Cp=> "Specific heat Capacity (J/kg•K)",
    :k => "Thermal conductivity (W/m•K)",
    :Tc=> "Ar closure temperature (K)",
    :tχα=> "α bombardment onset (Myₛₛ)",
    :τχα=> "α bombardment ℯ-folding time (My)",
    :Fχα=> "α initial impactor flux (My⁻¹)",
    :tχβ=> "β bombardment onset (Myₛₛ)",
    :τχβ=> "β bombardment ℯ-folding time (My)",
    :Fχβ=> "β initial impactor flux (My⁻¹)",
    :tχγ=> "γ bombardment onset (My after CAIs)",
    :τχγ=> "γ bombardment ℯ-folding time (My)",
    :Fχγ=> "γ initial impactor flux (My⁻¹)",            )

const ParamVars = Dict( :tss => "tₛₛ", :rAlo=>"²⁶Alₒ", :Tm=> "Tₘ", :R => "R", :ta=> "tₐ", :cAl=> "[Al]", :ρ=> "ρ", :Cp=> "Cₚ", :k => "K", :Tc=> "Tc",:tχα=> "tₒα", :τχα=> "τα", :Fχα=> "Fₒα", :tχβ=> "tₒβ", :τχβ=> "τβ", :Fχβ=> "Fₒβ", :tχγ=> "tₒγ", :τχγ=> "τγ", :Fχγ=> "Fₒγ")


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
    c = Vector{eltype(a)}(undef, 2*(length(a)))
    interleave!(c,a)
    return c
end

"""
```julia
interleave!(c,a)
```
In-place version of `interleave`, that fills `c` with pairs of each index in `a`.
"""
function interleave!(c::AbstractVector,a::AbstractVector)
    @assert length(c) == 2*length(a)
    i = 0
    @inbounds for x in 1:length(a)
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
    c = Vector{eltype(a)}(undef, 2*(length(a)-1))
    binweave!(c,a)
    return c
end


"""
```julia
binweave!(c,a)
```
In-place version of `binweave`, that overwrites `c` with woven bin edges.

See also: `binweave`
"""
function binweave!(c::AbstractVector, a::AbstractVector)   
    @assert length(c) == 2*(length(a)-1)
    i = 0
    @inbounds for x in 1:length(a)-1
        c[i += 1] = a[x]
        c[i += 1] = a[x+1]
    end
    c 
end

"""
```julia
distmeans(d::Dict;stdev=false)
```

Calculate means of each array in a `Dict` of `Array{<:Number}`.

Returns a `NamedTuple` of each key and its mean (default) or a `Tuple` of `(μ,σ)` if `stdev=true`

"""
function distmeans(d::Dict;stdev::Bool=false)
    if stdev
        dₒᵤₜ = Dict{Symbol,Tuple{Number,Number}}()
        for i ∈ keys(d)
            if eltype(d[i]) <: AbstractFloat
                dₒᵤₜ[i]= (vmean(d[i]), vstd(d[i]))
            end
        end
    else
        dₒᵤₜ = Dict{Symbol,Number}()
        for i ∈ keys(d)
            if eltype(d[i]) <: AbstractFloat
                dₒᵤₜ[i]= vmean(d[i])
            end
        end
    end
    return (; dₒᵤₜ...)
end


"""

```julia
distmedians(d::Dict;ci)
```
Calculate medians of each array in a `Dict` of `Array{<:Number}`.

Returns a `NamedTuple` of each key and its median (default) or a `Tuple` of `(median,lower,upper)` where `lower` and `upper` are the lower and upper bounds of the credible interval provied for `ci ∈ [0,1]`.

"""
function distmedians(d::Dict;ci::Number=0.)
    @assert 0 <= ci <= 1
    if ci>0
        α=1-ci
        dₒᵤₜ = Dict{Symbol,Tuple}()
        for i ∈ keys(d)
            if eltype(d[i]) <: AbstractFloat
                if length(d[i])>1
                    println(i)
                    med = vmedian(d[i])
                    "median was ok"
                    lower = med - vquantile(d[i],α/2)
                    upper = vquantile(d[i],1-α/2) - med
                else
                    med = d[i]
                    lower = upper = zero(med)
                end
                dₒᵤₜ[i]= (med,lower,upper)
            end
        end
    else
        dₒᵤₜ = Dict{Symbol,Number}()
        for i ∈ keys(d)
            if eltype(d[i]) <: AbstractFloat
                length(d[i])>1 ? dₒᵤₜ[i]= vmedian(d[i]) : dₒᵤₜ[i] = d[i]
            end
        end
    end
    return (; dₒᵤₜ...)
end

"""

```julia
cleanhist(x; nbins=50, scooch_nbins=4)
```
Calculates a histogram with extra (0 count) bins to buffer the edges and make it look nice and clean. 🧼

Optionally specify the number of histogram `bins` (default: 64 bins) and the number of buffering bins `scooch_nbins`. (Total bins = `nbins + scoochbins`)

Returns a `NamedTuple` with `x` and `y` values of histogram.
"""
function cleanhist(x::AbstractArray; nbins::Int=64, scooch_nbins::Int=2)

    x_scooch = scooch_nbins*(maximum(x)-minimum(x))/(nbins)
    binedges = LinRange(minimum(x)-x_scooch,maximum(x)+x_scooch,nbins+2*scooch_nbins+1)
    y=histcounts(x,binedges) ./ (length(x)*step(binedges))
    return (x=binweave(binedges), y=interleave(y),)
end


"""

```julia
summedpdfhist(x,d,ds; bins)
```

Create a histogram of the summed probability density functions (PDFs) of data with means `d` and standard deviations `ds`, respectively, calculated over domain `x`. Optionally specify the number of histogram `bins` (default: 64 bins).

Returns a `NamedTuple` with `x` and `y` values of histogram.

"""
function summedpdfhist(x::AbstractRange,d::AbstractVector,ds::AbstractVector; bins::Int=64)
# Ensure downscaling will work correctly
    rr = length(x)%bins # range residual when divied up by bins.
    if !iszero(rr)
        x = LinRange(first(x),last(x) - step(x)*rr,length(x)-rr) #first(x) : step(x) : (last(x) - step(x)*(length(x) % bins))
        @warn "input range adjusted for binning to x=$x. Use this for additional histograms. "
    end

# Prepare output range of x
    downscale_factor = div(length(x),bins)
    Δbin = downscale_factor*step(x)
    bincenters = sum(x[1:downscale_factor])/downscale_factor : Δbin : last(x)
    binnedpdf = Vector{eltype(d)}(undef,length(bincenters))

# Calculate pdf and downscale into bins
    summedpdf = sumpdfs(x,d,ds)
    ImpactChron.downscale!(binnedpdf,summedpdf)

    return (x=binweave(rangemidbounds(bincenters)), y=interleave(binnedpdf), etc=(; summedpdf, binnedpdf, bincenters, x))
end


"""

```julia
summedpdfhist!(A,d,ds)
```

In-place version of `summedpdfhist` that takes the output of a previously-run `summedpdfhist` operation (`A::NamedTuple`) as its input along with the observations `d` and `ds`. 

See also: `summedpdfhist`

"""
function summedpdfhist!(A::NamedTuple,d::AbstractVector,ds::AbstractVector)
    
# Calculate pdf and downscale into bins
    sumpdfs!(A.etc.summedpdf, A.etc.x, d,ds)
    ImpactChron.downscale!(A.etc.binnedpdf,A.etc.summedpdf)

    (x=A.x, y=interleave!(A.y,A.etc.binnedpdf), etc=A.etc)
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
function sumpdfs(z::AbstractVector,x::AbstractVector,δx::AbstractVector)
    ρ = Vector{float(eltype(x))}(undef,length(z))
	if isempty(x)
        ρ .= 0
    else
        sumpdfs!(ρ,z,x,δx)
    end
    ρ	
end


"""

```julia
sumpdfs!(ρ,z,x,δx)
```
In-place version of `sumpdfs` that overwrites a vector of probability densities `ρ` corresponding to values in `z`.

"""
function sumpdfs!(ρ::AbstractVector, z::AbstractVector, x::AbstractVector, δx::AbstractVector)
    @assert length(ρ) == length(z)
    ρ .= zero(float(eltype(ρ)))

        @tturbo for j in eachindex(x)
            xj = x[j]
            δxj = δx[j]
            for i in eachindex(z)
                ρ[i] += normdens(z[i],xj,δxj)
            end
        end
# Calculate the stepsize of z for normalization
	isa(z,AbstractRange) ? Δz = step(z) : Δz = z[2]-z[1]
# flatten and normalize z
	ρ ./= (vsum(ρ)*Δz)
end


"""

```julia
pdfsample(x::AbstractVector, p::AbstractVector; n::Int=1)
```

Draw a sample (optionally `n` samples) from a pdf corresponding to values `x` and probabilities `p`. Calculates a CDF from `p` and linearly interpolates the values of `x` from a `rand()` on the CDF.

"""
function pdfsample(x::AbstractVector, p::AbstractVector; n::Int=1)
    samples = Vector{float(eltype(x))}(undef,n)
    cp = copy(p)
    pdfsample!(samples,x,cp)
	ifelse(isone(length(samples)), samples[1], samples)
end

"""

```julia
pdfsample!(samples::AbstractVector, x::AbstractVector, p::AbstractVector)
```

In-place version of `pdfsample`. Overwrites `samples` and `p`. 

see also: `pdfsample`

"""
function pdfsample!(samples::AbstractVector, x::AbstractVector, cp::AbstractVector)
    @assert length(samples) > 1
	cumsum!(cp,cp)
	cp ./= last(cp) # ensure the cdf max = 1
	@batch for i = eachindex(samples)
		r = rand() # draw a random value from the cdf
		ir = searchsortedfirst(cp,r)
# Since searchsortedfirst returns the index ≥ the cumulative probability, interpolate between the returned step and the preceding step where the sample likely lies.
		if ir>1 # If no preceding timestep, just set to age of first timestep
			samples[i] =  x[ir-1] +  (r - cp[ir-1]) * (x[ir]-x[ir-1]) / (cp[ir]-cp[ir-1])
		else			
			samples[i] = x[ir]
		end
	end
end



"""

```julia
histhist(ht, t, d, ages, ages_sig; xbins=2^6, ybins=200, model_draws=100)
```

Calculate a histogram of the densities of the modelled cooling age distributions at each timestep. Calculates the histogram of each modeled distributions with `model_draws` resamplings, assuming uncertainties equal to the mean of `ages_sig`.

Outputs a `x`, `y`, and `z` (wrapped in a `NamedTuple`) to be used in any heatmap-style plot, with each vector of densities for a given `x` normalized to its maximum value. Also includes the `prior` with keys `x` and `y` to plot over heatmap.

| Inputs | Description | Type |
| --- | --- | --- | 
| `ht` | histogram timeseries | `Range` |
| `t` | input time parameters to model (`Δt`, `tmax`, `downscale`) | `NamedTuple` |
| `d` | logged `agedist` values from MCMC run | `Matrix` |
| `ages`, `ages_sig` | mean and 1σ of measured ages (prior) | `Vector` |

"""
function histhist(ht::AbstractRange, t::NamedTuple, d::AbstractMatrix, ages::AbstractVector, ages_sig::AbstractVector; xbins::Int=2^6, ybins::Int=200, model_draws::Int=100)

    @assert length(ages) == length(ages_sig)

    gt = 4567.3 .- ImpactChron.timemanagement(t.Δt,t.tmax,t.downscale)[2]

    prior = summedpdfhist(ht,ages,ages_sig, bins=xbins) #add bins
    m = deepcopy(prior)

    μσ = fill(sum(ages_sig)/length(ages_sig), model_draws) # mean of 1σ uncertainties on ages.
    model_samples = Vector{eltype(d)}(undef,model_draws)

    i1 = findfirst(!isnan, view(d,1,:)) # start at first accepted and logged step.

    di = Vector{eltype(d)}(undef,size(d,1))
    hists = Matrix{eltype(d)}(undef,size(m.etc.binnedpdf,1),size(d,2)-i1+1)

    printstyled("Resampling posterior date distributions: \n", color=:green)
    for i in axes(hists,2)
        i_ = i + i1 - 1 # calculate index in d
        if !isnan(d[1,i_]) 
            di .= view(d,:,i_)
            pdfsample!(model_samples,gt,di)
            summedpdfhist!(m,model_samples,μσ)
        end
        hists[:,i] .= m.etc.binnedpdf
        if iszero(i % (size(d,2)÷10)) 
            print(string("  ",100*i÷size(d,2),"% complete \r"))
            flush(stdout)
        end
    end

    ybinedges = LinRange(0,maximum(hists),ybins+1)
    ybincenters = rangemidpoints(ybinedges)
    counts = zeros(size(hists,1), size(ybincenters,1))
    printstyled("Almost there... counting each timestep: \n", color=:blue)
    @inbounds for i in axes(counts,1)
        cv = view(counts,i,:)
        histcounts!(cv, view(hists,i,:),ybinedges)
        cv ./= maximum(cv)
    end 

    (x=m.etc.bincenters, y=ybincenters, z=counts, prior)
end


println("""
Data Management Functions Loaded Successfully:

Parameter titles/variables locked in `ParamTitles`/`ParamVars`

Summary statistics: `distmeans`, `distmedians`

Visualization prep: `interleave(!)`, `binweave(!)`, `cleanhist`, `summedpdfhist(!)`, `normdens`, `sumpdfs(!)`, `pdfsample(!)`, `histhist`
""")
