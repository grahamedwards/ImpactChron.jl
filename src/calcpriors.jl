## (Re)Calculate mean ages or recalibrate ages with different decay constants

"""

```julia
mcmean(x,xsig,n=10_000_000,fullpost=false)
```

Calculate mean of age(s) `x` with 1σ uncertainties `xsig` by Monte Carlo method. A 10M resample (default) returns consistent results at the 10ka level.

Returns a tuple of `(mean,1σ)` by default. 

If `fullpost=true`, returns posterior samples rather than summary statistics.

"""
function mcmean(x::AbstractArray,xsig::AbstractArray;n::Int=10_000_000,fullpost::Bool=false)

    out = Vector{float(eltype(x))}(undef,n)
    nₓ = length(x)
    for i in 1:n
        mcx = zero(eltype(out))
        for j in 1:nₓ
            mcx += x[j] + randn() * xsig[j]
        end
        out[i] = mcx/nₓ # Save mean of resampled dates.
    end
    mu = vmean(out)
    sig = vstd(out,mean=mu)
    return fullpost ? out : (mu, sig)
end



"""

```julia
agerecal(x,sig;monitor_age,n, KAr=false)
```
Recalculate the age and uncertainty of an Ar-Ar age `x`±`sig` (1σ) Ma with the decay constants of Steiger & Jager (1977). Quantifies uncertainty by resampling `n` times (default: 1M). Optionally accepts `monitor_age`. If not given, resamples from a range of absolute (not necessarily K-Ar constrained) monitor ages ∈ [0,3] Ga.

If `KAr=true`, this recalculates the K-Ar age with the new decay constants, without correcting for a monitor.

"""
function agerecal(x::Number,sig::Number;monitor_age::Number=0,n::Int64=1_000_000, KAr::Bool=false)
    t, s = x * 1e6, sig * 1e6 # Convert input ages to years
    λᵦ₁, λₑ₁, K₁ = 4.72e-10, .585e-10, 0.0119 # Husain 1974 and Turner 1969 (T69 uses λₑ=0.58[4]e-10), which is near enough at unc levels.
    λᵦ₂, λₑ₂, K₂ = 4.962e-10, 0.581e-10, 0.01167 # Steiger & Jager, 1977
    l₁= λᵦ₁ + λₑ₁
    l₂ = λᵦ₂ + λₑ₂
    tc = Vector{float(eltype(t))}(undef,n)

    if KAr # Recalculate the K-Ar date with the new decay constants.
        for i in 1:n
            to = t + s*randn()
            tc[i] = log((l₂/λₑ₂) * ((K₁/K₂) * λₑ₁/l₁) * (exp(to*l₁)-1) + 1) /l₂ 
        end

    elseif iszero(monitor_age) #No defined monitor age for Ar-Ar age, assume monitor age ∈ [0,3] Ga

        for i in 1:n
            mtr = 3e9 * rand() # Draw monitor age from
            #kj = (l₂/λₑ₂) * ((K₁/K₂) * λₑ₁/l₁) * (exp.(mtr*l₁)-1) / (exp(l₁*mtr)-1)
            tᵢ = t + s*randn()
            kj = (exp(l₂*mtr)-1) / (exp(l₁*mtr)-1)
            tc[i] = (1/l₂)*log(1+(kj*exp(l₁*tᵢ)-1))
        end
       
    else # Ar-Ar age with specified, K-Ar constrained monitored age.
        mtr = monitor_age
        mtr₂ = (l₂/λₑ₂) * ((K₁/K₂) * λₑ₁/l₁) * (exp(mtr*l₁)-1)
            # The adjusted monitor age = log(mtr₂ +1)/l₂
            # This is skipped, as it is cancelled out by the exp(l₂*t)-1 in the kj calculation
        kj = mtr₂ / (exp(l₁*mtr)-1)
        for i in 1:n
            tᵢ = t + s*randn()
            tc[i] = (1/l₂)*log(1+(kj*exp(l₁*tᵢ)-1))
        end
    end

    μ = vmean(tc)
    σ = vstd(tc, mean=μ)
    return μ*1e-6 ,σ*1e-6
end



"""

```julia
lognorm(x)
```

Calculate the lognormal distribution of `x`, a collection of `Number`s. Returns a `lNrm` type. Conveniently accepts `NamedTuple`s.

See also: `lognorm`

"""
function lognorm(y)
    y isa NamedTuple ? x = log.(getproperty.(Ref(y),eachindex(y))) : x=log.(y)
    μ = vmean(x)
    σ = vstd(x)
    lNrm(μ,σ)
end


"""

```julia
draw(x::Union{Number,ImpactChron.PriorDistribution})
```

Make a random draw from  `x`, which may be a `Tuple` or any subtype of `PriorDistribution``. If `x` is a `Number`, it simply returns `x`.

Used exclusively in support of `lognormMC`.

See also: `Nrm`, `lNrm`, `Unf`, `lognormMC`

"""
draw(x::Nrm; rng=Random.Xoshiro()) = x.μ + x.σ*randn(rng)
draw(x::Unf; rng=Random.Xoshiro()) = x.a+ (x.b-x.a) * rand(rng)
draw(x::Number; rng=Random.Xoshiro()) = x
draw(x::Tuple; rng=Random.Xoshiro()) = rand(rng,x)


"""

```julia
lognormMC(x ; n)
```

Calculate the lognormal distribution of a collection `x` by resampling the entire collection `n`` times (default=1M).

`x` may contain data in the form of `Tuple`s, `PriorDistribution` subtypes, or `Number`s.

Returns a `lNrm` type.

See also: `lognorm`, `ImpactChron.draw`

---
Just in case, the function has a (very slow) safety net to prevent it from trying to calculate the `log` of any negative resamples. (This has never happened for the data I used)
"""
function lognormMC(x;n::Int=1_000_000,rng=Random.Xoshiro())
    xdraws = Vector{Float64}(undef,length(x)*n)
    xinds = eachindex(x)
    @inbounds for i ∈ eachindex(xdraws)
        d = draw(x[rand(rng,xinds)],rng=rng)
        d = ifelse(d>0,d,NaN)
        isnan(d)  &&  @warn "There was a negative log! But I've ignored it for you. 😎"
        xdraws[i] = log(d)
    end
    good_draws = xdraws[.!isnan.(xdraws)]
    μ = vmean(good_draws)
    σ = vstd(good_draws)
    lNrm(μ,σ)
end