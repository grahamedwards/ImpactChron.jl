## Parameters & parameter functions for planetesimal thermal model & metropolis code.
    # PriorDistribution
        # Nrm
        # lNrm
        # Unf
    # ImpactSiteShape
        # Cone
        # Parabola
        # Hemisphere
    # ImpactSite
    # PetroTypes
        #TempProp
    # AsteroidHistory
        # timemanagement
    # perturb

## Distribution structs
"""

```julia
PriorDistribution
```

Supertype of [`Nrm`](@ref), [`lNrm`](@ref), [`Unf`](@ref)

"""
abstract type PriorDistribution end

"""

```julia
Nrm(μ::Float64,σ::Float64)
```

Immutable `struct` to describe normally distributed data, reported as mean (`μ`) and 1σ (`σ`)

"""
struct Nrm <: PriorDistribution
    μ::Float64
    σ::Float64
end

"""

```julia
lNrm(μ::Float64,σ::Float64)
```

Immutable `struct` to describe lognormally distributed data, reported as natural-log-space mean (`μ`) and 1σ (`σ`)

"""
struct lNrm <: PriorDistribution
    μ::Float64
    σ::Float64
end

"""

```julia
Unf(a::Float64,b::Float64)
```

Immutable `struct` to describe uniformly distributed data, reported as minimum (`a`) and maximum (`b`).

"""
struct Unf <: PriorDistribution
    a::Float64
    b::Float64
end


## Structure support for radius_at_depth
"""

```julia
ImpactSiteShape
```

Supertype of `Cone`, `Parabola`, and `Hemisphere`.
Each contains a field `z` (depth) and `r` (radius), both in meters.

see also: [`ImpactSite`](@ref)

"""
abstract type ImpactSiteShape end

struct Cone{T<:Number} <: ImpactSiteShape
    z::T
    r::T
end

struct Parabola{T<:Number} <: ImpactSiteShape
    z::T
    r::T
end

struct Hemisphere{T<:Number} <: ImpactSiteShape
    z::T
    r::T
end

Hemisphere(x::Number) = Hemisphere(x,x)

struct ImpactSite{T<:ImpactSiteShape, N<:Number}
    heat::T
    #eject::ImpactSiteShape # an ejecta volume may be added some day.
    C::N
end


"""

```julia
ImpactSite(heat<:ImpactSiteShape, C<:Number)
```

`struct` describing the scale and shape of the simulated volume of impact heating.

---

CONSTRUCTOR FUNCTION
====================

```julia
ImpactSite(shape, impactor_diameter)
```

Providing an `impactor_diameter` (`<:Number`) calculates impact parameters based on approximate heat distribution modeled in Davison+ 2012 (GCA, http://doi.org/10.1016/j.gca.2012.08.001).

```julia
ImpactSite(shape; r, C)
```

If values of `r` and `C` are provided,  prepares an `ImpactSite` that extends to the center of an asteroid with radius `r` (km) and has a surface diameter `C` times the asteroid circumference (`C ∈ [0,1]`, `C = 0.01` by default).
If no `r` is provided this seeds an `ImpactSite` with zeroed `ImpactSiteShape` and a `C` value. 

"""
function ImpactSite(::Type{T}, impactor_diameter::Number) where {T<:ImpactSiteShape}

# Proportions of ejection site.
    #ejection_diameter=10 * impactor_diameter 
    #ejection_depth = 2 * impactor_diameter

    heated_diameter = 5 * impactor_diameter
    heated_depth = (2/3) * heated_diameter

    ImpactSite(T(heated_depth,heated_diameter/2),zero(float(impactor_diameter)))
end

function ImpactSite(::Type{T}; r::N=0., C::N=0.01 ) where {T<:ImpactSiteShape, N<: Number}
    @assert 0 ≤ C ≤ 1
    ImpactSite(T(r, 2π * r * C),C)
end



## Petrologic Type Temperatures and Proportions

struct TempProp{T<:AbstractFloat}
    T::T
    p::T
end


"""

```julia
PetroTypes(temps, samples)
```

`struct` containing fields of `type3`–`type6`, reflecting petrologic types. Each field contains a `ImpactChron.TempProp` with subfields `T` (maximum temperature in K) and `p` (proportion among the chondrite record). Note that `p`s sum to unity.

*Only build `PetroTypes` with its constructor function*

Constructor takes a `NamedTuple` containing maximum temperatures as fields `T3`-`T6` and a `Vector{String}` containing the petrologic types corresponding to chondrites used as a prior.

e.g. `PetroTypes( (T3 = 873, T4 = 973, T5 = 1073, T6 = 1223), ["4", "6", "3", "5,6", "im"])`

If no argument given --  `PetroTypes()` -- returns zeroed `type_` fields and `weight=false`, which short-circuits weighting by petrologic type. 

"""
struct PetroTypes{T<:AbstractFloat}
    type3::TempProp{T}
    type4::TempProp{T}
    type5::TempProp{T}
    type6::TempProp{T}
    weight::Bool
end

function PetroTypes(temps::NamedTuple,samples::Vector{String})
    n3 = count(contains.(samples,"3"))
    n4 = count(contains.(samples,"4"))
    n5 = count(contains.(samples,"5"))
    n6 = count(contains.(samples,"6"))

    n = n3+n4+n5+n6

    t3 = ImpactChron.TempProp(float(temps.T3),n3/n)
    t4 = ImpactChron.TempProp(float(temps.T4),n4/n)
    t5 = ImpactChron.TempProp(float(temps.T5),n5/n)
    t6 = ImpactChron.TempProp(float(temps.T6),n6/n)

    @assert t3.p+t4.p+t5.p+t6.p ≈ 1
    return PetroTypes(t3,t4,t5,t6,true)
end

function PetroTypes()
    x = ImpactChron.TempProp(0.,0.)
    PetroTypes(x,x,x,x,false)
end

"""

```julia
AsteroidHistory(typeseed<:Number; nnodes, Δt, tmax, downscale_factor)
```

`struct` containing `Array`s that record the evolution of a (bombarded) asteroid.

Constructor function takes a parameter from the proposals to seed type and requires the number of radial nodes `nnodes` (`::Int`), the timestep used `Δt`, the full time duration `tmax`, and the `downscale_factor` (`::Int`).

Fields in AsteroidHistory:

| Field | Description |
| :---- | :---------- |
|`Vfrxn` | volume fractions of each radial shell |
|`peakT` | peak temperature of each radial shell | 
|`cooltime` | indices in `t` of the primary cooling date |
|`impacts` | number of impacts at each timestep |
|`txr` | time x radius array of proportional cooling ages |
|`agedist` | distribution of ages corresponding to `t` |
|`agedist_downscaled` | distribution of ages corresponding to `t_downscaled` |
|`t` |  timesteps of full-scale model |
|`t_downscaled` | timesteps of downscaled model output |

"""
struct AsteroidHistory{T<:Number}

# Vectors output from planetesimal cooling model
    Vfrxn::Vector{T}
    peakT::Vector{T}
    cooltime::Vector{Int}
# Arrays recording impacts and reheating
    impacts::Vector{T}
    txr::Matrix{T}
# Age distributions
    agedist::Vector{T}
    agedist_downscaled::Vector{T}
# Time
    t::AbstractRange{T}
    t_downscaled::AbstractRange{T}
end
    
function AsteroidHistory(::T; nnodes::Int, Δt::Number, tmax::Number, downscale_factor::Int ) where {T<:Number}

# Vectors that correspond to radial nodes:
    Vfrxn = Vector{T}(undef,nnodes) # volumetric fraction of each concentric node in asteroid
    peakT = Vector{T}(undef,nnodes) # peak temperature of each concentric node
    cooltime = Vector{Int}(undef,nnodes) # tracker of indices of primary cooling date in timeseries

timerange, time_downscaled = timemanagement(Δt,tmax,downscale_factor)

# Vectors that correspond to timesteps:
    impacts = Vector{T}(undef,length(timerange)) # tracker of # of impacts at each timestep
    txr = Matrix{T}(undef,length(timerange),nnodes) # time x radial position array to be used in impact resetting scheme
    agedist = Vector{T}(undef,length(timerange)) # distribution of ages relative to timeseries (time_r, time_v)
    agedist_downscaled = Vector{T}(undef,length(time_downscaled)) # distribution of ages relative to downscaled timeseries

    return AsteroidHistory(Vfrxn, peakT, cooltime, impacts, txr, agedist, agedist_downscaled, timerange, time_downscaled)
end


"""

```julia
ImpactChron.timemanagement(Δt, tmax, downscale_factor::Int)
```

Calculate (and adjust if necessary) the model timescales for a given `tmax` (My after CAIs), `Δt` timestep (in My), and a `downscale_factor`. 
May overwrite `tmax` to ensure downscaling works properly.

Returns a `Tuple` containing the timescale and the downscaled timescale.
"""
function timemanagement(Δt::Number, tmax::Number, downscale_factor::Int)
    downscale_adj = length(0:Δt:tmax)%downscale_factor
    tmax = tmax - downscale_adj * Δt
    iszero(downscale_adj) || @warn "time range adjusted for downscale to 0:Δt(=$Δt):$tmax"
    timerange = 0.0:Δt:tmax
    time_downscaled = sum(timerange[1:downscale_factor])/downscale_factor : Δt*downscale_factor : tmax
    timerange, time_downscaled
end

## Function(s)
"""
```julia
perturb(p::NamedTuple,k::Symbol,n::Number)
```
Return a NamedTuple identical to `p`,
with one field (key `k`) changed to the value of `n`.
Note that `==` identity is preserved only if the
order of fields in `p` is: `tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ,tχγ,τχγ,Fχγ`

"""
function perturb(p::NamedTuple,k::Symbol,n::Number)
    tχγ  = ifelse(k==:tχγ,n,p.tχγ)
    τχγ  = ifelse(k==:τχγ,n,p.τχγ)
    Fχγ  = ifelse(k==:Fχγ,n,p.Fχγ)
    tχβ  = ifelse(k==:tχβ,n,p.tχβ)
    τχβ  = ifelse(k==:τχβ,n,p.τχβ)
    Fχβ  = ifelse(k==:Fχβ,n,p.Fχβ)
    tχα  = ifelse(k==:tχα,n,p.tχα)
    τχα  = ifelse(k==:τχα,n,p.τχα)
    Fχα  = ifelse(k==:Fχα,n,p.Fχα)
    tss  = ifelse(k==:tss,n,p.tss)
    rAlo = ifelse(k==:rAlo,n,p.rAlo)
    R    = ifelse(k==:R,n,p.R)
    ta   = ifelse(k==:ta,n,p.ta)
    cAl  = ifelse(k==:cAl,n,p.cAl)
    Tm   = ifelse(k==:Tm,n,p.Tm)
    Tc   = ifelse(k==:Tc,n,p.Tc)
    ρ    = ifelse(k==:ρ,n,p.ρ)
    Cp   = ifelse(k==:Cp,n,p.Cp)
    k    = ifelse(k==:k,n,p.k)
    (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ,tχγ,τχγ,Fχγ)
end
