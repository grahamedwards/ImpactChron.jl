# # # # # # # # # # # # # # # # # # # # # # # # #
#            ImCh_thermal.jl                    #
#                                               #
# Thermal codes describing primary planetesimal #   
# cooling and impact reheating.                 #
#                                               #   
#   planetesimal_temperature                    #   
#   planetesimal_cooling_dates[!]               #    
#   planetesimal_cooling_timestep!              #   
#   impact_reset_array!                         #
#   radius_at_depth                             #
#   asteroid_agedist                            #
# # # # # # # # # # # # # # # # # # # # # # # # #


"""

```julia
ImpactChron.planetesimal_temperature(time::AbstractArray, radii::AbstractArray; To, Ao, λ, K, κ)
```

Calculates the evolution of temperature over `time` steps at a range of `radii` defined
for a conductively cooling sphere with thermal conductivity `K` and thermal diffusivity `κ`,
given ambient temperature `To`, initial heat production `Ao`, and heat-production decay constant `λ`.

`last(radii)` defines the radius of the sphere.

Adapted from:
Carlslaw & Jäger (1959, ISBN-13: 978-0198533689)
and
[Hevey & Sanders (2006)](http://doi.org/10.1111/j.1945-5100.2006.tb00195.x)

"""
function planetesimal_temperature(time::AbstractArray,radii::AbstractArray;
    To::Float64,
    Ao::Float64,
    λ::Float64,
    K::Float64,
    κ::Float64)

    T=Array{Float64}(undef,length(radii),length(time))
    R = last(radii)
    n=1:300

    @inbounds for i = eachindex(radii)
        @inbounds for j = eachindex(time)
            r = radii[i]
            t = time[j]

            #Σ = sum( @. ( ((-1)^n) / (n*((n^2)-( λ*(R^2)/(κ*π^2) ) ) ) ) * sin(n*π*r/R) * exp(-κ*(n^2)*(π^2)*t/(R^2)) )
            Σ = zero(Float64)
            @tturbo for nᵢ ∈ n
                α = ifelse(isodd(nᵢ), -1.0, 1.0)
                β = nᵢ*((nᵢ^2)-( λ*(R^2)/(κ*π^2) ) )
                γ = sin(nᵢ*π*r/R)
                δ = exp(-κ*(nᵢ^2)*(π^2)*t/(R^2))
                Σ += (α / β ) * γ * δ
            end

            T[i,j] = To +
            (κ*Ao/(K*λ)) * exp(-λ*t) *
            ( ( R*sin(r*(λ/κ)^0.5) / (r*sin(R*(λ/κ)^0.5)) ) - 1. ) +
            (2(R^3)*Ao/(r*K*π^3)) * Σ
        end
    end
    return T
end



## Simulate the Ar-Ar cooling dates and their abundances for a planetesimal
"""

```julia
planetesimal_cooling_dates(p::NamedTuple; nᵣ, Δt, tmax, Tmax, Tmin)
```

Returns an NTuple containing (in this order) thermochronologic cooling dates (My after CAIs) and corresponding volume fractions, radial depths (km from center), and peak temperatures (K) of `nᵣ` evenly spaced nodes in a spherical body. 

Physical and environmental parameters are described in `p`. Alternatively, these parameters may be individually listed in lieu of `p`. These parameters are outlined in the table below. 

**Note:** several of these parameters need to be entered as the natural logarithm of the value for easy compatibility with the inversion function.

`Δt` gives the timestep (in My), `tmax`` describes the duration of the model (My after CAIs = Myₛₛ), and `Tmax` and `Tmin` define the maximum and minimum temperatures (K) allowed for chondritic material in the body. Default values are only given for `tmax` (2000 Myₛₛ), `Tmax` (1500 K), and `Tmin` (0 K).


    | Parameter                 | log?  | `NmTpl`| `func`  |
    | ------------------------- | ----  | ------ | -------- |
    | solar system age (Ma)     | no    | `tss`  | `tₛₛ`    |
    | initial ²⁶Al/²⁷Al         | no    | `rAlo` | `rAlo`   |
    | closure temperature (K)   | yes   | `Tc`   | `Tc`     |
    | body radius (m)           | yes   | `R`    | `R`      |
    | accretion date (Myₛₛ)       | yes   | `ta`   | `tₐ`     |
    | disk temperature (K)      | yes   | `Tm`   | `To`     |
    | [Al] (g/g)                | yes   | `cAl`  | `Al_conc`|
    | density (kg/m³)           | yes   | `ρ`    | `ρ`      |
    | thermal diffusivity       | yes   | `k`    | `K`      |
    | specific heat capacity    | yes   | `Cp`   | `Cₚ`     |
    | ------------------------- | ----  | ------ | -------- |

see also: [`planetesimal_cooling_dates!`](@ref)
"""
function planetesimal_cooling_dates(;
            nᵣ::Integer,          # Number of simulated radial distances
            Δt::Number,    # absolute timestep in Ma
            tmax::Number = 2000.,  # maximum time allowed to model
            Tmax::Number = 1500.,  # maximum temperature (solidus after 1200C max solidus in Johnson+2016)
            Tmin::Number=0.,
            Tc::Number,
            tₛₛ::Number,
            tₐ::Number,         # accretion time
            R::Number,            # Body radius
            To::Number,           # Disk temperature (K)
            Al_conc::Number,      # Fractional abundance of Al (g/g)
            rAlo::Number,         # initial solar ²⁶Al/²⁷Al
            ρ::Number,            # rock density
            K::Number,            # Thermal Conductivity
            Cₚ::Number)           # Specific heat capacity

    p = (; tss=tₛₛ,rAlo=rAlo,R=R,ta=tₐ,cAl=Al_conc,Tm=To,Tc=Tc,ρ=ρ,Cp=Cₚ,k=K)
    ages=Array{float(typeof(tₛₛ))}(undef,nᵣ)
    Vfrxn=Array{float(typeof(R))}(undef,nᵣ)
    peakT=Array{float(typeof(To))}(undef,nᵣ)
    radii = LinRange(0.5*exp(R)/nᵣ,exp(R)*(1-0.5/nᵣ),nᵣ)
    planetesimal_cooling_dates!(ages,Vfrxn,peakT,p,Tmax=Tmax,Tmin=Tmin,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
    return ages,Vfrxn,radii,peakT
end

function planetesimal_cooling_dates(p::NamedTuple;
            nᵣ::Integer,            # Number of simulated radial distances
            Δt::Number,             # absolute timestep, default 10 ka
            tmax::Number = 2000.,   # maximum time allowed to model
            Tmax::Number = 1500.,   # maximum temperature (solidus after 1200C max solidus in Johnson+2016)
            Tmin::Number = 0.)      # minimum temperature (K)
    ages=Array{float(typeof(p.tss))}(undef,nᵣ)
    Vfrxn=Array{float(typeof(p.R))}(undef,nᵣ)
    peakT=Array{float(typeof(p.Tm))}(undef,nᵣ)
    R=exp(p.R)
    radii = LinRange(0.5*R/nᵣ,R*(1-0.5/nᵣ),nᵣ)
    planetesimal_cooling_dates!(ages,Vfrxn,peakT,p,Tmax=Tmax,Tmin=Tmin,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
    return ages,Vfrxn,radii,peakT
end



"""

```julia
planetesimal_cooling_dates!(ages, Vfrxn, peakT, p::NamedTuple; nᵣ, Δt, tmax, Tmax, Tmin)
```

In-place `planetesimal_cooling_dates` that updates Arrays `ages`, `Vfrxn`, and `peakT`.

see also: [`planetesimal_cooling_dates`](@ref)
"""
function planetesimal_cooling_dates!(ages::AbstractArray, #pre-allocated vector for cooling dates
    Vfrxn::AbstractArray, # pre-allocated vector for volume fraction of each date
    peakT::AbstractArray, # pre-allocated vector for each date's peak temperature
    p::NamedTuple;
    nᵣ::Integer,          # Number of simulated radial distances
    Δt::Number,           # absolute timestep in Ma
    tmax::Number = 2000., # maximum time allowed to model
    Tmax::Number = 1500., # maximum temperature (K, solidus after 1200C max solidus in Johnson+2016)
    Tmin::Number=0.)      # minimum temperature (K)

    Tc = p.Tc            # closure temperature
    tₛₛ = p.tss           # age of CAIs
    rAlo = p.rAlo        # initial solar ²⁶Al/²⁷Al

    Cₚ = exp(p.Cp)        # specific heat capacity
    tₐ = exp(p.ta)       # accretion time
    R = exp(p.R)         # body radius
    To = exp(p.Tm)       # disk temperature (K) (lognormally distributed)
    Al_conc = exp(p.cAl) # fractional abundance of Al (g/g) (lognormally distributed)
    ρ = exp(p.ρ)         # rock density (lognormally distributed)
    K = exp(p.k)         # thermal conductivity (lognormally distributed)
    Tc = exp(p.Tc)       # Ar closure temperature (K) (lognormally distributed)


    κ = K / (ρ*Cₚ)
    s_a  = 3.155692608e7 # seconds per annum, for Physics™!

# Assume ²⁶Al is primary heat producer
    λ = 3.0634557591238076e-14 # = log(2) / 7.17e5 / s_a   # ²⁶Al decay constant in s⁻¹
    H = 0.355     # Specific power production of ²⁶Al (W/kg; Castillo-Rogez+2009)

    shells= LinRange(zero(R),R,nᵣ+1)
    radii = LinRange(0.5*R/nᵣ,R*(1-0.5/nᵣ),nᵣ)

# Calculate proportional volumes of shells around each radial node
    Vbody = R^3
    Vo = 0.
    @inbounds for z ∈ 1:nᵣ
        Vz = shells[z+1]^3
        Vfrxn[z] = ( Vz - Vo ) / Vbody
        Vo=Vz
    end

# Time Management
    tₐ_ = ceil(tₐ/Δt) * Δt # Make the first timestep (after accretion) a multiple of Δt.
    tₐadj = tₐ_-tₐ
    time  = (tₐadj : Δt : tmax - tₐ_ ) * 1e6 * s_a # time in s (after accretion) adjusted for rounding in tₐ_

# All parts of the body are at least accretion-aged (chondrule formation was high temperature)
    fill!(ages,tₐ_) # Vector{Float64}(undef,length(radii))

# Initial ²⁶Al heat production
    Aₒ = ρ * Al_conc * rAlo * H * exp(-λ * tₐ * 1e6 * s_a )

    n=1:300 # Σ is an infinite summation, but get good returns on n=300

    @batch for i = 1:nᵣ
        Tᵢ = T = Tₚₖ = zero(To)
        HotEnough = false
        @inbounds for j = eachindex(time)

            r = radii[i]
            t = time[j]
            Σ = zero(Float64)

            @turbo for nᵢ ∈ n # tturbo -> turbo if use @batch above.
                α = ifelse(isodd(nᵢ), -1.0, 1.0)
                β = nᵢ*((nᵢ^2)-( λ*(R^2)/(κ*π^2) ) )
                γ = sin(nᵢ*π*r/R)
                δ = exp(-κ*(nᵢ^2)*(π^2)*t/(R^2))
                Σ += (α / β ) * γ * δ
            end

            T = To +
            (κ*Aₒ/(K*λ)) * exp(-λ*t) *
            ( ( R*sin(r*(λ/κ)^0.5) / (r*sin(R*(λ/κ)^0.5)) ) - 1. ) +
            (2.0*(R^3)*Aₒ/(r*K*π^3)) * Σ

# While the shell is warming...
            if T > Tᵢ
# Record warmest temperature yet:
                Tₚₖ = T
# Ensure the shell remains chondritic (does not melt):
                T > Tmax && (ages[i]=NaN; break)
# Check whether the shell gets `HotEnough` (hotter than Tmin)
                HotEnough = ifelse(T > Tmin,true,false)
# compare T to Tc only if cooling (T < Tᵢ) & it got `HotEnough`
            elseif (T <= Tc) & HotEnough
                ages[i] = tₐ_ + Δt*(j-1)  # log time only when T falls below Tc (calculated from Δt and timestep j)
                break               # kill loop
            end
            Tᵢ = T
        end # of j (time) loop
        peakT[i] = Tₚₖ
    end     # of i (radius) loop
    ages
end

"""

```julia
ImpactChron.planetesimal_cooling_timestep!(solartime::AbstractRange, time_i::Vector, Vfrxn::Vector, peakT::Vector, p::NamedTuple; nᵣ, Tmax, Tmin)
```

Returns (overwrites) thermochronologic cooling dates in `time_i` as indices of `solartime`, along with corresponding volumetric fractions (`Vfrxn`) and peak temperatures in K (`peakT`) for `nᵣ` nodes in a body with model parameters given in `p`. 
`Tmax` and `Tmin` respectively describe the maximum and minimum peak temperatures allowed. 
Failing to exceed `Tmin` returns the date of accretion, and exceeding `Tmax` sets the volumetric fraction to zero (achondritic).

see also: [`planetesimal_cooling_dates`](@ref), [`planetesimal_cooling_dates!`](@ref)

"""
function planetesimal_cooling_timestep!(solartime::AbstractRange,
    time_i::AbstractVector,
    Vfrxn::AbstractVector,
    peakT::AbstractVector,
    p::NamedTuple;
    nᵣ::Integer,          # Number of simulated radial distances
    Tmax::Number,  # maximum temperature (K, solidus after 1200C max solidus in Johnson+2016)
    Tmin::Number=0.)       # minimum temperature (K)


# Parameters
    Tc = p.Tc            # closure temperature
    tₛₛ = p.tss           # age of CAIs
    rAlo = p.rAlo        # initial solar ²⁶Al/²⁷Al

    Cₚ = exp(p.Cp)        # specific heat capacity
    tₐ = exp(p.ta)       # accretion time
    R = exp(p.R)         # body radius
    To = exp(p.Tm)       # disk temperature (K) (lognormally distributed)
    Al_conc = exp(p.cAl) # fractional abundance of Al (g/g) (lognormally distributed)
    ρ = exp(p.ρ)         # rock density (lognormally distributed)
    K = exp(p.k)         # thermal conductivity (lognormally distributed)
    Tc = exp(p.Tc)       # Ar closure temperature (K) (lognormally distributed)

    κ = K / (ρ*Cₚ)
    s_a  = 3.155692608e7 # seconds per annum, for Physics™!

# Assume ²⁶Al is primary heat producer
    λ = 3.0634557591238076e-14 # = log(2) / 7.17e5 / s_a   # ²⁶Al decay constant in s⁻¹
    H = 0.355     # Specific power production of ²⁶Al (W/kg; Castillo-Rogez+2009)

# Time Management
    tₐ_i = searchsortedfirst(solartime,tₐ)
    tₐ_ = solartime[tₐ_i]  # The first timestep (after accretion) that is a multiple of Δt.
    tₐadj = tₐ_-tₐ
    planetesimal_time  = (tₐadj : step(solartime) : last(solartime) - tₐ_ ) * 1e6 * s_a # time in s (after accretion) adjusted for rounding in tₐ_

# All parts of the body are at least accretion-aged (chondrule formation was high temperature)
    fill!(time_i,tₐ_i)

# Divide up the body into concentric `shells` with radial midpoints `radiii`
    shells= LinRange(zero(R),R,nᵣ+1)
    radii = LinRange(0.5*R/nᵣ,R*(1-0.5/nᵣ),nᵣ)

# Calculate proportional volumes of shells around each radial node
    Vbody = R^3
    Vo = 0.
    @inbounds for z ∈ 1:nᵣ
        Vz = shells[z+1]^3
        Vfrxn[z] = ( Vz - Vo ) / Vbody
        Vo=Vz
    end


# Initial ²⁶Al heat production
    Aₒ = ρ * Al_conc * rAlo * H * exp(-λ * tₐ * 1e6 * s_a )

    n=1:300 # Σ is an infinite summation, but get good returns on n=300

    @batch for i = 1:nᵣ
        Tᵢ = T = Tₚₖ = zero(To)
        HotEnough = false
        @inbounds for j = eachindex(planetesimal_time)

            r = radii[i]
            t = planetesimal_time[j]
            Σ = zero(Float64)

            @turbo for nᵢ ∈ n # tturbo -> turbo if use @batch above.
                α = ifelse(isodd(nᵢ), -1.0, 1.0)
                β = nᵢ*((nᵢ^2)-( λ*(R^2)/(κ*π^2) ) )
                γ = sin(nᵢ*π*r/R)
                δ = exp(-κ*(nᵢ^2)*(π^2)*t/(R^2))
                Σ += (α / β ) * γ * δ
            end

            T = To +
            (κ*Aₒ/(K*λ)) * exp(-λ*t) *
            ( ( R*sin(r*(λ/κ)^0.5) / (r*sin(R*(λ/κ)^0.5)) ) - 1. ) +
            (2.0*(R^3)*Aₒ/(r*K*π^3)) * Σ

# While the shell is warming...
            if T > Tᵢ
# Record warmest temperature yet:
                Tₚₖ = T
# Check whether the shell gets `HotEnough` (hotter than Tmin)
                HotEnough = ifelse(T > Tmin,true,false)
# compare T to Tc only if cooling (T < Tᵢ) & it got `HotEnough`
            elseif (T <= Tc) & HotEnough
                time_i[i] = j + tₐ_i - 1 # Adjust index to solartime
                break               # kill loop
            end
            Tᵢ = T
        end # of j (time) loop
# Record peak temperature.
        peakT[i] = Tₚₖ
# If the shell melted, exclude the achondritic volume from the body.
        Vfrxn[i] = ifelse(Tₚₖ > Tmax,zero(eltype(Vfrxn)),Vfrxn[i])
    end     # of i (radius) loop
# Renormalize volume
    chond_vol = vsum(Vfrxn)
    Vfrxn ./= chond_vol
end

"""

```julia
function impact_reset_array!(tₓr::AbstractArray,solartime::AbstractArray,tcoolₒ::AbstractArray,Vfrxn::AbstractArray,
                                impacts::AbstractArray,
                                p::NamedTuple,c::NamedTuple;
                                nᵣ::Integer,Δt::Number)
```
Simulates an impact history from -χ parameters in `p` (below), and resets
primary planetesimal cooling dates (indices of dates in `solartime` in `tcoolₒ`) and fractional volumes (`Vfraxn`)
based on impact/crater properties described in `c`.

Depth-cooling age (relative) abundances are tracked in array `tₓr` (time x radial depth),
with dimensions `(time,radius) = (length(solartime),nᵣ)`,
where `nᵣ` describes the number of radial nodes, as in the `planetesimal_cooling_dates` function.

`impacts` and `tcoolₒ` are pre-allocated vectors that respectively record
the number of impacts at each time step and the index of the primary cooling date in `solartime`.

Impact site geometries are codified by types (`Cone`,`Parabola`,`Hemisphere`)
defined in `c` and calculated in the `radius_at_depth`.

Impact flux follows an exponential decay described by parameters in `p`:
\np.tχ ~ instability start time
\np.τχ ~ e-folding timescale of impact flux
\np.Fχ ~ initial impact flux

"""
function impact_reset_array!(tₓr::AbstractArray,solartime::AbstractArray,tcoolₒ::AbstractArray,Vfrxn::AbstractArray,
                            impacts::AbstractArray,
                            p::NamedTuple,c::ImpactSite;
                            nᵣ::Integer,Δt::Number)

# Declare variables from input
    tₛₛ = p.tss
    R = exp(p.R) # m | asteroid radius
    tᵅ, tᵝ, tᵞ = p.tχα, p.tχβ, p.tχγ # Ma after CAIs
    Fᵅ, Fᵝ, Fᵞ = p.Fχα, p.Fχβ, p.Fχγ # Initial impactor flux Ma⁻¹
    λᵅ, λᵝ, λᵞ = 1/p.τχα, 1/p.τχβ, 1/p.τχγ  # Ma⁻¹ | decay constant of impact flux

# Identify some useful variables that will be reused in reheating section.
    radii = LinRange(0.5*R/nᵣ,R*(1-0.5/nᵣ),nᵣ)
    ntimes = length(solartime) # full length of time columns in timeXdepth array
    Δr = step(radii)
    Vbody = (4/3) * R^3 #note: π cancels out in I_Vfraxnᵣ calculation
    Vwhole = zero(eltype(Vfrxn)) # Sum of included elements of Vfrxn, unity if no NaN elements in `dates`

# Set all cells in tₓr to zero
    @tturbo for i ∈ eachindex(tₓr) #faster than fill!(tₓr,0.)
        tₓr[i] = zero(eltype(tₓr))
    end

# Populate each shell with primary cooling date, unless its a melted (NaN-date) layer, then fill no primary cooling date (0 instead of 1)
    @inbounds for i in 1:nᵣ
        #Vwhole += Vfrxn[i]
        tₓr[tcoolₒ[i],i] = Vfrxn[i]
    end

# Calculate "number" of impacts at each timestep
    @tturbo for i ∈ eachindex(solartime)
        Iᵅ = ifelse(solartime[i] >= tᵅ, Fᵅ*exp(-λᵅ * (solartime[i]-tᵅ) ), zero(Fᵅ) )
        Iᵝ = ifelse(solartime[i] >= tᵝ, Fᵝ*exp(-λᵝ * (solartime[i]-tᵝ) ), zero(Fᵝ) )
        Iᵞ = ifelse(solartime[i] >= tᵞ, Fᵞ*exp(-λᵞ * (solartime[i]-tᵞ) ), zero(Fᵞ) )
        impacts[i] = Δt * (Iᵅ + Iᵝ + Iᵞ)
    end

# Time to reheat this planetesimal:

### Ignore excavation for now/ever.
# At excavated depths, material is removed
#    r_baseₑ = searchsortedfirst(radii,R-c.excavate_shape.z) # deepest excavated radius index
#    @batch for r ∈ r_baseₑ:nᵣ # For each cratered radial node
#        x = ImpactChron.radius_at_depth(radii[r],R,c.excavate_shape) # Excavated radius at this depth
#        iVfrxn = x * x * Δr / (Vbody*Vfrxn[r]) # Fractional volume of shell of each excavation. note: π removed for cancelling out Vbody
#        tₒ = tcoolₒ[r] # Time index of primary cooling date
# For each timestep after primary cooling...
#        primdate = tₓr[tₒ,r]
#        @turbo for t ∈ (tₒ+1):ntimes
#            primdate -= iVfrxn * impacts[t] # Subtract crater ejecta for that layer
#        end
#        tₓr[tₒ,r]= ifelse(primdate<0, zero(iVfrxn),primdate) # Ensure that tₓr[tₒ,r] ≥ 0 (i.e. non-negative volume)
#    end

# Identify and perturb reheated depths
    r_baseₕ = searchsortedfirst(radii,R-c.heat.z) # deepest reheated radius index

    @batch for r ∈ r_baseₕ:nᵣ # radial node ||| upper limit = (r_baseₑ-1) if excavation enabled.
        x = ImpactChron.radius_at_depth(radii[r],R,c.heat) # calculate radius of reheating at this depth
        iVfrxn = x * x * Δr / (Vbody*Vfrxn[r]) # Fractional volume *of layer* reset per impact. note: π removed for cancelling out Vbody
        tₒ = tcoolₒ[r] # Time index of primary cooling date

        @inbounds for t ∈ (tₒ+1):ntimes # Model impact thermal history after primary cooling date
            if !iszero(impacts[t]) # see if there is an impact at time `t`
                reheat = impacts[t] * iVfrxn # reheated fraction of layer
                reheat = ifelse(reheat>1,one(reheat),reheat) # ensure the `reheat`ed fraction of layer does not exceed the volume of layer
                tₓr[t,r] =  reheat*Vfrxn[r] # scale reheated proportion of layer to proportion of body.

                @turbo for i ∈ tₒ:(t-1) # for each preceding timestep
                    tₓr[i,r] *= (1-reheat) # subtract the proportion `reheat`ed
                end
            end
        end
    end
# NOT NEEDED AS LONG AS ϵ ERRORS STAY SMALL, OR IF EXCAVATION IS INCORPORATED
# Finally, re-normalize everything to the body volume.
#    @tturbo for r ∈ 1:nᵣ
#        Vfᵣ = Vfrxn[r]/(Vwhole * Δt)
#        for t ∈ 1:ntimes
#            tₓr[t,r] *= Vfᵣ
#        end
#    end
end


## Impact shape functions (struct support in ImCh_parameters.jl)

"""
```julia
ImpactChron.radius_at_depth(rᵢ, R, x<:ImpactSiteShape)
```

Calculates the radius of the circle traced by a `x`-shaped region at a distance of `rᵢ` from the center of a body of radius `R`.
Note that for  `x::Hemisphere`, only its `x.r` is used.

see also: [`ImpactSiteShape`](@ref)

"""
radius_at_depth(rᵢ::Number, R::Number, x::Cone) = (rᵢ + x.z - R) * x.r / x.z # Conical approximation
radius_at_depth(rᵢ::Number, R::Number, x::Parabola) = x.r * sqrt( (rᵢ+ x.z -R) / x.z ) # Parabolic approximation
radius_at_depth(rᵢ::Number, R::Number, x::Hemisphere) = sqrt( x.r*x.r - (rᵢ-R)*(rᵢ-R) ) # Hemispheric approximation, assumes z=r


"""

```julia
asteroid_agedist!(a::AsteroidHistory, p::NamedTuple, petrotypes::PetroTypes, impactsite::ImpactSite; 
    nᵣ,Tmax,Tmin, melt_reject)
```

Calculates and updates an asteroid thermochronologic history, contained in `a`, from parameters defined in `p` for an impact-heating morphology of `impactsite`. 
Weights the abundance of ages in each radial shell by `petrotypes`.

If any petrologic types are missing (see [`ImpactChron.weight_petro_types!`](@ref)) or the radial fraction of the body exceeding Tmax `≥ melt_reject` (`0.1` by default), the age distribution is zeroed, which [`ll_dist_params`](@ref) rejects.

In short, this is a convenient wrapper for a series of functions used in the weighted thermochronologic model.

see also: [`planetesimal_cooling_timestep!`](@ref), [`planetesimal_cooling_dates`](@ref),[`ImpactChron.weight_petro_types!`](@ref),[`impact_reset_array!`](@ref), [`ImpactChron.downscale!`](@ref) 

"""
function asteroid_agedist!(a::AsteroidHistory, p::NamedTuple, petrotypes::PetroTypes, impactsite::ImpactSite{T,N}; nᵣ::Integer, Tmax::Number, Tmin::Number, melt_reject::Number=0.1) where {T,N}
# Calculate cooling history 
    planetesimal_cooling_timestep!(a.t, a.cooltime,a.Vfrxn, a.peakT, p, nᵣ=nᵣ, Tmax=Tmax, Tmin=Tmin)
# If petrologic type temperatures and abundances are included, weight accordingly.
    petrotypes.weight && ImpactChron.weight_petro_types!(a.Vfrxn,a.peakT,petrotypes)
    
# If ≥10% of interior radius melts, reject proposal
    if iszero(a.Vfrxn[ceil(Int,nᵣ*melt_reject)])
        printstyled("(meltdown) rejected\n"; color=:light_magenta);flush(stdout)
        a.agedist .= zero(eltype(a.agedist))
    else
        impactsite = ifelse( iszero(impactsite.C), impactsite, ImpactSite(T,r=exp(p.R),C=impactsite.C))
        impact_reset_array!(a.txr, a.t, a.cooltime, a.Vfrxn, a.impacts, p, impactsite, nᵣ=nᵣ,Δt=step(a.t))
        a.agedist .= vsum(a.txr,dims=2)
    end
# Downscale age distribution
    ImpactChron.downscale!(a.agedist_downscaled, a.agedist)
end