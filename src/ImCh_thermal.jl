## Thermal codes describing primary planetesimal cooling and impact reheating.
    # plntsml_Tz
    # PlntsmlAr
    # ImpactResetAr
    # Impact shape functions: cone, pbla, hemi

"""
```julia
function plntsml_temperature(time::AbstractArray, radii::AbstractArray; To::Float64, Ao::Float64, λ::Float64, K::Float64, κ::Float64 )
```

Calculates the evolution of temperature at a range of depths defined
for a conductively cooling sphere with thermal conductivity `K` and thermal diffusivity `κ`,
given ambient temperature `To`, initial heat production `Ao`, and heat-producing-element decay constant 'λ'.

Adapted from:
Carlslaw & Jäger (1959)
and
[Hevey & Sanders (2006)](http://doi.wiley.com/10.1111/j.1945-5100.2006.tb00195.x)

"""
function plntsml_temperature(time::AbstractArray,radii::AbstractArray;
    To::Float64,
    Ao::Float64,
    λ::Float64,
    K::Float64,
    κ::Float64)

    T=Array{Float64}(undef,length(radii),length(time))
    R = last(radii)
    n=1:300

    @inbounds for i = 1:length(radii)
        @inbounds for j = 1:length(time)
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

function PlntsmlAr(;
            nᵣ::Integer,          # Number of simulated radial distances
            Δt::Number = 0.01,    # absolute timestep, default 10 ka
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
    radii = LinRange(0.5*R/nᵣ,R*(1-0.5/nᵣ),nᵣ)
    PlntsmlAr!(ages,Vfrxn,peakT,p,Tmax=Tmax,Tmin=Tmin,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
    return ages,Vfrxn,radii
end

function PlntsmlAr(p::NamedTuple;
            nᵣ::Integer,          # Number of simulated radial distances
            Δt::Number = 0.01,    # absolute timestep, default 10 ka
            tmax::Number = 2000.,  # maximum time allowed to model
            Tmax::Number = 1500.,  # maximum temperature (solidus after 1200C max solidus in Johnson+2016)
            Tmin::Number = 0.)       # minimum temperature (K)
    ages=Array{float(typeof(p.tss))}(undef,nᵣ)
    Vfrxn=Array{float(typeof(p.R))}(undef,nᵣ)
    peakT=Array{float(typeof(p.Tm))}(undef,nᵣ)
    R=exp(p.R)
    radii = LinRange(0.5*R/nᵣ,R*(1-0.5/nᵣ),nᵣ)
    PlntsmlAr!(ages,Vfrxn,peakT,p,Tmax=Tmax,Tmin=Tmin,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
    return ages,Vfrxn,radii,peakT
end

"""
```julia
function PlntsmlAr!(ages::AbstractArray, Vfrxn::AbstractArray,peakT::AbstractArray, p::NamedTuple;
    nᵣ::Integer, Δt::Number, tmax::Number, Tmax::Number, Tmin::Number)
```

In-place `PlntsmlAr` that updates `ages`, `Vfrxn`, and `peakT`.

see also `PlntsmlAr`
"""

function PlntsmlAr!(ages::AbstractArray, #pre-allocated vector for cooling dates
    Vfrxn::AbstractArray, # pre-allocated vector for volume fraction of each date
    peakT::AbstractArray, # pre-allocated vector for each date's peak temperature
    p::NamedTuple;
    nᵣ::Integer,          # Number of simulated radial distances
    Δt::Number = 0.1,    # absolute timestep, default 10 ka
    tmax::Number = 2000.,  # maximum time allowed to model
    Tmax::Number = 1500.,  # maximum temperature (K, solidus after 1200C max solidus in Johnson+2016)
    Tmin::Number=0.)       # minimum temperature (K)

    Tc = p.Tc            # closure temperature
    tₛₛ = p.tss           # age of CAIs
    rAlo = p.rAlo        # initial solar ²⁶Al/²⁷Al
    Cₚ = p.Cp            # specific heat capacity

    tₐ = exp(p.ta)       # accretion time
    R = exp(p.R)         # body radius
    To = exp(p.Tm)       # disk temperature (K) (lognormally distributed)
    Al_conc = exp(p.cAl) # fractional abundance of Al (g/g) (lognormally distributed)
    ρ = exp(p.ρ)         # rock density (lognormally distributed)
    K = exp(p.k)         # thermal conductivity (lognormally distributed)


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

    # possibly: using Polyester: @batch
    @inbounds for i = 1:nᵣ
        Tᵢ = T = Tₚₖ = zero(To)
        HotEnough = false
        @inbounds for j = 1:length(time)

            r = radii[i]
            t = time[j]
            Σ = zero(Float64)

            @tturbo for nᵢ ∈ n # tturbo -> turbo if use @batch above.
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
# Ensure the shell remains chondritic (does not melt):
                T > Tmax && (ages[i]=NaN; break)
# Check whether the shell gets `HotEnough` (hotter than Tmin)
                HotEnough = ifelse(T > Tmin,true,false)
# Record warmest temperature yet:
                Tₚₖ = T
# compare T to Tc only if cooling (T < Tᵢ) & it got `HotEnough`
            elseif (T <= Tc) & HotEnough
                ages[i] = tₐ_ + Δt*(j-1)  # log time only when T falls below Tc (calculated from Δt and timestep j)
                peakT[i] = Tₚₖ
                break               # kill loop
            end
            Tᵢ = T
        end # of j (time) loop
    end     # of i (radius) loop
end


## Impact (Re)Heater v2
"""

```julia
impact_reset_array!(tₓr::AbstractArray, solartime::AbstractRange, impacts::AbstractArray, tcoolₒ::AbstractArray,
                            dates::AbstractArray,Vfrxn::AbstractArray,
                            p::NamedTuple,c::NamedTuple;
                                nᵣ::Integer,Δt::Number)
```
Simulates an impact history from -χ parameters in `p` (below), and resets Ar-Ar
primary planetesimal cooling `dates` and fractional volumes (`Vfraxn`)
based on impact/crater properties described in `c`.

Depth-cooling age (relative) abundances are tracked in array `tₓr` (time x radial depth),
with dimensions (time,radius) = (length(solartime),nᵣ),
where `nᵣ` describes the number of radial nodes, as in the `PlntsmlAr` function.

`impacts` and `tcoolₒ` are pre-allocated vectors that respectively record
the number of impacts at each time step and the index of the primary cooling date in `solartime`.

Impact site geometries are codified by types (`Cone`,`Parabola`,`Hemisphere`)
defined in `c` and calculated in the `radius_at_depth`.

Impact flux follows an exponential decay described by parameters in `p`:
\np.tχ ~ instability start time
\np.τχ ~ e-folding timescale of impact flux
\np.Fχ ~ initial impact flux

"""
function impact_reset_array!(tₓr::AbstractArray,solartime::AbstractArray,impacts::AbstractArray,tcoolₒ::AbstractArray,
                            dates::AbstractArray,Vfrxn::AbstractArray,
                            p::NamedTuple,c::NamedTuple;
                            nᵣ::Integer,Δt::Number)

# Declare variables from input
    tₛₛ = p.tss
    R = exp(p.R) # m | asteroid radius
    tᵅ = p.tχα # Ma after CAIs
    Fᵅ = p.Fχα #Initial impactor flux Ma⁻¹
    λᵅ = 1/p.τχα  # Ma⁻¹ | decay constant of impact flux
    tᵝ = p.tχβ # Ma after CAIs
    Fᵝ = p.Fχβ #Initial impactor flux Ma⁻¹
    λᵝ = 1/p.τχβ  # Ma⁻¹ | decay constant of impact flux

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
    for i in 1:nᵣ
        tᵢ = searchsortedfirst(solartime, dates[i]) # if dates[i] is NaN, tᵢ > length(solartime) and loops below won't run.
        NaNdate = isnan(dates[i])
        Vwhole += ifelse(NaNdate, zero(eltype(Vfrxn)), Vfrxn[i] )
        tᵢ = ifelse(NaNdate,ntimes,tᵢ)
        tcoolₒ[i]= tᵢ # Save index for excavation/reheating loops below
        tₓr[tᵢ,i] = ifelse(NaNdate, zero(eltype(tₓr)), one(eltype(tₓr)) ) #Place 1 at primary cooling date in timeXradius array. This will be corrected for volume later. If NaNdate, then set to 0, excluded volume.
    end

# Calculate "number" of impacts at each timestep
    @tturbo for i ∈ eachindex(solartime)
        Iᵅ = ifelse(solartime[i] >= tᵅ, Fᵅ*exp(-λᵅ * (solartime[i]-tᵅ) ), zero(Fᵅ) )
        Iᵝ = ifelse(solartime[i] >= tᵝ, Fᵝ*exp(-λᵝ * (solartime[i]-tᵝ) ), zero(Fᵝ) )
        impacts[i] = Δt * (Iᵅ + Iᵝ)
    end

# Time to reheat this planetesimal:

### Remove excavation for now.
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
    r_baseₕ = searchsortedfirst(radii,R-c.reheat_shape.z) # deepest reheated radius index

    @batch for r ∈ r_baseₕ:nᵣ # radial node ||| upper limit = (r_baseₑ-1) if excavation enabled.
        x = ImpactChron.radius_at_depth(radii[r],R,c.reheat_shape) # calculate radius of reheating at this depth
        iVfrxn = x * x * Δr / (Vbody*Vfrxn[r]) # Fractional volume *of layer* removed per impact. note: π removed for cancelling out Vbody
        tₒ = tcoolₒ[r] # Time index of primary cooling date

        @inbounds for t ∈ (tₒ+1):ntimes # Model impact thermal history after primary cooling date
            if !iszero(impacts[t]) # see if there is an impact at time `t`
                reheat = impacts[t] * iVfrxn # reheated fraction of layer
                tₓr[t,r] = reheat = ifelse(reheat>1,one(reheat),reheat) # ensure the `reheat`ed fraction of layer does not exceed the volume of layer (may be an issue for large Δt)

                @turbo for i ∈ tₒ:(t-1) # for each preceding timestep
                    tₓr[i,r] *= (1-reheat)
                end
            end
        end
    end
# Finally, re-noramlize everything to the body volume.
    @tturbo for r ∈ 1:nᵣ
        Vfᵣ = Vfrxn[r]/(Vwhole * Δt)
        for t ∈ 1:ntimes
            tₓr[t,r] *= Vfᵣ
        end
    end
end


## Impact shape functions (struct support in ImCh_parameters.jl)

"""
```julia
radius_at_depth(rᵢ::Number, R::Number, x::T) where T<:{Cone,Parabola,Hemisphere}
```

Calculates the radius of a circular area at a radial distance of rᵢ from the center of a body with radius `R`,
where the volume of the region is approximated by a cone (x::Cone), paraboloid of rotation (x::Parabola),
or hemisphere (x::Hemisphere). `x` includes a maximum depth (x.z) and surface radius (x.r).

"""
radius_at_depth(rᵢ::Number, R::Number, x::Cone) = (rᵢ + x.z - R) * x.r / x.z # Conical approximation
radius_at_depth(rᵢ::Number, R::Number, x::Parabola) = x.r * sqrt( (rᵢ+ x.z -R) / x.z ) # Parabolic approximation
radius_at_depth(rᵢ::Number, R::Number, x::Hemisphere) = sqrt( x.z*x.z - (rᵢ-R)*(rᵢ-R) ) # Hemispheric approximation, assumes r = z


"""
This is the old ImpactRestAr. Needs to be largely rewritten to a framework similar to ImpactResetArray.
Most importantly, replace probabilistic maths with number of hits at each timestep maths.

CURRENTLY NOT FUNCTIONAL


FORMER DOCUMENTATION:

ImpactResetAr ~ reheat volumes for an exponential impact flux
                described by parameters in p {::Proposal}
                p.tχ ~ instability start time
                p.τχ ~ e-folding timescale of impact flux
                p.Fχ ~ initial impact flux

ImpactResetAr(  dates::AbsractArray,
                Vfrxn::AbstractArray,
                p::NamedTuple,
                c::NamedTuple
                Δt::Number,tmax::Number,nᵣ::Integer)

Simulates an impact history from _χ parameters in `p`, and resets Ar-Ar
cooling `dates` and fractional volumes (`Vfraxn`)
based on impact/crater properties described in `c`.
Impact site morphologies may be described by a conical (`cone`),
parabolic (`pbla`), or hemispheric (`hemi`) approximation.

`Δt`, `tmax`, and `nᵣ` respectively define the timestep, model duration,
and radial nodes, as in `PlntsmlAr` function.

Note that Vfrxn is overwritten.

function impact_reset_quench( dates::AbstractArray,Vfrxn::AbstractArray,p::NamedTuple,c::NamedTuple;
                        Δt::Number,tmax::Number,nᵣ::Integer)

# Declare variables from input
    tₛₛ = p.tss
    tₒ = ceil(p.ta/Δt) * Δt #Set initial time equal to that in planetesimal thermal code.
    R = p.R # m | asteroid radius
    tᵅ = p.tχα # Ma after CAIs
    Fᵅ = p.Fχα #Initial impactor flux Ma⁻¹
    λᵅ = 1/p.τχα  # Ma⁻¹ | decay constant of impact flux
    tᵝ = p.tχβ # Ma after CAIs
    Fᵝ = p.Fχβ #Initial impactor flux Ma⁻¹
    λᵝ = 1/p.τχβ  # Ma⁻¹ | decay constant of impact flux

# Step 1: Impact Events
# Draw impact dates from flux distribution
    Itime = tₒ : Δt : tmax # timeseries for potential impacts in Ma after CAIs
    Ilog = BitVector(undef,length(Itime))

    @inbounds @batch for i ∈ eachindex(Itime)

        if (Itime[i]>tᵝ) && (Itime[i]>tᵅ) # tᵝ and tᵅ are both active
            Itᵅ = Itime[i] - tᵅ
            Itᵝ = Itime[i] - tᵝ

###########################
### THE FOLLOWING PROBABILITIES SECTION IS MATHEMATICALLY INCORRECT AND NEEDS TO BE UPDATED ###
###########################
# Calculate the union of the probablities for both impact fluxes
    # P(A∪B) = P(A)+P(B) - P(A∩B) = -(P(A)-1)*(P(B)-1)+1, as long as A,B are independent
            p_hit = Δt * ( -( Fᵝ*exp(-λᵝ*Itᵝ)-1 ) * ( Fᵅ*exp(-λᵅ*Itᵅ)-1 ) + 1 )
            Ilog[i] = ifelse( rand() < p_hit, true,false)

        elseif (Itime[i]>tᵝ) && (Itime[i]<tᵅ) # only tᵝ is active
            Itᵝ = Itime[i] - tᵝ
            p_hit = Δt * Fᵝ * exp(-λᵝ*Itime[i])
            Ilog[i] = ifelse( rand() < p_hit, true,false)

        elseif (Itime[i]<tᵝ) && (Itime[i]>tᵅ) # only tᵅ is active
            Itᵅ = Itime[i] - tᵅ
            p_hit = Δt * Fᵅ * exp(-λᵅ*Itime[i])
            Ilog[i] = ifelse( rand() < p_hit, true,false)
        else
            Ilog[i]=false
        end
    end
# Create vector of absolute impact dates from drawing
    impacts = Itime[Ilog]
    impacts .*= -1.
    impacts .+= tₛₛ # equivalent to impacts = tₛₛ .- Itime[Ilog] but one less allocation
# Step 2: Resetting Ar-Ar on the body
# Define depths for each radial node.
    radii = LinRange(0.5*R/nᵣ,R*(1-0.5/nᵣ),nᵣ)
# Identify base of impact-affected zone
    r_baseₕ = searchsortedfirst(radii,R-c.zₕ) # deepest reheated radius index
    r_baseₑ = searchsortedfirst(radii,R-c.zₑ) # deepest excavated radius index
# Calculate important node thicknesses and body volume.
    Δr = step(radii)
    Vbody = (4/3) * R^3 #note: π cancels out in I_Vfraxnᵣ calculation
# Create impact fractional volume vector
    iVfrxn = zeros(eltype(Vfrxn),length(impacts))
# Establish functions describing the shape of impact processes.
    @inbounds for imp ∈ eachindex(impacts)
# Crater excavation and volume removal
        @inbounds for r ∈ r_baseₑ:nᵣ
            x = ImpactChron.radius_at_depth(radii[r],R,c.excavate_shape) # calculate radius of excavation at this depth
            iVfrxnᵣ = x * x * Δr / Vbody # note: π removed for cancelling out Vbody
# Only remove material that's still there:
            lost = ifelse(Vfrxn[r] > iVfrxnᵣ, iVfrxnᵣ, Vfrxn[r])
# Subtract reheated fraction from Vfrxn.
            Vfrxn[r] -= lost
        end
# Sub-crater impact site reheating and Ar-Ar resetting.
        @inbounds for r ∈ r_baseₕ:(r_baseₑ-1)
            x = ImpactChron.radius_at_depth(radii[r],R,c.reheat_shape) # calculate radius of reheating at this depth
            iVfrxnᵣ = x * x * Δr / Vbody # note: π removed for cancelling out Vbody
# Only reset material that reflects primary cooling.
            reheated = ifelse(Vfrxn[r] > iVfrxnᵣ, iVfrxnᵣ, Vfrxn[r])
# Add reheated fractions to total reheated volume
            iVfrxn[imp] += reheated
# Subtract reheated fraction from Vfrxn.
            Vfrxn[r] -= reheated
        end
    end
    return vcat(dates,impacts),vcat(Vfrxn,iVfrxn)
end
"""
