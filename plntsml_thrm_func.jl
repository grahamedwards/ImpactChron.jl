
using LoopVectorization,Polyester
using Random; #rng = MersenneTwister()
using Distributions

"""
Reference info:
t_acr = 2.13 # Ma after CAIs
    Al_conc = 1.18 /100  # Concentration of Al in chondrite (kg/kg)
    ρ  = 3210   # kg/m3 | density
    cp = 950   # J/(kg K) | specific heat capacity
    k  = 4      # W/(mK) | H chond upperbound, Yomogida & Matsui (1983)
    #kappa  = K/(ρ*cp) # m²/s | Thermal diffusivity
    rAlo = 5.23e-5 # initial solar ²⁶Al/²⁷Al | Jacobsen et al. (2008)
    To = 250;   # K, ambient temp Woolum & Cassen (1999)

const s_a  = 365.2422*24*60*60 # seconds per annum
const λ_Al = (log(2)/7.17e5)  / s_a # s⁻¹ | decay constant of Al in s
const H_Al = 0.355 # W/kg (Castillo-Rogez et al., 2009) hard


s_a  = 365.2422*24*60*60 # seconds per annum
λ_Al = (log(2)/7.17e5)  / s_a # s⁻¹ | decay constant of Al in s
H_Al = 0.355 # W/kg (Castillo-Rogez et al., 2009)

radius = 150e3 # m | Radius
Δr = 10e3 # m | radial distance increment

t_run = 100. # Ma after CAIs
Δt    = 0.1 # Ma

#### temporary
Al_conc = 1.18 /100
rAlo = 5.11e-5
ρ = 3210
t_acr = 2.13
####

time_Ma = collect( t_acr : Δt : t_run )
time_s  = time_Ma .* 1e6 .* s_a

r_rng = collect( 1000. : Δr : radius )

Aₒ = ρ*Al_conc*rAlo*H_Al*exp(-λ_Al*time_s[1]) # W/m³ |
    # Volumetric heat production at instant of accretion

T = plntsml_Tz(time_s,r_rng,
                To = 250., Ao = Aₒ, λ = λ_Al,
                K = 4., κ = 4/(ρ*950.) )

##

plot_radii = plot(time_Ma,T[1,:], label = string(Int(r_rng[1]/1000), " km"))
    for i = 2:length(r_rng)
        plot!(plot_radii,time_Ma,T[i,:], label = string(Int(r_rng[i]/1000), " km"))
    end
    plot(plot_radii, xlabel = "time (Ma after CAIs)", ylabel = "Temperature (K)")
"""
## Functions!

"""
histogramify
~~~~~~~~~~~~
histogramify(domain::AbstractRange,x::AbstractVector,y::AbstractVector;Δd::Number)

Constructs histogram over centers of `domain` from model outputs in x with corresponding
abundances in y. Requires that each value of x falls within the bounds of `domain`.



Normalizes the output, such that for output `dist` ∑ dist(xᵢ) * Δx (for each xᵢ in the bincenters of x)

Returns only the histogram masses, bincenters must be calculated externally.

REMOVE THE Δd input. This is just asking for mistakes. It takes <2 ns.

"""
function histogramify(domain::AbstractRange,x::AbstractVector,y::AbstractVector)
    dist = Vector{float(eltype(y))}(undef,length(domain)-1)
    histogramify!(dist,domain,x,y)
    return dist
end

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
# Identify indices of domain that bound all values of x
        xmin = searchsortedfirst(domain,first(xₛ)) - 1
        xmax = searchsortedlast(domain,last(xₛ)) + 1
# Ensure that xmin and xmax are defined in domain
        if iszero(xmin) || xmax > length(domain)
            printstyled("CAUTION: "; color=:yellow)
            println("Time domain bounds exceeded. Proposal rejected.")
            flush(stdout)
        elseif (xmax - xmin) > 1 # if only 1 bin filled, xmax-xmin=1
            @inbounds for i ∈ (xmin):(xmax) # 1 step outward of xmin,xmax to get peripheral values
                l = searchsortedfirst(xₛ , domain[i]) # lower index
                u = searchsortedlast(xₛ , domain[i+1]) # upper index
# Ensure values of (xₛ,yₛ) fall within bounds (if not, searchsortedfirst/searchsortedlast return l > u)
                u >= l && ( dist[i] = vreduce(+,yₛ[l:u])/Δd )
            end
        elseif (xmax - xmin) == 1
            dist[xmin] = 1.0 / Δd
        end
    end

    ∫distdx = vreduce(+,dist) * Δd
    vmap!(x -> x/∫distdx,dist,dist)

    return dist
end

function plntsml_Tz(time::AbstractArray,radii::AbstractArray;
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

    p = (; tss=tₛₛ,rAlo=rAlo,R=R,ta=tₐ,cAl=Al_conc,Tm=To,Tc=Tc,ρ=ρ,Cp=Cₚ,k=K,tχ=0.,τχ=0.,Fχ=0.)
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
    radii = LinRange(0.5*p.R/nᵣ,p.R*(1-0.5/nᵣ),nᵣ)
    PlntsmlAr!(ages,Vfrxn,peakT,p,Tmax=Tmax,Tmin=Tmin,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
    return ages,Vfrxn,radii,peakT
end

function PlntsmlAr!(ages::AbstractArray, #pre-allocated vector for cooling dates
    Vfrxn::AbstractArray, # pre-allocated vector for volume fraction of each date
    peakT::AbstractArray, # pre-allocated vector for each date's peak temperature
    p::NamedTuple;
    nᵣ::Integer,          # Number of simulated radial distances
    Δt::Number = 0.01,    # absolute timestep, default 10 ka
    tmax::Number = 2000.,  # maximum time allowed to model
    Tmax::Number = 1500.,  # maximum temperature (K, solidus after 1200C max solidus in Johnson+2016)
    Tmin::Number=0.)       # minimum temperature (K)

    Tc = p.Tc            # closure temperature
    tₛₛ = p.tss           # age of CAIs
    tₐ = p.ta            # accretion time
    R = p.R              # body radius
    rAlo = p.rAlo        # initial solar ²⁶Al/²⁷Al
    Cₚ = p.Cp            # specific heat capacity

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

    fill!(ages,NaN) # Vector{Float64}(undef,length(radii))

    time_Ma = tₐ : Δt : tmax  # time in Ma (after CAIs)
    time  = (0. : Δt : tmax - tₐ) * 1e6 * s_a # time in s (after accretion)

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
                T > Tmax && break
# Check whether the shell gets `HotEnough` (hotter than Tmin)
                HotEnough = ifelse(T > Tmin,true,false)
# Record warmest temperature yet:
                Tₚₖ = T
# compare T to Tc only if cooling (T < Tᵢ) & it got `HotEnough`
            elseif (T <= Tc) & HotEnough
                ages[i] = tₛₛ - time_Ma[j]  # log time only when T falls below Tc
                peakT[i] = Tₚₖ
                break               # kill loop
            end
            Tᵢ = T
        end # of j (time) loop
    end     # of i (radius) loop
end
"""
## rmNaN // Remove NaN code
    if !rmNaN #remove NaN option.
        dates = tₛₛ .- ages
        volumes = shell_vol/last(vols)
        radii_out = radii
    else rmNaN
        Xnan = .!isnan.(ages)
        dates = tₛₛ .- ages[Xnan]
        volumes = shell_vol[Xnan]/last(vols)
        radii_out = radii[Xnan]
    end

## ~histogramify~ WITHIN PlntsmlAr
    if time_domain == 0:0
        return dates,volumes,radii_out
    elseif length(dates) > 1 # Ensure more than 1 radius
# Return histogram of data binned by `time_domain`
        bincenters,dist = histogramify(time_domain,dates,volumes)
        return bincenters,dist,time_domain
    else
        Δtd = step(time_domain)
        bincenters = (first(time_domain) + 0.5Δtd) : Δtd : (last(time_domain) - 0.5Δtd)
        dist = zeros(length(bincenters))
        return bincenters,dist,time_domain
    end
"""

## Naive Resampler
"""
function PlntsmlRsmpl(N,P::ResampleParams;
            Δt = 0.1,      # absolute timestep, default 10 ka
            tmax = 1000.,     # maximum time allowed to model
            nᵣ = 100,        # # radial nodes
            ArDates::Function=PlntsmlAr)

    Tₒ(x::AbstractVector) = 100. * rand() + x[rand(1:40)] # sample from histogram bins of x

    dtₛₛ = Normal(P.tₛₛ.μ,P.tₛₛ.σ)       #solar system age, Ma
    drAlₒ = Normal(P.rAlₒ.μ,P.rAlₒ.σ)     # initial solar ²⁶Al/²⁷Al
    # skip P.Tm since this is done with custom function Tₒ
    dR = Uniform(P.R.a,P.R.b)   # Body radius
    dtₐ = Uniform(P.tₐ.a,P.tₐ.b)      # Accretion date, My after CAIs
    dcAl = Uniform(P.cAl.a,P.cAl.b)  # Fractional abundance of Al (g/g)
    dTc = Normal(P.Tc.μ,P.Tc.σ)       # Ar closure temperature, K
    dρ = Uniform(P.ρ.a,P.ρ.b)       # rock density, kg/m³
    dCₚ = Uniform(P.Cₚ.a,P.Cₚ.b) # Specific Heat Capacity
    dk = Uniform(P.k.a,P.k.b)          # Thermal Conductivity


    ages = Array{Float64}(undef,nᵣ,N)

    for i ∈ 1:N

        ages[:,i] = ArDates(  tₛₛ = rand(dtₛₛ),       #solar system age, Ma
                rAlo = rand(drAlₒ),     # initial solar ²⁶Al/²⁷Al
                tₐ = rand(dtₐ),      # accretion time, My after CAIs
                R = rand(dR),      # Body radius
                To = Tₒ(P.Tm),       # Disk temperature @ 2.5 au, K
                Al_conc = rand(dcAl),   # Fractional abundance of Al (g/g)
                Tc = rand(dTc),       # Ar closure temperature, K
                ρ = rand(dρ),       # rock density, kg/m³
                K = rand(dk),          # Thermal Conductivity
                Cₚ = rand(dCₚ),     # Specific Heat Capacity
                Δt = Δt,      # absolute timestep, default 10 ka
                tmax = tmax,     # maximum time allowed to model
                nᵣ = nᵣ,        # radial nodes
                rmNaN=false)[1] # allow NaNs to remain
    end
    return ages

end
"""

## Impact Heater v0.0
"""
ImpactResetAr ~ reheat volumes for an exponential impact flux
                described by parameters in p {::Proposal}
                p.tχ ~ instability start time
                p.τχ ~ e-folding timescale of impact flux
                p.Fχ ~ initial impact flux

ImpactResetAr(  dates::AbsractArray,
                Vfrxn::AbstractArray,
                p::Proposal;
                Δt::Number,tmax::Number,nᵣ::Integer)

Assumes impact heating depth of 20 km & impact reheating zone with z/D = 0.3
"""
function ImpactResetAr( dates::AbstractArray,Vfrxn::AbstractArray,p::NamedTuple;
                        #IzD::Number,   # depth/diameter ratio of impact reheating region
                        #Iz::Number,    # m | depth of impact reheating
                        Δt::Number,tmax::Number,nᵣ::Integer)
#
# FOR NOW
    IzD = 0.3
    Iz = 20e3 # m | impact heating depth
    rIz = 0.5/IzD # ratio of impact crater radius / depth

# Declare variables from input
    tₛₛ = p.tss
    tₐ =p.ta
    R = p.R # m | asteroid radius
    tᵅ = p.tχα # Ma after CAIs
    Fᵅ = p.Fχα #Initial impactor flux Ma⁻¹
    λᵅ = 1/p.τχα  # Ma⁻¹ | decay constant of impact flux
    tᵝ = p.tχβ # Ma after CAIs
    Fᵝ = p.Fχβ #Initial impactor flux Ma⁻¹
    λᵝ = 1/p.τχβ  # Ma⁻¹ | decay constant of impact flux



# Step 1: Impact Events
# Draw impact dates from flux distribution
    Itime = tₐ : Δt : tmax # timeseries for potential impacts in Ma after CAIs
    Ilog = BitVector(undef,length(Itime))

    @inbounds @batch for i ∈ eachindex(Itime)

        if (Itime[i]>tᵝ) && (Itime[i]>tᵅ) # tᵝ and tᵅ are both active
            Itᵅ = Itime[i] - tᵅ
            Itᵝ = Itime[i] - tᵝ
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
    impacts = tₛₛ .- Itime[Ilog]

# Step 2: Resetting Ar-Ar on the body
# Define depths for each radial node.
    radii = LinRange(0.5*R/nᵣ,R*(1-0.5/nᵣ),nᵣ)
# Identify base of impact-affected zone
    I_r_baseᵢ = searchsortedfirst(radii,R-Iz) # deepest reheated radius index
    I_r_base = radii[I_r_baseᵢ] #radial node at base of reheating zone.
# Calculate important node thicknesses and body volume.
    Δr = step(radii)
    Vbody = (4/3) * R^3 #note: π cancels out in I_Vfraxnᵣ calculation
# Create impact fractional volume vector
    I_Vfrxn = zeros(eltype(Vfrxn),length(impacts))

    @inbounds for imp ∈ eachindex(impacts)
        @inbounds for r ∈ I_r_baseᵢ:nᵣ
            x = (radii[r]+Iz-R) * rIz
            I_Vfrxnᵣ = x^2 * Δr / Vbody # note: π removed for cancelling out Vbody
# Only reset material that reflects primary cooling.
            I_Vfrxnᵣ = ifelse(Vfrxn[r] > I_Vfrxnᵣ, I_Vfrxnᵣ, Vfrxn[r])
# Add reheated fractions to total reheated volume
            I_Vfrxn[imp] += I_Vfrxnᵣ
# Subtract reheated fraction from Vfrxn.
            Vfrxn[r] -= I_Vfrxnᵣ
        end
    end
    return vcat(dates,impacts),vcat(Vfrxn,I_Vfrxn)
end

"...func load successful"
