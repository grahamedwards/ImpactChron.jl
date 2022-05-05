
using LoopVectorization
using Random; rng = MersenneTwister()
using Distributions
#using Polyester?

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
converts fractional abundances `y` of model outputs `x`
into a histogram over centers of bins in `domain`

Assumes that Σy=1:
This allows us to skip calculating this value to normalize values of dist,
    so that ∫ dist dx = 1

"""
function histogramify(domain::AbstractRange,x::AbstractVector,y::AbstractVector;Δd::Number=0.)
    iszero(Δd) && (Δd=step(domain))
    dist = Vector{float(eltype(y))}(undef,length(domain)-1)
    histogramify!(dist,domain,Δd,x,y)
    return dist
end

function histogramify!(dist::AbstractVector,domain::AbstractRange,Δd::Number,x::AbstractVector,y::AbstractVector)
# Declare distribution vector
    fill!(dist,zero(eltype(dist)))
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
            for i ∈ (xmin):(xmax) # 1 step outward of xmin,xmax to get peripheral values
                l = searchsortedfirst(xₛ , domain[i]) # lower index
                u = searchsortedlast(xₛ , domain[i+1]) # upper index
# Ensure values of (xₛ,yₛ) fall within bounds (if not, searchsortedfirst/searchsortedlast return l > u)
                u >= l && ( dist[i] = vreduce(+,yₛ[l:u])/Δd )
            end
        elseif (xmax - xmin) == 1
            dist[xmin] = 1.0 / Δd
        end
    end
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
            Tmin::Number=0,
            Tc::Number,
            tₛₛ::Number,
            tₐ::Number,         # accretion time
            Δt::Number = 0.01,    # absolute timestep, default 10 ka
            tmax::Number = 2000,  # maximum time allowed to model
            R::Number,            # Body radius
            nᵣ::Integer,          # Number of simulated radial distances
            To::Number,           # Disk temperature (K)
            Al_conc::Number,      # Fractional abundance of Al (g/g)
            rAlo::Number,         # initial solar ²⁶Al/²⁷Al
            ρ::Number,            # rock density
            K::Number,            # Thermal Conductivity
            Cₚ::Number,           # Specific heat capacity
            rmNaN::Bool=false)     # Remove NaNs (never warms) from cooling history.

    ages=Array{typeof(tₛₛ)}(undef,nᵣ)
    Vfrxn=Array{typeof(R)}(undef,nᵣ)
    radii = LinRange(zero(R),R,nᵣ)
    PlntsmlAr!(ages,Vfrxn,radii,Tmin=Tmin,Tc=Tc,tₛₛ=tₛₛ,tₐ=tₐ,Δt=Δt, tmax=tmax,R=R,nᵣ=nᵣ,To=To,Al_conc=Al_conc,rAlo=rAlo,ρ=ρ,K=K,Cₚ=Cₚ,rmNaN=rmNaN)
    return ages,Vfrxn,radii
end

function PlntsmlAr!(ages::AbstractArray, #pre-allocated vector for cooling dates
    Vfrxn::AbstractArray,
    radii::AbstractRange;
    #time_domain::AbstractRange=0:0; # Ouptut time domain
    Tmin::Number=0,
    Tc::Number,
    tₛₛ::Number,
    tₐ::Number,         # accretion time
    Δt::Number = 0.01,    # absolute timestep, default 10 ka
    tmax::Number = 2000,  # maximum time allowed to model
    R::Number,            # Body radius
    nᵣ::Integer,          # Number of simulated radial distances
    To::Number,           # Disk temperature (K)
    Al_conc::Number,      # Fractional abundance of Al (g/g)
    rAlo::Number,         # initial solar ²⁶Al/²⁷Al
    ρ::Number,            # rock density
    K::Number,            # Thermal Conductivity
    Cₚ::Number,           # Specific heat capacity
    rmNaN::Bool=false)     # Remove NaNs (never warms) from cooling history.

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
        Tᵢ = 0.0 #
        T = 0.0
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

            T > Tᵢ && T>Tmin && (HotEnough = true) # Only becomes true if T exceeds Tmin while warming.
# compare T to Tc only if cooling & got hot enough
            if T < Tᵢ && T <= Tc && HotEnough
                ages[i]= tₛₛ - time_Ma[j]  # log time only when T falls below Tc
                break               # kill loop
            else
                Tᵢ = T
            end
        end
    end
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

"...func load successful"
