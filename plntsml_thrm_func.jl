
using LoopVectorization
#using Polyester
using Plots; gr()
using Random; rng = MersenneTwister()
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

function plntsml_Tz(time::AbstractArray,radii::AbstractArray;
    To::Float64,
    Ao::Float64,
    λ::Float64,
    K::Float64,
    κ::Float64)

    T=Array{Float64}(undef,length(radii),length(time))
    R = last(radii)
    n=1:1000

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


function PlntsmlAr(;
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
    rmNaN::Bool=true)     # Remove NaNs (never warms) from cooling history.

    κ = K / (ρ*Cₚ)
    s_a  = 365.2422 * 24.0 * 60.0 * 60.0 # seconds per annum, for Physics™!

    # Assume ²⁶Al is primary heat producer
    λ=log(2) / 7.17e5 / s_a   # ²⁶Al decay constant in s⁻¹
    H=0.355     # Specific power production of ²⁶Al (W/kg; Castillo-Rogez+2009)

    shells = LinRange(0,R,nᵣ+1)
    radii = (0.5 * R / nᵣ) .+ shells[1:end-1]

    vols = [(4.0 * π / 3.0) * z^3 for z ∈ shells]
    shell_vol = vols[2:end] .- vols[1:end-1]

    ages = fill(NaN,length(radii)) # Vector{Float64}(undef,length(radii))

    time_Ma = tₐ : Δt : tmax  # time in Ma (after CAIs)
    time  = (time_Ma .- tₐ) .* 1e6 .* s_a # time in s (after accretion)

    Aₒ = ρ * Al_conc * rAlo * H * exp(-λ * tₐ * 1e6 * s_a )

    n=1:1000 # Σ is an infinite summation, but get good returns on n=1000


    # possibly: using Polyester: @batch
    @inbounds for i = 1:length(radii)
        Tᵢ = 0.0 #
        T = 0.0
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

            if T < Tᵢ && T <= Tc    # compare T to Tc only if cooling
                ages[i]=time_Ma[j]  # log time only when T falls below Tc
                break               # kill loop
            else
                Tᵢ = T
            end
        end
    end
    if rmNaN #remove NaN option.
        Xnan = .!isnan.(ages)
        return (tₛₛ .- ages[Xnan] , shell_vol[Xnan]/last(vols) , radii[Xnan])
    else
        return (tₛₛ .- ages , shell_vol/last(vols) , radii)
    end
    #return (tₛₛ .- ages) # Return ages in geologic time (Ma)
end

## Testing plntsml_Ar output
Ar_ages = PlntsmlAr( Tc = 550,       # Ar closure temperature, K
            tₛₛ = 4567.4,    #solar system age, Ma
            tₐ = 2.13,      # accretion time, My after CAIs
            Δt = 0.1,      # absolute timestep, default 10 ka
            tmax = 1000.,     # maximum time allowed to model
            R = 150e3,      # Body radius
            nᵣ = 100,         # radial nodes
            To = 250.,       # Disk temperature @ 2.5 au, K
            Al_conc = 0.0118,   # Fractional abundance of Al (g/g)
            rAlo = 5.11e-5, # initial solar ²⁶Al/²⁷Al
            ρ = 3210,       # rock density, kg/m³
            K = 4.,          # Thermal Conductivity
            Cₚ = 950.,     # Thermal diffusivity
            rmNaN = true)

## Naive Resampler

function PlntsmlRsmpl(N,acrn::accretion_params,thrm::thermal_params;
            Δt = 0.1,      # absolute timestep, default 10 ka
            tmax = 1000.,     # maximum time allowed to model
            nᵣ = 100)         # # radial nodes

    Tₒ(x) = 100. * rand() + x[rand(1:40)] # sample from histogram bins of x
# Accretion Parameters: create distributions
    dtₛₛ = Normal(acrn.tₛₛ.μ,acrn.tₛₛ.σ)       #solar system age, Ma
    drAlₒ = Normal(acrn.rAlₒ.μ,acrn.rAlₒ.σ)     # initial solar ²⁶Al/²⁷Al
    dR = Uniform(acrn.R.a,acrn.R.b)   # Body radius
    dtₐ = Uniform(acrn.tₐ.a,acrn.tₐ.b)      # Accretion date, My after CAIs
    dcAl = Uniform(acrn.cAl.a,acrn.cAl.b)  # Fractional abundance of Al (g/g)
    # skip acrn.Tm since this is done with custom function (above)

# Thermal Parameters
    dTc = Normal(thrm.Tc.μ,thrm.Tc.σ)       # Ar closure temperature, K
    dρ = Uniform(thrm.ρ.a,thrm.ρ.b)       # rock density, kg/m³
    dCₚ = Uniform(thrm.Cₚ.a,thrm.Cₚ.b) # Specific Heat Capacity
    dk = Uniform(thrm.k.a,thrm.k.b)          # Thermal Conductivity


    ages = Array{Float64}(undef,nᵣ,N)

    for i ∈ 1:N

        ages[:,i] = PlntsmlAr(  tₛₛ = rand(dtₛₛ),       #solar system age, Ma
                rAlo = rand(drAlₒ),     # initial solar ²⁶Al/²⁷Al
                tₐ = rand(dtₐ),      # accretion time, My after CAIs
                R = rand(dR),      # Body radius
                To = Tₒ(acrn.Tm),       # Disk temperature @ 2.5 au, K
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


#many_ages = PlntsmlRsmpl(100,accret,therm, nᵣ=100)

#histogram(vec(many_ages))
