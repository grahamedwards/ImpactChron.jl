## Parameters for planetesimal thermal model
# Normally distributed data, reported as mean and 1σ
struct Nrm
    μ::Float64
    σ::Float64
end
# Normally distributed data, reported as mean and 1σ
struct lNrm
    μ::Float64
    σ::Float64
end

# Uniformly distributed data
struct Unf
    a::Float64
    b::Float64
end

mutable struct Proposal #Planetesimal Proposal
    # Background Conditions
    tss::Float64   # age of CAIs (Ma)
    rAlo::Float64 # initial solar ²⁶Al/²⁷Al
    # Accretion Event
    R::Float64 # Body radius
    ta::Float64 # Instantaneous accretion date, My after CAIs
    cAl::Float64 # Fractional abundance of Al (g/g)
    Tm::Float64
    # Thermal Parameters
    Tc::Float64  # Ar-Ar closure temperature (K)
    ρ::Float64 # Bulk density not a thermal property proper, I know)
    Cp::Float64  # Specific heat capacity (Cₚ, J kg⁻¹ K⁻¹)
    k::Float64 # Thermal Conductivity (k, W m⁻¹ K⁻¹)
end

# And a little function support for the Proposal mutable struct
Base.copy(s::Proposal)= Proposal(s.tss,s.rAlo,s.R,s.ta,s.cAl,s.Tm,s.Tc,s.ρ,s.Cp,s.k)

function Base.copyto!(s2::Proposal,s1::Proposal)
    s2.tss  = s1.tss
    s2.rAlo = s1.rAlo
    s2.R    = s1.R
    s2.ta   = s1.ta
    s2.cAl  = s1.cAl
    s2.Tm   = s1.Tm
    s2.Tc   = s1.Tc
    s2.ρ    = s1.ρ
    s2.Cp   = s1.Cp
    s2.k    = s1.k
    return
end

#vars =[:tss,:rAlo,:R,:ta,:cAl,:Tm,:Tc,:ρ,:Cp,:k]

## Resample struct: convert these to namedtuples

mutable struct ResampleParams
# Background Conditions
    tₛₛ::NamedTuple{(:μ,:σ), Tuple{Float64, Float64}}   # age of CAIs (Ma)
    rAlₒ::NamedTuple{(:μ,:σ), Tuple{Float64, Float64}}# initial solar ²⁶Al/²⁷Al
    Tm::Vector{Float64}
# Planetesimal Characteristics
    R::NamedTuple{(:a,:b), Tuple{Float64, Float64}} # Body radius
    tₐ::NamedTuple{(:a,:b), Tuple{Float64, Float64}} # Instantaneous accretion date, My after CAIs
    cAl::NamedTuple{(:a,:b), Tuple{Float64, Float64}} # Fractional abundance of Al (g/g)
# Thermal Parameters
    Tc::NamedTuple{(:μ,:σ), Tuple{Float64, Float64}}  # Ar-Ar closure temperature (K)
    ρ::NamedTuple{(:a,:b), Tuple{Float64, Float64}} # Bulk density not a thermal property proper, I know)
    Cₚ::NamedTuple{(:a,:b), Tuple{Float64, Float64}}  # Specific heat capacity (Cₚ, J kg⁻¹ K⁻¹)
    k::NamedTuple{(:a,:b), Tuple{Float64, Float64}} # Thermal Conductivity (k, W m⁻¹ K⁻¹)
end
## Epsilon skootch
rϵ = 1. + eps() # relative epsilon factor to scale "constant uniforms"

## Age of CAIs & initial ²⁶Al/²⁷Al Jacobsen+2008
    # Composite of CAI mineral Mg and Pb data
tₛₛ_J08 = (; μ= 4567.44, σ= 0.34*0.5)
rAlₒ_J08 = (; μ= 5.11e-5, σ= (0.14e-5) *0.5)
# T_midplane @ 2.5 AU
    # from distribution of Woolum & Cassen, 1999
Tm2d5 = Vector{Float64}(undef,40)
        Tm2d5[1] = 0.
        Tm2d5[2:19] .= 100.
        Tm2d5[20:32] .= 200.
        Tm2d5[33:37] .= 300.
        Tm2d5[38] = 400.
        Tm2d5[39] = 500.
        Tm2d5[40] = 600.
    #formerly Allende bulk CAI isochron value (5.23 ± 0.13) e-5
Radius = (; a= 150e3, b= 150e3 * rϵ) # Body radius ~ [min , max]
t_accr = (; a= 2.13, b= 2.13 * rϵ) # Instantaneous Accretion Date, My after CAIs
Al_conc = (; a= 0.011, b= 0.011 * rϵ) #Fractional abundance of Al (g/g)

#Thermal parameters
# Ar-Ar closure temperature (K)
T_olg = (; μ= 550., σ= 10.)    # Trieloff+2003, oligoclase feldspar
# Bulk density (kg/m³ | not a thermal property proper, I know)
ρ = (; a= 3160., b= 3360.) # OC range, see below
# Specific heat capacity (Cₚ, J kg⁻¹ K⁻¹)
Cₚ = (; a=750. , b= 922.)  # Range from Wach+2013 @ 873 K
# Thermal Conductivity (k, W m⁻¹ K⁻¹)
k = (; a=4. , b= 4. * rϵ)       # OC estimate, see below...


plntsml_params = ResampleParams(tₛₛ_J08,rAlₒ_J08,Tm2d5,Radius,t_accr,Al_conc,T_olg,ρ,Cₚ,k)


#####
# MAY NEED TO MAKE FINITE DIFFERENCE FUNCTION TO TRACK CHANGING ρ, k, Cp
    # Yomogida+Matsui1984, HeveySanders2006
#####
#

## Data Compilation & old code formats.

"""

Nrm(4567.60,0.36)   # Jacobsen2008 CAI AJEF Pb-Pb
Nrm(4567.4,0.34) # Jacobsen2008 AJEF-A34 Pb-Pb
Nrm(5.11e-5,0.14e-5) # Jacobsen2008 AJEF lowerbound

# normally distributed data, reported at mean and 1σ
mutable struct nrm
    μ::Float64
    σ::Float64
end

# range of data, on assumed uniformd distribution
mutable struct unf
    a::Float64
    b::Float64
end


## Accretion time ~~ all uniform distributions
    accr = unf(1.8, 2.3) # Composite of inner solar system data

## Parent Body Size
    iss_R = unf(115e3,210e3)

## Al content, Lodders & Fegley 1999
    Al_LL = 1.18 / 100.
    Al_L  = 1.16 / 100.
    Al_H  = 1.06 / 100.
    Al_EH = 0.82 / 100.
    Al_EL = 1.00 / 100.
    Al_R  = 1.06 / 100.
    Al_all = vcat(Al_LL,Al_L,Al_H,Al_EH,Al_EL)

## Parent Body Size
    LL_R = [150,200] *1e3 # min: Edwards&Blackburn2020 || max: Henke/Gail
    L_R  = [115,160] *1e3 # Gail+2019
    H_R  = [160,200] *1e3 #Henke+ 2013(upper),2016(lower)
    E_R = [120,210] *1e3 # Trieloff+2022 (EL)

# Thermal Conductivity (k, W m⁻¹ K⁻¹)

# OCs compiled in Opeil+ 2012, (@200 K)
Op12_comp = [1.25,3.05,0.82,1.90,0.45,1.15,3.15,2.72,2.26, # new Opeil+ 2012
    1.88,1.47, # Opeil+ 2010
    3.53,0.75,3.60,2.16,2.35,3.85,1.54,1.15,0.55,1.20,.73,.85,2.31, 1.03,2.14,1.86,.4,.47,1.54,.78,1.24,.97] # Yomogida+Matsui 1983
Op12excluded = (Piillistfer=5.51, Abee=5.33)


k_E = 5 #Enstatite chondrites (Opeil+2012) for T > 350 K
k_E(T) = 4.11 + 248/T

k(P) = P^(-0.96796) * (1-P) * exp(-2.0942)
    k_OC = k.([0.03,0.06])
 # Consolmagno+2008 ~ baseline porosity of 6%, use this as in situ upper bound
 # Explore down from here down to 3%


# Alternatively, assume [3,4] interval for OCs and ~5 for ECs, using upperbound
# values from Opeil+2010,

# Unnecessary and likely inaccurate Carbonaceous
k_CM(T) = 0.26 + 0.0013 * T # Modeled over 100 - 300 K Opeil+2010
k_CK(T) = 1.26 + 0.0011 * T # Modeled over 100 - 300 K Opeil+2010

# k_CM(T) = -0.0254 + ( 0.00563 * T ) - ( 2.07e-5 * T^2 ) + ( 3.11e-8 * T^3 )
    # # Constrained over 6 - 300 K (Opeil+2010), not good for >300 K

## Density ~ bulk densities (Consolmagno+,2008)


# OC Falls || Macke 2010, via Flynn+ 2018 (also see Ostrowski+Bryson 2019)
ρ_H  = [3.35 , 0.01]
ρ_L  = [3.30 , 0.01]
ρ_LL = [3.18 , 0.02]

# # EC Falls || Macke 2010, via Flynn+ 2018
ρ_EH = [3.58 , 0.05] # ±1σ
ρ_EL = [3.48 , 0.05]  # ±1σ

# RC || Macke 2010, via Ostrowski+Bryson 2019
ρ_R = [3.14 , 0]    # extracted w/o uncertainty from text.
                    # Actually look at p. 274 to account for wXing.
"""

"...params load successful"
