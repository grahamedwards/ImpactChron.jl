## Parameters for planetesimal thermal model


## Ar-Ar closure temperature

T_olg = [550,10]    # Trieloff+2003, oligoclase feldspar

## Parent Body Size

iss_R = [115,210] *1e3
"""
    LL_R = [150,200] *1e3 # min: Edwards&Blackburn2020 || max: Henke/Gail
    L_R  = [115,160] *1e3 # Gail+2019
    H_R  = [160,200] *1e3 #Henke+ 2013(upper),2016(lower)
    E_R = [120,210] *1e3 # Trieloff+2022 (EL)
"""

## Accretion time ~~ all uniform distributions

accr = [1.8, 2.3] # Composite of inner solar system data

## Al concentration

"""
# Lodders & Fegley 1999
Al_LL = 1.18 / 100.
Al_L  = 1.16 / 100.
Al_H  = 1.06 / 100.
Al_EH = 0.82 / 100.
Al_EL = 1.00 / 100.
Al_R  = 1.06 / 100.
    #Al_all = vcat(Al_LL,Al_L,Al_H,Al_EH,Al_EL)
"""



## Age of CAIs & initial ²⁶Al/²⁷Al Jacobsen+2008
    # Composite of CAI mineral Mg and Pb data

tₛₛ_J08 = [4567.4, 0.34*0.5]
rAlₒ_J08 = [5.11e-5 , 0.14 *0.5]
    #formerly Allende bulk CAI isochron value (5.23 ± 0.13) e-5


mutable struct accretion_params
# Background Conditions

    tₛₛ = tₛₛ_J08   # [μ , σ] ~ age of CAIs (Ma)
    rAlₒ = rAlₒ_J08 # [μ , σ] ~ initial solar ²⁶Al/²⁷Al
# Accretion Event
    R = 150e3 * ones(2) # Body radius ~ [min , max]
    tₐ = 2.13 * ones(2) # [min , max] ~ accretion time, My after CAIs
    Al_conc = 0.011 * ones(2) # [min , max] ~ Fractional abundance of Al (g/g)

# T_midplane @ 2.5 AU
    # from distribution of Woolum & Cassen, 1999
    Tm2d5 = Vector{Float64}(undef,40)
        Tm2d5[1] = 0
        Tm2d5[2:19] .= 100
        Tm2d5[20:32] .= 200
        Tm2d5[33:37] .= 300
        Tm2d5[38] = 400
        Tm2d5[39] = 500
        Tm2d5[40] = 600
end


## Thermal parameters

mutable struct thermal
# Ar-Ar closure temperature (K)
    Tc = T_olg
# Bulk density (kg/m³ | not a thermal property proper, I know)
    ρ = [3160 , 3360] # OC range, see below
# Specific heat capacity (cₚ, J kg⁻¹ K⁻¹)
    cₚ = [750 , 922]  # Range from Wach+2013 @ 873 K
# Thermal Conductivity (k, W m⁻¹ K⁻¹)
    k = [4 , 4]       # OC estimate, see below...
end

#####
# MAY NEED TO MAKE FINITE DIFFERENCE FUNCTION TO TRACK CHANGING ρ, k, Cp
    # Yomogida+Matsui1984, HeveySanders2006
#####
#


# Thermal Conductivity (k, W m⁻¹ K⁻¹)
"""
k_E = 5 #Enstatite chondrites (Opeil+2012) for T > 350 K
k_E(T) = 4.11 + 248/T

k(P) = P^(-0.96796) * (1-P) * exp(-2.0942)
    k_OC = k.([0.03,0.06])
 # Consolmagno+2008 ~ baseline porosity of 6%, use this as in situ upper bound
 # Explore down from here down to 3%
"""

# Alternatively, assume [3,4] interval for OCs and ~5 for ECs, using upperbound
# values from Opeil+2010,
"""
# Unnecessary and likely inaccurate Carbonaceous
k_CM(T) = 0.26 + 0.0013 * T # Modeled over 100 - 300 K Opeil+2010
k_CK(T) = 1.26 + 0.0011 * T # Modeled over 100 - 300 K Opeil+2010

# k_CM(T) = -0.0254 + ( 0.00563 * T ) - ( 2.07e-5 * T^2 ) + ( 3.11e-8 * T^3 )
    # # Constrained over 6 - 300 K (Opeil+2010), not good for >300 K
"""
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
