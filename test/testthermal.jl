## Test functions in ImCh_thermal.jl
    # radius_at_depth
    #

"""
    radius_at_depth(rᵢ::Number, R::Number, x::Cone) = (rᵢ + x.z - R) * x.r / x.z # Conical approximation
    radius_at_depth(rᵢ::Number, R::Number, x::Parabola) = x.r * sqrt( (rᵢ+ x.z -R) / x.z ) # Parabolic approximation
    radius_at_depth(rᵢ::Number, R::Number, x::Hemisphere) = sqrt( x.r*x.r - (rᵢ-R)*(rᵢ-R) ) # Hemispheric approximation, assumes z=r
    nan_regolith!(d::AbstractArray,T::AbstractArray,Tmin::Number)
"""

using ImpactChron # safetestset

tss=4567.3
rAlo=5.23e-5
R=log(150e3)
ta=log(2.13)
cAl=log(0.011)
Tm=log(250)
Tc=500.
ρ=log(3210)
Cp=log(900)
k=log(3)

ϕ = (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k)

nodes = 10
temp_min = 800.
temp_max = 1194.9
time_max = 600.
timestep = 1.

d,v,r,T = planetesimal_cooling_dates(Δt = timestep,tmax = time_max,nᵣ = nodes,Tmin=temp_min,Tmax=temp_max,tₛₛ=tss,rAlo=rAlo,R=R,tₐ=ta,Al_conc=cAl,To=Tm,Tc=Tc,ρ=ρ,Cₚ=Cp,K=k,)
dₙₜ,vₙₜ,rₙₜ,Tₙₜ =  planetesimal_cooling_dates(ϕ,Δt = timestep,tmax = time_max,nᵣ = nodes,Tmin=temp_min,Tmax=temp_max)

@test isnan(d[1])
@test isnan(dₙₜ[1])
@test d[2:end] ≈ dₙₜ[2:end] ≈ [142.0, 137.0, 130.0, 119.0, 106.0, 88.0, 64.0, 36.0, 3.0]
@test v ≈ vₙₜ ≈ [0.001, 0.007, 0.019, 0.037, 0.061, 0.091, 0.127, 0.169, 0.217, 0.271]
@test r ≈ rₙₜ ≈ LinRange(0.5*exp(R)/length(d),exp(R)*(1-0.5/length(d)),length(d))
@test T ≈ Tₙₜ ≈ [1194.9083889204526, 1194.895971089976, 1194.847287519626, 1194.6756232947826, 1194.1652222467533, 1192.4078388987657, 1185.867523082089, 1161.4349530706934, 1077.7819961534767, 732.8385434324932]

t☼ = 0:timestep:time_max
tcool = Vector{Int}(undef,nodes)
Tₜₛ = fill(NaN,nodes)
vₜₛ = copy(Tₜₛ)

ImpactChron.planetesimal_cooling_timestep!(t☼,tcool,vₜₛ,Tₜₛ,ϕ, nᵣ=nodes,Tmax=temp_max,Tmin=temp_min)

@test t☼[tcool] ≈ [144.0, 142.0, 137.0, 130.0, 119.0, 106.0, 88.0, 64.0, 36.0, 3.0]
@test vₜₛ ≈ [0.0, 0.007007007007007008, 0.019019019019019014, 0.037037037037037056, 0.06106106106106105, 0.09109109109109105, 0.12712712712712712, 0.1691691691691693, 0.2172172172172171, 0.27127127127127126]
@test Tₜₛ ≈ [1194.9091253070367, 1194.895971089976, 1194.847287519626, 1194.6756232947826, 1194.1652222467533, 1192.4078388987657, 1185.867523082089, 1161.4349530706934, 1077.7819961534767, 732.8385434324932]

## Test impact reheating

# Use a new timestep. 
t☼_ = 0:20:400

# Add impact parameters to proposal.
ϕᵢ = (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k, tχα=60., τχα=50., Fχα=3., tχβ=0., τχβ=63., Fχβ=9.)

# Preallocate time x radius matrix, impact log.
A = fill(NaN,length(t☼_),nodes)
impact_log = similar(t☼_)

# Crater parameters:
dᵢ = 15_000.; dᵪ = dᵢ *5 ; zᵪ = dᵪ * 2/3; reheat_shape = ImpactChron.Parabola(zᵪ,dᵪ/2)
crater= ImpactSite(Parabola,dᵢ) #(; reheat_shape)

# Overwrite primary cooling with new timesteps
ImpactChron.planetesimal_cooling_timestep!(t☼_,tcool,vₜₛ,Tₜₛ,ϕ, nᵣ=nodes,Tmax=temp_max,Tmin=temp_min)

impact_reset_array!(A, t☼_, tcool, vₜₛ, impact_log, ϕᵢ, crater, nᵣ=nodes, Δt=step(t☼_))

@test isapprox(vec(ImpactChron.vsum(A,dims=2)), [0.0, 0.0, 0.0, 0.0, 0.03904, 0.18619, 0.23539, 0.14912, 0.09822, 0.07544, 0.06407, 0.04467, 0.03404, 0.02867, 0.0147, 0.01508, 0.00764, 0.00773, 0.0, 0.0, 0.0], atol = 2e-5)


# Combining everything into `asteroid_agedist!`
AH = AsteroidHistory(ϕᵢ.R, nnodes=nodes, Δt=timestep,tmax=600, downscale_factor=50)

# Melted Condition
asteroid_agedist!(AH,ϕᵢ,PetroTypes(),crater; nᵣ=nodes, Tmax=temp_max, Tmin=temp_min, melt_reject=0.1)

@test AH.agedist_downscaled == zero(AH.agedist_downscaled)

# Impact Reset
asteroid_agedist!(AH,ϕᵢ,PetroTypes(),crater; nᵣ=nodes, Tmax=1300, Tmin=temp_min, melt_reject=0.1)

@test isapprox(AH.agedist_downscaled, [4.1411e-5, 0.0042482, 0.0079689, 0.0035466, 0.0021881, 0.0010999, 0.00050767, 0.00022672, 0.00010016, 4.4173e-5, 1.9515e-5, 8.6465e-6], rtol=1e-4)

# One Impact Flux
ϕ₁ = perturb(ϕᵢ,:Fχα,0.)
asteroid_agedist!(AH,ϕ₁,PetroTypes(),crater; nᵣ=nodes, Tmax=1300, Tmin=temp_min, melt_reject=0.1)

@test isapprox(AH.agedist_downscaled, [0.00030836, 0.0057676, 0.0079282, 0.002803, 0.0016347, 0.00083097, 0.00039634, 0.00018361, 8.3942e-5, 3.8147e-5, 1.7288e-5, 7.8255e-6], rtol=1e-4) 

# No Impacts
ϕ₂ = perturb(ϕ₁,:Fχβ,0.)
asteroid_agedist!(AH,ϕ₂,PetroTypes(),crater; nᵣ=nodes, Tmax=1300, Tmin=temp_min, melt_reject=0.1)

@test isapprox(AH.agedist_downscaled, [0.00976, 0.00592, 0.00432, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], rtol=1e-4)