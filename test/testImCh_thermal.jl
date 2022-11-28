## Test functions in ImCh_thermal.jl
    # radius_at_depth
    # Plntsml_Tz

"""
    radius_at_depth(rᵢ::Number, R::Number, x::Cone) = (rᵢ + x.z - R) * x.r / x.z # Conical approximation
    radius_at_depth(rᵢ::Number, R::Number, x::Parabola) = x.r * sqrt( (rᵢ+ x.z -R) / x.z ) # Parabolic approximation
    radius_at_depth(rᵢ::Number, R::Number, x::Hemisphere) = sqrt( x.r*x.r - (rᵢ-R)*(rᵢ-R) ) # Hemispheric approximation, assumes z=r
    nan_regolith!(d::AbstractArray,T::AbstractArray,Tmin::Number)
"""


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

# Use a new timestep for a straightforward `downscale!`-ing.
t☼_ = 0:timestep:(time_max-timestep)

# Add impact parameters to proposal.
ϕᵢ = (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k, tχα=60., τχα=50., Fχα=3., tχβ=0., τχβ=63., Fχβ=9.)

# Preallocate time x radius matrix, impact log.
A = fill(NaN,length(t☼_),nodes)
impact_log = similar(t☼_)

# Crater parameters:
dᵢ = 15_000.; dᵪ = dᵢ *5 ; zᵪ = dᵪ * 2/3; reheat_shape = ImpactChron.Parabola(zᵪ,dᵪ/2)
crater=(; reheat_shape)


impact_reset_array!(A, t☼_, tcool, vₜₛ, impact_log, ϕᵢ, crater, nᵣ=nodes, Δt=timestep)
# Downscale for a manageable answer to compare to.
downscale_factor=50
downscaled_dist = similar(A,length(t☼_)÷downscale_factor)
downscaled_t = ImpactChron.vmean(t☼[1:1+downscale_factor]):timestep*downscale_factor:last(t☼_)
ImpactChron.downscale!(downscaled_dist,ImpactChron.vsum(A,dims=2))

@test downscaled_dist ≈ [0.0001279600254134624, 0.009810062887162012, 0.021223876609076663, 0.015796896586235336, 0.009586235188802486, 0.004785997034385818, 0.0022026129373045127, 0.000982439732953423, 0.0004337960428320089, 0.00019126749954453417, 8.449129494635327e-5, 3.743315398276359e-5]