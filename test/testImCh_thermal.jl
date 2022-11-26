## Test functions in ImCh_thermal.jl
    # radius_at_depth
    # Plntsml_Tz


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
