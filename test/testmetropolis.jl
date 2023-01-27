## Test functions in ImCh_metropolis.jl
    # thermochron_metropolis

using ImpactChron
using StableRNGs


ϕ = (tss=4567.3,rAlo=5.23e-5,R=log(150e3),ta=log(2.13),cAl=log(0.011),Tm=log(250),Tc=500.,ρ=log(3210),Cp=log(900),k=log(3),tχα=60., τχα=50., Fχα=3., tχβ=0., τχβ=63., Fχβ=9.)

ϕσ = (tss=.08, rAlo=0.065e-5, Tm=0.47,R=0.16, ta=.07, cAl=0.13, ρ=0.05, Cp=0.08, k=0.6, Tc=20.,tχα=15., τχα=10., Fχα=1., tχβ=1., τχβ=10., Fχβ=1.)

paramdist = (
    tss = Nrm(4567.3,.08),
    rAlo= Nrm(5.23e-5,0.065e-5),
    R   = lNrm(11.927,0.1617),
    ta  = lNrm(0.6994,0.0726),
    cAl = lNrm(-4.5665, 0.1316),
    Tm  = lNrm(5.3517,0.4691),
    Tc  = Nrm(500,75),
    ρ   = lNrm(8.112,0.0507),
    Cp  = lNrm(6.74,0.085),
    k   = lNrm(0.3319,0.6335),
    tχα  = Unf(0,799),
    τχα  = Unf(0,600),
    Fχα  = Unf(0,60),
    tχβ  = Unf(0,799),
    τχβ  = Unf(0,600),
    Fχβ  = Unf(0,60) ) 

dᵢ = 15_000.
dₑ = 10*dᵢ ; zₑ = 2*dᵢ ; ejection_shape = ImpactChron.Parabola(zₑ,dₑ/2)
dₕ = dᵢ *5 ; zₕ = dₕ * 2/3; reheat_shape = ImpactChron.Parabola(zₕ,dₕ/2)
crater= ImpactSite(Parabola,dᵢ) # (; reheat_shape)
vars = (:R,:ta,:cAl,:Tm,:Tc,:ρ,:Cp,:k,:tχα,:τχα,:Fχα,:τχβ,:Fχβ)


ages = [4548.0, 4544.0, 4541.0, 4538.0, 4533.0, 4533.0, 4532.0, 4530.0, 4530.0, 4522.0, 4520.0, 4520.0, 4520.0, 4517.0, 4514.0, 4514.0, 4511.0, 4505.0, 4505.0, 4503.0, 4500.0, 4500.0, 4500.0, 4497.0, 4495.0, 4494.0, 4490.0, 4490.0, 4490.0, 4483.0, 4480.0, 4480.0, 4480.0, 4480.0, 4480.0, 4477.0, 4470.0, 4470.0, 4469.0, 4461.0, 4460.0, 4460.0, 4454.0, 4452.0, 4450.0, 4450.0, 4450.0, 4450.0, 4450.0, 4444.0, 4440.0, 4440.0, 4435.0, 4433.0, 4430.0, 4430.0, 4430.0, 4430.0, 4420.0, 4420.0, 4411.0, 4400.0, 4400.0, 4383.0, 4380.0, 4370.0, 4360.0, 4351.0, 4350.0, 4350.0, 4340.0, 4330.0, 4313.0, 4300.0, 4249.0, 4240.0, 4230.0, 4200.0, 4180.0, 4090.0, 4005.0]
ages_1σ = [30.0, 18.0, 41.0, 13.0, 6.0, 8.0, 16.0, 20.0, 20.0, 8.0, 80.0, 10.0, 30.0, 11.0, 48.0, 20.0, 11.0, 10.0, 15.0, 52.0, 30.0, 30.0, 2.0, 9.0, 11.0, 46.0, 70.0, 30.0, 20.0, 14.0, 30.0, 30.0, 30.0, 8.0, 14.0, 20.0, 30.0, 20.0, 6.0, 8.0, 20.0, 10.0, 6.0, 9.0, 50.0, 30.0, 30.0, 30.0, 30.0, 17.0, 30.0, 40.0, 5.0, 4.0, 30.0, 30.0, 10.0, 40.0, 30.0, 20.0, 5.0, 30.0, 30.0, 10.0, 20.0, 10.0, 120.0, 8.0, 13.0, 10.0, 20.0, 40.0, 14.0, 70.0, 13.0, 20.0, 30.0, 50.0, 60.0, 40.0, 80.0]

mettest1 = thermochron_metropolis(ϕ, ϕσ, vars, ages, ages_1σ,crater,plims=paramdist, petrotypes=PetroTypes(), burnin=10, nsteps=10,  Δt= 1., downscale=10,Tmin=0.,Tmax=1373., tmax=999., nᵣ=200, updateN=10_000, archiveN=0,rng=StableRNG(4567))

@test mettest1[:rAlo]  === ϕ.rAlo
@test isapprox(mettest1[:ll],[-490.492, -475.391, -475.391, -475.391, -475.395, -475.395, -474.41, -474.41, -474.41, -474.41],atol=0.01)

# petrotypes on
test_petrotemps = ( T3=600+273., T4=700+273., T5=750+273., T6=1000+273.)
test_sample_types = vcat(fill("3",6),fill("4",9),fill("5",12),fill("6",31))
test_petrotypes =PetroTypes(test_petrotemps,test_sample_types)


mettest2 = thermochron_metropolis(ϕ, ϕσ, vars, ages, ages_1σ,crater,plims=paramdist, petrotypes=test_petrotypes, burnin=10, nsteps=10,  Δt= 1., downscale=10,Tmin=0.,Tmax=1373., tmax=999., nᵣ=200, updateN=10_000, archiveN=0,rng=StableRNG(4567))

@test isapprox(mettest2[:ll],[-473.243, -473.243, -473.243, -473.243, -473.082, -473.7, -472.512, -472.512, -472.512, -472.512], atol=0.01)