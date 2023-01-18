
# Test: pertub
tss=rAlo=R=ta=cAl=Tm=Tc=ρ=Cp=k=tχα=τχα=Fχα=tχβ=τχβ=Fχβ=1.
ϕ = (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ)
ϕₜₛₜ =(; tss,rAlo,R=2.,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ)
@test perturb(ϕ,:R,2.) === ϕₜₛₜ

# Test ImpactSite Constructors
testimpactor = 6.
@test ImpactSite(Cone,testimpactor) == ImpactSite(Cone(testimpactor *(10/3),testimpactor*(5/2)), 0.)
@test ImpactSite(Cone,r=10.,C=0.1) == ImpactSite(Cone(10.,2π),0.1)
@test ImpactSite(Cone) == ImpactSite(Cone(0.,0.),0.01)