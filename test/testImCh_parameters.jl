
# Test: pertub
tss=rAlo=R=ta=cAl=Tm=Tc=ρ=Cp=k=tχα=τχα=Fχα=tχβ=τχβ=Fχβ=1.
ϕ = (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ)
ϕₜₛₜ =(; tss,rAlo,R=2.,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ)
@test perturb(ϕ,:R,2.) === ϕₜₛₜ
