## Parameters & parameter functions for planetesimal thermal model & metropolis code.

"""
```julia
Nrm(μ::Float64,σ::Float64)
```
Immutable `struct` to describe normally distributed data,
    reported as mean (`μ`) and 1σ (`σ`)
"""
struct Nrm
    μ::Float64
    σ::Float64
end

"""
```julia
lNrm(μ::Float64,σ::Float64)
```
Immutable `struct` to describe lognormally distributed data,
    reported as log-space mean (`μ`) and 1σ (`σ`)
"""
struct lNrm
    μ::Float64
    σ::Float64
end

"""
```julia
Unf(a::Float64,b::Float64)
```
Immurable `struct` to describe uniformly distributed data,
    reported as minimum (`a`) and maximum (`b`).
"""
struct Unf
    a::Float64
    b::Float64
end

"""
```julia
perturb(p::NamedTuple,k::Symbol,n::Number)
```
Return a NamedTuple identical to `p`,
with one field (key `k`) changed to the value of `n`.
Note that `==` identity is preserved only if the
order of fields in `p` is as below

Fields: `tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ`
"""
function perturb(p::NamedTuple,k::Symbol,n::Number)
    tχβ  = ifelse(k==:tχβ,n,p.tχβ)
    τχβ  = ifelse(k==:τχβ,n,p.τχβ)
    Fχβ  = ifelse(k==:Fχβ,n,p.Fχβ)
    tχα  = ifelse(k==:tχα,n,p.tχα)
    τχα  = ifelse(k==:τχα,n,p.τχα)
    Fχα  = ifelse(k==:Fχα,n,p.Fχα)
    tss  = ifelse(k==:tss,n,p.tss)
    rAlo = ifelse(k==:rAlo,n,p.rAlo)
    R    = ifelse(k==:R,n,p.R)
    ta   = ifelse(k==:ta,n,p.ta)
    cAl  = ifelse(k==:cAl,n,p.cAl)
    Tm   = ifelse(k==:Tm,n,p.Tm)
    Tc   = ifelse(k==:Tc,n,p.Tc)
    ρ    = ifelse(k==:ρ,n,p.ρ)
    Cp   = ifelse(k==:Cp,n,p.Cp)
    k    = ifelse(k==:k,n,p.k)
    (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ)
end
