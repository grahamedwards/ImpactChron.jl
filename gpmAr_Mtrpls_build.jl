## gpmAr Metropolis functions: construction site

## Parameters & parameter functions for metropolis code.

"""
```julia
lNrm(μ::Float64,σ::Float64)
```
Immutable struct to describe normally distributed data,
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
Immutable struct to describe lognormally distributed data,
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
Immurable struct to describe uniformly distributed data,
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


## Plot Evolution of proposals
function plotproposals(d::Dict,plims::NamedTuple,cols::Integer;vars::Tuple=(),
                        ll::Bool=true,bounds::Bool=true)
    isempty(vars) ? v=keys(plims) : v=vars
    ll && (v=tuple(:ll,v...))
    nᵥ=length(v)

#Convert to strings if necessary
    isequal(eltype(keys(d)),String) ? (v= String.(v); llₛ = "ll" ; acpt = "accept") : (llₛ = :ll; acpt = :accept)
# Calculate number of rows needed to accomodate all variables in `cols` columns.
    rows = Int(ceil(nᵥ/cols,digits=0))

    panels = Vector{Any}(nothing,nᵥ)
    for i ∈ 1:nᵥ
        k = v[i]
        y = d[k]
        x = 1:length(y)
        if k == llₛ
            panels[i] = plot(x,y,xticks=[],ylabel="ll",linecolor=:black) #use \scrl eventually
            r = 100 * round(sum(d[acpt])/length(d[acpt]),digits=3)
            xlabel!("acceptance = $r %",xguidefontsize=6)
            #annotate!(last(x), (y[end]+y[1])/2, text("acceptance = $r %", :black,:bottomleft,6))
        elseif isnan(last(y))
            panels[i] = plot([1,last(x)],fill(y[1],2),xticks=[],ylabel="$k",linecolor=:black)
        else
            panels[i] = plot(x,y,xticks=[],linecolor=:black)
        end

        if bounds && k != llₛ
            B = plims[Symbol(k)]
            if isa(B,Unf)
                plot!([1,last(x)],fill(B.a,2),ylabel="$k",linecolor=:grey,linestyle=:solid)
                plot!([1,last(x)],fill(B.b,2),linecolor=:grey,linestyle=:solid)
            elseif isa(B,Nrm)
                plot!([1,last(x)],fill(B.μ+B.σ,2),ylabel="$k",linecolor=:grey,linestyle=:dash)
                plot!([1,last(x)],fill(B.μ-B.σ,2),linecolor=:grey,linestyle=:dash)
            elseif isa(B,lNrm)
                plot!([1,last(x)],fill(B.μ+B.σ,2),ylabel="log[" * "$k" * "]",linecolor=:grey,linestyle=:dashdot)
                plot!([1,last(x)],fill(B.μ-B.σ,2),linecolor=:grey,linestyle=:dashdot)
            end
        end
    end
    sbplts=rows*cols
    Δplts = sbplts-length(panels)
    if Δplts > 0
        blnkplt = plot(legend=false,grid=false,foreground_color_subplot=:white)
        [ push!(panels,blnkplt) for j ∈ 1:Δplts]
    end

    plot(panels...,layout=grid(rows,cols),labels="")
end

## Log-likelihood calculation

"""
ll_dist(x::AbstractVector,dist::AbstractVector,mu::AbstractVector,sigma::AbstractVector)

where `x` contains the bincenters of a normalized histogram `dist`,
and the vectors mu and sigma respectively contain the mean and 1σ of the observations.

"""
function ll_dist(   x::AbstractRange, dist::AbstractVector,#p_dist::Tuple{AbstractVector,AbstractVector},    # Proposed distribution (ages,proportions)
                    mu::AbstractVector,      # 1D Array/Vector of observed μ's (sorted)
                    sigma::AbstractVector)   # 1D Array/Vector of observed σ's (sorted)

# Define some frequently used variables
    ll = zero(float(eltype(dist)))
    nₘ = length(mu)            # n of "Measured" data
    nₚ = length(x)     # n of bincenters
    nbtwns = nₚ - 1   # n of spaces between bin centers.
    Δx = step(x)
    xᵣ = abs(last(x) - first(x)) # range of x
    #tkr = zeros(nₘ)

# Cycle through each datum in (mu,sigma)
    @inbounds for j ∈ 1:nₘ #might batch speed this up?
# Find index of μ in the `dist` array
        iₓ = ceil(Int, (mu[j]-x[1]) / xᵣ * nbtwns + 1)
                # formerly iₓ = searchsortedfirst(x,mu[j]) # x[iₓ] ≥ mu[j]
# If possible, prevent aliasing problems by interpolation
        if (iₓ>1) && (iₓ<=nₚ) && ( (2sigma[j]) < (x[iₓ]-mu[j]) ) && ( (2sigma[j])<(mu[j]-x[iₓ-1]) )
            # && (sigma[j] < (x[iₓ]-x[iₓ-1]) ) # original threshold, see notes.
# Interpolate corresponding distribution value, note: Δx cancels in second term
            likelihood = dist[iₓ] * Δx - (x[iₓ]-mu[j]) * (dist[iₓ]-dist[iₓ-1])
            #likelihood = 6 * sigma[j] * ( dist[iₓ] - (dist[iₓ]-dist[iₓ-1]) * (x[iₓ]-mu[j]) / Δx)
                    #alternate likelihood calculation that only integrates some "width"(e.g. 3σ) of mu distribution...
# Otherwise, sum contributions from Gaussians at each point in distribution
        else
            likelihood = zero(float(eltype(dist)))
            @turbo for i ∈ 1:nₚ     # @turbo faster than @tturbo
# Likelihood curve follows a Gaussian PDF.
                likelihood += ( dist[i] / (sigma[j] * sqrt(2*π)) ) *
                        exp( - (x[i]-mu[j])*(x[i]-mu[j]) / (2*sigma[j]*sigma[j]) )
            end
            likelihood*=Δx
        end
        ll += log(likelihood)
# Normalize by total area under curve for intercomparability among proposals.
    end
    return ll
end

ll_param(x::Number,D::Nrm) = -(x-D.μ)*(x-D.μ)/(2*D.σ*D.σ)
ll_param(x::Number,D::lNrm) = -(x-D.μ)*(x-D.μ)/(2*D.σ*D.σ)
    #lnx = log(x); return -lnx-(lnx-D.μ)*(lnx-D.μ) / (2*D.σ*D.σ)
ll_param(x::T,D::Unf) where T<:Number = zero(T)

function ll_params(p::NamedTuple,d::NamedTuple)
    ll = 0.
    ll += ll_param(p.Cp, d.Cp)
    ll += ll_param(p.R, d.R)
    ll += ll_param(p.Tc, d.Tc)
    ll += ll_param(p.Tm, d.Tm)
    ll += ll_param(p.cAl, d.cAl)
    ll += ll_param(p.k, d.k)
    ll += ll_param(p.rAlo, d.rAlo)
    ll += ll_param(p.ta, d.ta)
    ll += ll_param(p.tss, d.tss)
    ll += ll_param(p.ρ, d.ρ)
    return ll
end

# Status function to keep user updated...
    function MetropolisStatus(p::NamedTuple,vars::Tuple,ll::Number,stepI::Integer,stepN::Integer,stage::String,t::Number;accpt::AbstractVector=[])
        println("---------------------------")
        stepI != 0 && println("Step $stepI of $stepN in $stage. \n")
        println("run time: ",round((time()-t)/60.,digits=2)," minutes \n")
        isempty(accpt) || println("acceptance rate =",reduce(+,accpt)/stepI)
        println("ll=$ll \n")
        for v ∈ vars
            println(v," → ",p[v])
        end
        println("---------------------------")
    end

function MetropolisAr(  time_domain::AbstractRange,
                        p::NamedTuple,   # Parameter proposal
                        pσ::NamedTuple, # proposed σ for pertrubations.
                        pvars::Tuple, # Variable parameters in proposal
                        mu::AbstractArray,  # Observed means
                        sigma::AbstractArray;# Observed 1σ's
                        plims::NamedTuple=(g=(),), # Paramter distributions.
                        burnin::Int=0,      # Burn-in iterations
                        nsteps::Int=10000,  # Post burn-in iterations
                        Δt::Number= 0.1,    # Time-step (Ma)
                        tmax::Number=2000,  # Max model duration (Ma, starts at CAIs)
                        Tmax::Number=1500,  # maximum temperature (K, solidus after 1200C max solidus in Johnson+2016)
                        Tmin::Number=0,     # minimum temperature (K)
                        nᵣ::Integer=100,    # Radial nodes
                        updateN::Integer=50) # Frequency of status updates (every `updateN` steps)
# Prepare output Distributions
    acceptanceDist = falses(nsteps)
    nᵥ = length(pvars)
    llDist = Array{float(eltype(mu))}(undef,nsteps) # Vector to track loglikelihood
    pDist = Array{float(eltype(mu))}(undef,nsteps,nᵥ) # Array to track proposal evolutions
    prt = similar(acceptanceDist,Symbol) # Vector to track proposed perturbations

# Calculate bincenters for time_domain
    Δd = step(time_domain) #calculate time_domain step
    bincenters = LinRange(first(time_domain) + 0.5Δd, last(time_domain) - 0.5Δd, length(time_domain)-1)

# If no plims given, set ranges to
    plims[1] == () && ( plims = (;zip(pvars,fill((-Inf,Inf,:U),length(pvars)))...) )

    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize acceptance probability at 50%
    stepfactor = 2.9

    # Sort the dataset from youngest to oldest
    sI = sortperm(mu)
    mu_sorted = mu[sI] # Sort means
    sigma_sorted = sigma[sI] # Sort uncertainty

    # These quantities will be used more than once
    datarows = length(mu_sorted)
    pₚ = p
    #p = copy(p)
    #pₚ = copy(p)
    #step_σ=copy(pσ)

    # Calculate initial proposal distribution
    dates,Vfrxn = PlntsmlAr(pₚ, Δt=Δt, tmax=tmax, nᵣ=nᵣ, Tmax=Tmax, Tmin=Tmin)
# Convert thermal code output into a binned histogram
    if iszero(pₚ.Fχ)
        distₚ = histogramify(time_domain,dates,Vfrxn)
    else
        Iages,Ivols = ImpactResetAr(dates,Vfrxn,pₚ,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
        distₚ = histogramify(time_domain,Iages,Ivols)
    end

# Log likelihood of initial proposal
    ll = llₚ = ll_dist(bincenters, distₚ, mu_sorted, sigma_sorted) + ll_params(p,plims)

# Start the clock
    start = time()

    # Burnin
    #@inbounds
    for i = 1:burnin
# Start with fresh slate of parameters
        #copyto!(pₚ, p)

# Adjust one parameter
        k = rand(pvars)
        δpₖ = pσ[k] * randn() #getproperty(step_σ,k)*randn()
        pₚ = perturb(p,k,p[k]+δpₖ)    #setproperty!(pₚ,k,getproperty(pₚ,k)+δpₖ)

# Calculate log likelihood for new proposal, ensuring bounds are not exceeded
        if !isa(plims[k], Unf) || plims[k].a < getproperty(pₚ,k) < plims[k].b
# Calculate cooling history if  pₚ[k] ∈ ( plims[k][1] , plims[k][2] )
            PlntsmlAr!(dates, Vfrxn, pₚ, Δt=Δt, tmax=tmax, nᵣ=nᵣ, Tmax=Tmax, Tmin=Tmin)
            #k == problem && println("problem"); flush(stdout)

# If >10% of interior radius melts, reject proposal
            if isnan(dates[div(nᵣ,10)])
                printstyled("meltdown rejected\n"; color=:light_magenta);flush(stdout)
                fill!(distₚ,zero(eltype(distₚ)))
# Only calculate Impact Resetting if flux is nonzero
            elseif iszero(pₚ.Fχ)
                histogramify!(distₚ,time_domain,dates,Vfrxn)
            else
                Iages,Ivols = ImpactResetAr(dates,Vfrxn,pₚ,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
                histogramify!(distₚ,time_domain,Iages,Ivols)
            end
# Ensure the returned distribution is nonzero
            if vreduce(+,distₚ) > 0 # actually faster than iszero() when there's lots of zeros
                llₚ = ll_dist(bincenters, distₚ , mu_sorted, sigma_sorted) + ll_params(pₚ,plims)
            else
                llₚ=-Inf
            end
# Reject proposal if propposal exceeds uniform bounds: pₚ[k] ∉ ( plims[k][1] , plims[k][2] )
        else
            llₚ = -Inf
        end

# Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
# Record new step sigma
            pσ = perturb(pσ,k,abs(δpₖ)*stepfactor) #setproperty!(step_σ,k,abs(δpₖ)*stepfactor)
# Record new parameters
            p = pₚ  #copyto!(p, pₚ)
# Record new log likelihood
            ll = llₚ
        end
        i%updateN == 0 && MetropolisStatus(p,pvars,ll,i,burnin,"Burn In",start); flush(stdout)
    end

# Hooray, we finished the burn-in, let's tell someone!
    println("===  BURN IN COMPLETE  ===\n\n")
    println("Post-Burn-In Status:")
    MetropolisStatus(p,pvars,ll,0,0,"",start)
    println("== == == == == == == == ==")
    println("Now beginning $nsteps Markov chain iterations...")
    flush(stdout)

    # Step through each of the N steps in the Markov chain
    #@inbounds
    for i=1:nsteps

        # Start with fresh slate of parameters
        #copyto!(pₚ, p)

        # Adjust one parameter
        k = rand(pvars)
        prt[i]=k
        δpₖ = pσ[k] * randn() #getproperty(step_σ,k)*randn()
        pₚ = perturb(p,k,p[k]+δpₖ)    #setproperty!(pₚ,k,getproperty(pₚ,k)+δpₖ)

# Calculate log likelihood for new proposal, ensuring bounds are not exceeded
        if !isa(plims[k], Unf) || plims[k].a < getproperty(pₚ,k) < plims[k].b
# Calculate cooling history if  pₚ[k] ∈ ( plims[k][1] , plims[k][2] )
            PlntsmlAr!(dates, Vfrxn, pₚ, Δt=Δt, tmax=tmax, nᵣ=nᵣ, Tmax=Tmax, Tmin=Tmin)
            #k == problem && println("problem"); flush(stdout)

# If >10% of interior radius melts, reject proposal
            if isnan(dates[div(nᵣ,10)])
                printstyled("meltdown rejected\n"; color=:light_magenta); flush(stdout)
                fill!(distₚ,zero(eltype(distₚ)))
# Only calculate Impact Resetting if flux is nonzero
            elseif iszero(pₚ.Fχ)
                histogramify!(distₚ,time_domain,dates,Vfrxn)
            else
                Iages,Ivols = ImpactResetAr(dates,Vfrxn,pₚ,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
                histogramify!(distₚ,time_domain,Iages,Ivols)
            end
# Ensure the returned distribution is nonzero
            if vreduce(+,distₚ) > 0
                llₚ = ll_dist(bincenters, distₚ , mu_sorted, sigma_sorted) + ll_params(pₚ,plims)
            else
                llₚ=-Inf
            end
# Reject proposal if propposal exceeds uniform bounds: pₚ[k] ∉ ( plims[k][1] , plims[k][2] )
        else
            llₚ = -Inf
        end

# Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
# Record new step sigma
            pσ = perturb(pσ,k,abs(δpₖ)*stepfactor) #setproperty!(step_σ,k,abs(δpₖ)*stepfactor)
# Record new parameters
            p = pₚ  #copyto!(p, pₚ)
# Record new log likelihood
            ll = llₚ
            acceptanceDist[i]=true
        end

        for j = 1:nᵥ
            pDist[i,j] = p[pvars[j]] #getproperty(p,pvars[j])
        end

        llDist[i] = ll
        i%updateN == 0 && MetropolisStatus(p,pvars,ll,i,nsteps,"Main Chain",start,accpt=acceptanceDist); flush(stdout)
    end
    MetOut = Dict{Symbol,Any}((pvars[i],pDist[:,i]) for i ∈ 1:length(pvars))
    for x ∈ keys(plims)
# Record proposal values of unvaried parameters
        in(x,pvars) || (MetOut[x]= vcat(getproperty(p,x),fill(NaN,nsteps-1)))
# Calculate exponent of lognormally distributed variables
        #isa(plims[x],lNrm) && vmapt!(exp,MetOut[x],MetOut[x])
    end
    MetOut[:ll] = llDist
    MetOut[:accept] = acceptanceDist
    MetOut[:prt] = prt
    return MetOut #(; MetOut...) # convert to NamedTuple
end
