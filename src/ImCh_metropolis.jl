## Metropolis Function

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

# The Metropolis algorithm applied to Ar-Ar measured thermal histories in meteorites
function MetropolisAr(  time_domain::AbstractRange,
                        p::NamedTuple,   # Parameter proposal
                        pσ::NamedTuple, # proposed σ for pertrubations.
                        pvars::Tuple, # Variable parameters in proposal
                        mu::AbstractArray,  # Observed means
                        sigma::AbstractArray, # Observed 1σ's
                        crater::NamedTuple; # crater/impact parameters
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
    plims[1] == () && ( plims = (;zip(pvars,fill(Unf(-Inf,Inf),length(pvars)))...) )

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
    dates,Vfrxn,radii,peakT = PlntsmlAr(pₚ, Δt=Δt, tmax=tmax, nᵣ=nᵣ, Tmax=Tmax, Tmin=Tmin)
# Convert thermal code output into a binned histogram
    if iszero(pₚ.Fχα) & iszero(pₚ.Fχβ)
        distₚ = histogramify(time_domain,dates,Vfrxn)
    else
        Iages,Ivols = ImpactResetAr(dates,Vfrxn,pₚ,crater,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
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
            PlntsmlAr!(dates, Vfrxn, peakT, pₚ, Δt=Δt, tmax=tmax, nᵣ=nᵣ, Tmax=Tmax, Tmin=Tmin)
            #k == problem && println("problem"); flush(stdout)

# If >10% of interior radius melts, reject proposal
            if isnan(dates[div(nᵣ,10)])
                printstyled("meltdown rejected\n"; color=:light_magenta);flush(stdout)
                fill!(distₚ,zero(eltype(distₚ)))
# Only calculate Impact Resetting if flux is nonzero
            elseif iszero(pₚ.Fχα) & iszero(pₚ.Fχβ)
                histogramify!(distₚ,time_domain,dates,Vfrxn)
            else
                Iages,Ivols = ImpactResetAr(dates,Vfrxn,pₚ,crater,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
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
        iszero(i%updateN) && ImpactChron.MetropolisStatus(p,pvars,ll,i,burnin,"Burn In",start); flush(stdout)
    end

# Hooray, we finished the burn-in, let's tell someone!
    println("===  BURN IN COMPLETE  ===\n\n")
    println("Post-Burn-In Status:")
    ImpactChron.MetropolisStatus(p,pvars,ll,0,0,"",start)
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
            PlntsmlAr!(dates, Vfrxn, peakT, pₚ, Δt=Δt, tmax=tmax, nᵣ=nᵣ, Tmax=Tmax, Tmin=Tmin)
            #k == problem && println("problem"); flush(stdout)

# If >10% of interior radius melts, reject proposal
            if isnan(dates[div(nᵣ,10)])
                printstyled("meltdown rejected\n"; color=:light_magenta); flush(stdout)
                fill!(distₚ,zero(eltype(distₚ)))
# Only calculate Impact Resetting if flux is nonzero
            elseif iszero(pₚ.Fχα) & iszero(pₚ.Fχβ)
                histogramify!(distₚ,time_domain,dates,Vfrxn)
            else
                Iages,Ivols = ImpactResetAr(dates,Vfrxn,pₚ,crater,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
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
        iszero(i%updateN) && ImpactChron.MetropolisStatus(p,pvars,ll,i,nsteps,"Main Chain",start,accpt=acceptanceDist); flush(stdout)
    end
    MetOut = Dict{Symbol,Any}((pvars[i],pDist[:,i]) for i ∈ 1:length(pvars))
    for x ∈ keys(plims)
# Record proposal values of unvaried parameters
        in(x,pvars) || (MetOut[x]= p[x])
# Calculate exponent of lognormally distributed variables
        #isa(plims[x],lNrm) && vmapt!(exp,MetOut[x],MetOut[x])
    end
    MetOut[:ll] = llDist
    MetOut[:accept] = acceptanceDist
    MetOut[:prt] = prt
    return MetOut
end
