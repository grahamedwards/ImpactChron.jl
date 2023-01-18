## Metropolis Function

# Status function to keep user updated...
function metropolis_status(p::NamedTuple,vars::Tuple,ll::Number,stepI::Integer,stepN::Integer,stage::String,t::Number;accpt::AbstractVector=[])
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
function thermochron_metropolis(  p::NamedTuple,   # Parameter proposal
                        pσ::NamedTuple, # proposed σ for pertrubations.
                        pvars::Tuple, # Variable parameters in proposal
                        mu::AbstractArray,  # Observed means
                        sigma::AbstractArray, # Observed 1σ's
                        impactsite::ImpactSite; # crater/impact parameters
                        plims::NamedTuple=(;), # Paramter distributions.
                        burnin::Int=0,      # Burn-in iterations
                        nsteps::Int,  # Post burn-in iterations
                        Δt::Number= 1.,    # Time-step (Ma)
                        tmax::Number=2000.,  # Max model duration (Ma, starts at CAIs)
                        Tmax::Number=1500.,  # maximum temperature (K, solidus after 1200C max solidus in Johnson+2016)
                        Tmin::Number=0.,     # minimum temperature (K)
                        nᵣ::Integer=100,    # Radial nodes
                        updateN::Integer=1_000, # Frequency of status updates (every `updateN` steps)
                        archiveN::Integer=0, # Save archive of output data every `archiveN` steps. Off (=0) by default.
                        downscale::Integer=1, # Downscale high-res timesteps to `downscale`-times fewer bins
                        petrotypes::PetroTypes=PetroTypes()) # petrologic types, each with max Temp and rel. abundances in record

# PREPARE OUTPUT DISTRIBUTIONS
    acceptanceDist = falses(nsteps)
    nᵥ = length(pvars)
    llDist = Array{float(eltype(mu))}(undef,nsteps) # Vector to track loglikelihood
    pDist = Array{float(eltype(mu))}(undef,nsteps,nᵥ) # Array to track proposal evolutions
    prt = similar(acceptanceDist,Symbol) # Vector to track proposed perturbations

# Prepare AsteroidHistory 

ah = AsteroidHistory(p.R, nnodes=nᵣ, Δt=Δt, tmax=tmax, downscale_factor=downscale)

# TIME MANAGEMENT
    # Ensure that age of CAIs (tₛₛ) is constant
    :tss ∈ pvars && error("tₛₛ must be a constant (:tss ∉ pvars) for time array framework to function properly")

# CHECK BOUNDS OF PROPOSAL VARIABLES
    # If no plims given, set infinite ranges to explore the studio space.
    isempty(plims) && ( plims = (;zip(pvars,fill(Unf(-Inf,Inf),length(pvars)))...) )
    # Check that uniformly distributed proposals do not violate their bounds.
    for i ∈ keys(plims)
        isa(plims[i],Unf) && ( plims[i].a <= p[i] <= plims[i].b || error("Initial proposal for $i exceeds permissible bounds ($(plims[i].a),$(plims[i].b))") )
    end


    stepfactor = 2.9 # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize acceptance probability at 50%

# SORT OBSERVED COOLING AGE DATASET (for ll_dist function)
    sI = sortperm(mu)
    mu_sorted = p.tss .- mu[sI] # Sort means as dates in Ma after CAIs
    sigma_sorted = sigma[sI] # Sort uncertainty

# Calculate initial proposal distribution
    pₚ = p # Use the "perturbed" version of `p`, pₚ, for consistancy.

# Calculate asteroid thermochronologic history
    asteroid_agedist!(ah, pₚ, petrotypes, impactsite,nᵣ=nᵣ, Tmax=Tmax,Tmin=Tmin) 
# Log likelihood of initial proposal
    ll = llₚ = ll_dist_params(ah,pₚ, plims,mu_sorted,sigma_sorted)
    
# Start the clock
    start = time()

# ~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~

## Burnin
    @inbounds for i = 1:burnin

# Adjust one parameter
        k = rand(pvars)
        δpₖ = pσ[k] * randn() #getproperty(step_σ,k)*randn()
        pₚ = perturb(p,k,p[k]+δpₖ)    #setproperty!(pₚ,k,getproperty(pₚ,k)+δpₖ)

# Calculate log likelihood for new proposal, ensuring bounds are not exceeded
        if isa(plims[k], Unf) && !(plims[k].a < pₚ[k] < plims[k].b)
# Reject proposal if propposal exceeds uniform bounds: pₚ[k] ∉ ( plims[k].a , plims[k].b )
            llₚ = -Inf
        else
# Calculate cooling history if  pₚ[k] ∈ ( plims[k][1] , plims[k][2] )
            asteroid_agedist!(ah, pₚ, petrotypes, impactsite,nᵣ=nᵣ, Tmax=Tmax,Tmin=Tmin) 
# Log likelihood of initial proposal
            llₚ = ll_dist_params(ah,pₚ, plims,mu_sorted,sigma_sorted)
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
        iszero(i%updateN) && ImpactChron.metropolis_status(p,pvars,ll,i,burnin,"Burn In",start); flush(stdout)
    end

# Hooray, we finished the burn-in, let's tell someone!
    println("===  BURN IN COMPLETE  ===\n\n")
    println("Post-Burn-In Status:")
    ImpactChron.metropolis_status(p,pvars,ll,0,0,"",start)
    println("== == == == == == == == ==")
    println("Now beginning $nsteps Markov chain iterations...")
    flush(stdout)

# ~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~

# Step through each of the nsteps in the Markov chain
    @inbounds for i=1:nsteps

        # Start with fresh slate of parameters
        #copyto!(pₚ, p)

        # Adjust one parameter
        k = rand(pvars)
        prt[i]=k
        δpₖ = pσ[k] * randn() #getproperty(step_σ,k)*randn()
        pₚ = perturb(p,k,p[k]+δpₖ)    #setproperty!(pₚ,k,getproperty(pₚ,k)+δpₖ)

# Calculate log likelihood for new proposal, ensuring bounds are not exceeded
        if isa(plims[k], Unf) && !(plims[k].a < pₚ[k] < plims[k].b)
# Reject proposal if propposal exceeds uniform bounds: pₚ[k] ∉ ( plims[k].a , plims[k].b )
            llₚ = -Inf
        else
# Calculate cooling history if  pₚ[k] ∈ ( plims[k][1] , plims[k][2] )
            asteroid_agedist!(ah, pₚ, petrotypes, impactsite,nᵣ=nᵣ, Tmax=Tmax,Tmin=Tmin) 
# Log likelihood of initial proposal
            llₚ = ll_dist_params(ah,pₚ, plims,mu_sorted,sigma_sorted)
        end

# Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
# Record new step sigma
            pσ = perturb(pσ,k,abs(δpₖ)*stepfactor) 
# Record new parameters
            p = pₚ  
# Record new log likelihood
            ll = llₚ
            acceptanceDist[i]=true
        end

        @inbounds for j = 1:nᵥ
            pDist[i,j] = p[pvars[j]] #getproperty(p,pvars[j])
        end

        llDist[i] = ll
        iszero(i%updateN) && ImpactChron.metropolis_status(p,pvars,ll,i,nsteps,"Main Chain",start,accpt=acceptanceDist); flush(stdout)
        iszero(archiveN) || iszero(i%archiveN) && Serialization.serialize("metropolis_archive_step_$i.js", (;acceptanceDist,llDist,pDist,prt) )
    end
    MetOut = Dict{Symbol,Any}((pvars[i],pDist[:,i]) for i ∈ eachindex(pvars))
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
