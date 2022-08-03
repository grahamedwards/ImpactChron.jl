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
function MetropolisAr(  p::NamedTuple,   # Parameter proposal
                        pσ::NamedTuple, # proposed σ for pertrubations.
                        pvars::Tuple, # Variable parameters in proposal
                        mu::AbstractArray,  # Observed means
                        sigma::AbstractArray, # Observed 1σ's
                        crater::NamedTuple; # crater/impact parameters
                        plims::NamedTuple=(g=(),), # Paramter distributions.
                        burnin::Int=0,      # Burn-in iterations
                        nsteps::Int,  # Post burn-in iterations
                        Δt::Number= 1.,    # Time-step (Ma)
                        downscale::Integer=0, # Downscale high-res timesteps to `downscale`-times fewer bins
                        tmax::Number=2000.,  # Max model duration (Ma, starts at CAIs)
                        Tmax::Number=1500.,  # maximum temperature (K, solidus after 1200C max solidus in Johnson+2016)
                        Tmin::Number=0.,     # minimum temperature (K)
                        nᵣ::Integer=100,    # Radial nodes
                        updateN::Integer=1_000, # Frequency of status updates (every `updateN` steps)
                        archiveN::Integer=0) # Save archive of output data every `archiveN` steps. Off (=0) by default.

# Prepare output Distributions
    acceptanceDist = falses(nsteps)
    nᵥ = length(pvars)
    llDist = Array{float(eltype(mu))}(undef,nsteps) # Vector to track loglikelihood
    pDist = Array{float(eltype(mu))}(undef,nsteps,nᵥ) # Array to track proposal evolutions
    prt = similar(acceptanceDist,Symbol) # Vector to track proposed perturbations

# Time Management:
# Declare age variable: the comprehensive timeseries of the model solar system history.
    # First ensure that age of CAIs (tₛₛ) is constant
    :tss ∉ pvars || error("tₛₛ must be a constant (:tss ∉ pvars) for time array framework to function properly")

    if iszero(downscale) # If a downscale value is declared (≠0), make sure length of time_r and time_ll are divisible.
        time_r = time_ll = 0:Δt:tmax
    else
        if (0 >= p.Fχα) & (0 >= p.Fχβ); error("Downscaling without the array framework is not supported") end
        downscale_adj = length(0:Δt:tmax)%downscale
        tmax_adj = tmax-downscale_adj*Δt
        iszero(downscale_adj) || @warn "time range adjusted for downscale to 0:Δt(=$Δt):$tmax_adj"
        time_r = 0:Δt:tmax_adj
        time_ll = vmean(time_r[1:1+downscale]):Δt*downscale:tmax_adj
    end

    time_v = collect(time_r)
# Calculate bound values for age_range
    time_bounds = rangemidbounds(time_r)
    #age_bounds = rangemidbounds(age)
# Deal with statistical bounds of proposal variables
# If no plims given, set infinite ranges to explore the studio space.
    plims[1] == () && ( plims = (;zip(pvars,fill(Unf(-Inf,Inf),length(pvars)))...) )
# Check that all proposed parameters do not violate uniform distribution bounds
    for i ∈ keys(plims)
        isa(plims[i],Unf) && ( plims[i].a <= p[i] <= plims[i].b || error("Initial proposal for $i exceeds permissible bounds ($(plims[i].a),$(plims[i].b))") )
    end

# standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize acceptance probability at 50%
    stepfactor = 2.9

# Sort the dataset from youngest to oldest
    sI = sortperm(mu)
    mu_sorted = p.tss .- mu[sI] # Sort means as dates in Ma after CAIs
    sigma_sorted = sigma[sI] # Sort uncertainty

# These quantities will be used more than once
    tₓr = Array{eltype(mu_sorted)}(undef,length(time_v),nᵣ) # time x radial position array to be used in impact resetting scheme
    impacts = Vector{Float64}(undef,length(time_v)) # tracker of # of impacts at each timestep
    tcoolₒ = Vector{Int64}(undef,nᵣ) # tracker of indices of primary cooling date in age and time columns of tₓr
    distₚ = Vector{eltype(tₓr)}(undef,length(time_v))
# Calculate initial proposal distribution
    pₚ = p # Use the "perturbed" version of `p`, pₚ, for consistancy.
    dates,Vfrxn,radii,peakT = PlntsmlAr(pₚ, Δt=Δt, tmax=tmax, nᵣ=nᵣ, Tmax=Tmax, Tmin=Tmin)
# Only calculate impact resetting if flux is positive and nonzero
    if (0 >= pₚ.Fχα) & (0 >= pₚ.Fχβ)
        histogramify!(distₚ,time_bounds,dates,Vfrxn)
    else
        impact_reset_array!(tₓr, time_v, impacts, tcoolₒ, dates, Vfrxn, pₚ, crater, nᵣ=nᵣ,Δt=Δt)
        if iszero(downscale)
            distₚ .= vsum(tₓr,dims=2)
        else
            distₚ = Vector{eltype(tₓr)}(undef,length(time_ll))
            ImpactChron.downscale!(distₚ,vsum(tₓr,dims=2))
        end
    end
# Log likelihood of initial proposal
    ll = llₚ = ll_dist(time_ll, distₚ, mu_sorted, sigma_sorted) + ll_params(p,plims)

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
# Only calculate impact resetting if flux is positive and nonzero
            elseif (0 >= pₚ.Fχα) & (0 >= pₚ.Fχβ)
                histogramify!(distₚ,time_bounds,dates,Vfrxn)
            else
                impact_reset_array!(tₓr, time_v, impacts, tcoolₒ, dates, Vfrxn, pₚ, crater, nᵣ=nᵣ,Δt=Δt)
                if iszero(downscale)
                    distₚ .= vsum(tₓr,dims=2)
                else
                    ImpactChron.downscale!(distₚ,vsum(tₓr,dims=2))
                end
            end
# Ensure the returned distribution is nonzero
            if vreduce(+,distₚ) > 0 # actually faster than iszero() when there's lots of zeros
                llₚ = ll_dist(time_ll, distₚ , mu_sorted, sigma_sorted) + ll_params(pₚ,plims)
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
# Only calculate impact resetting if flux is positive and nonzero
            elseif (0 >= pₚ.Fχα) & (0 >= pₚ.Fχβ)
                histogramify!(distₚ,time_bounds,dates,Vfrxn)
            else
                impact_reset_array!(tₓr, time_v, impacts, tcoolₒ, dates, Vfrxn, pₚ, crater, nᵣ=nᵣ,Δt=Δt)
                if iszero(downscale)
                    distₚ .= vsum(tₓr,dims=2)
                else
                    ImpactChron.downscale!(distₚ,vsum(tₓr,dims=2))
                end
            end
# Ensure the returned distribution is nonzero
            if vreduce(+,distₚ) > 0
                llₚ = ll_dist(time_ll, distₚ, mu_sorted, sigma_sorted) + ll_params(pₚ,plims)
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
        iszero(archiveN) || iszero(i%archiveN) && Serialization.serialize("metropolis_archive_step_$i.js", (;acceptanceDist,llDist,pDist,prt) )
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
