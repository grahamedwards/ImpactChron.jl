## Metropolis code and support functions
    # metropolis_status
    # prior_bounds
    # strict_priors
    # thermochron_metropolis

# Status function to keep user updated...
function metropolis_status(p::NamedTuple,vars::Tuple,ll::Number,stepI::Integer,stepN::Integer,stage::String,t::Number;accpt::AbstractVector=[])
    println("---------------------------")
    stepI != 0 && println("Step $stepI of $stepN in $stage. \n")
    println("run time: ",round((time()-t)/60.,digits=2)," minutes \n")
    isempty(accpt) || println("acceptance rate =",vreduce(+,accpt)/stepI)
    println("ll=$ll \n")
    for v ∈ vars
        println(v," → ",p[v])
    end
    println("---------------------------")
end

"""

```julia
ImpactChron.prior_bounds(x, p<:PriorDistribution)
```

Evaluates whether a proposal `x` falls within the permissible bounds of its prior `p <:` [`PriorDistribution`](@ref). 
Always returns `true` if p is of type `Nrm` or `lNrm`, tests if `x ∈ (p.a,p.b)` for `p::Unf`.

"""
prior_bounds(x::Number,p::Unf) = p.a < x < p.b
prior_bounds(x::Number,p::Nrm) = true
prior_bounds(x::Number,p::lNrm) = true

"""

```julia
ImpactChron.strict_priors(p, k, p_prior<:PriorDistribution)
```

Evaluates strict priors related to Uniform distributions and other rules for paramter proposal `p` (`::NamedTuple`) and perturbed variable `k` (`::Symbol`), where `p_prior` is a [`PriorDistribution`](@ref) of `p[k]`.

Returns `true` if all priors are satisfied. Returns `false` if any priors fail.

Currently includes:\n
1. Ensure bounds of uniform priors are not exceeded. See [`ImpactChron.prior_bounds`](@ref).\n
2. Ensure bombardment events α, β, γ begin in sequential order.\n
3. The ℯ-folding time of bombardment α must be longer than those of β and γ.

"""
function strict_priors(p::NamedTuple,k::Symbol,p_prior::PriorDistribution)

# Return false if p[k] falls outside of Uniform bounds
    bool = prior_bounds(p[k],p_prior)
# Return false if the bombardment onset times fall out of order.
    bool *= p.tχα <= ifelse(iszero(p.Fχβ),Inf,p.tχβ) # if β flux is off, always accept
    Bγ = ifelse(iszero(p.Fχγ),Inf,p.tχγ) # if γ flux is off, always accept
    bool *= p.tχα <= Bγ
    bool *= p.tχβ <= Bγ
# The τ of instability/scattering fluxes must be shorter than that of the background impactor flux.
    bool *= ifelse(iszero(p.Fχβ),-Inf,p.τχβ) <= p.τχα # accept if β flux is off
    bool *= ifelse(iszero(p.Fχγ),-Inf,p.τχγ) <= p.τχα # accept if γ flux is off

    bool
end


"""
```julia
function thermochron_metropolis(p, pσ, pvars, mu, sigma, impactsite; nsteps, plims, downscale, petrotypes, kwargs...)
```

Runs a Markov chain Monte Carlo (MCMC) routine that explores the parameter space of the variables in `pvars` (as type `Symbol`, see PARAMETERS table below), constrained by observed thermochronologic ages given in `mu` with corresponding 1σ uncertainties in `sigma`.
Requires an initial proposal `p` (`::NamedTuple`) for parameter values and a Gaussian jump size `pσ` (`::NamedTuple`), each containing all model parameters listed in the table below.
Note that many parameters are log-normally distributed and require a natural-log-space initial guess (see table).
`impactsite` (`::ImpactSite`) describes the shape of asteroidal volume reheated per impact. See [`ImpactSite`](@ref) for details.


Four `kwargs` are particularly important:

1. `nsteps ::Int` is the number of post-burn-in Markov chain steps.

2. `plims ::NamedTuple` -- The model relies on hierarchical priors, such that all non-bombardment variables (`χ` ∉ parameter name) are constrained by prior distributions. 
Field names are as in `p` and `pσ` and values are subtypes of [`PriorDistribution`](@ref). 

3. `downscale ::Integer` bins the modeled age distribution to smooth it. The distribution is "downscaled" by the factor `downcale`. It is set to `1` by default (off). I recommend a value of `10`. See [`ImpactChron.downscale!`](@ref) for details. 

4. `petrotypes ::PetroTypes` contains the weights and upperbound temperatures of different petrologic types. This weighting is turned off by default. See [`PetroTypes`](@ref) for details.

---

PARAMETERS:

| Description               | log?  | `field`|
| :------------------------ | :--:  | :----: |
| solar system age (Ma)     | no    | `tss`  |
| initial ²⁶Al/²⁷Al         | no    | `rAlo` |
| closure temperature (K)   | yes   | `Tc`   | 
| body radius (m)           | yes   | `R`    | 
| accretion date (Myₛₛ)      | yes   | `ta`   | 
| disk temperature (K)      | yes   | `Tm`   |
| Al abundance (g/g)        | yes   | `cAl`  |
| density (kg/m³)           | yes   | `ρ`    | 
| thermal diffusivity       | yes   | `k`    | 
| specific heat capacity    | yes   | `Cp`   |
| α bombardment onset (Myₛₛ) | no    | `tχα`  | 
| α initial flux (My⁻¹)     | no    | `Fχα`  | 
| α ℯ-folding time (My)     | no    | `τχα`  | 
| β bombardment onset (Myₛₛ) | no    | `tχβ`  | 
| β initial flux (My⁻¹)     | no    | `Fχβ`  | 
| β ℯ-folding time (My)     | no    | `τχβ`  | 
| γ bombardment onset (Myₛₛ) | no    | `tχγ`  | 
| γ initial flux (My⁻¹)     | no    | `Fχγ`  | 
| γ ℯ-folding time (My)     | no    | `τχγ`  | 

---

All the other `kwargs`...

| kwarg | Default | Description |
| :---- | :-----: | :------ |
| `burnin` | 0 | MCMC burn-in iterations |
| `Δt` | 1 | model timestep (My) | 
| `Tmax` | 1500 | maximum temperature (K) | 
| `Tmin` | 0 | minimum temperature (K) |
| `nᵣ` | 100 | number of radial nodes in simulated asteroid | 
| `stepfactor` | 2.9 | scales each accepted jump (`pσ *= stepfactor`, tuned to ~50% acceptance rate) | 
| `updateN` | 1000 | print a status update of Markov chain every `...N` steps | 
| `archiveN` | 0 | save an archive of Markov chain every `...N` steps (`0`-> off, see also `serial2dict`) |
|`archiveages` | false | `true` -> archives the downscaled age dist. for each step |
| `rng` | `Random.Xoshiro()` | pseudorandom number generator seed | 

"""
function thermochron_metropolis(  p::NamedTuple,   # Parameter proposal
                        pσ::NamedTuple, # proposed σ for pertrubations.
                        pvars::Tuple, # Variable parameters in proposal
                        mu::AbstractArray,  # Observed means
                        sigma::AbstractArray, # Observed 1σ's
                        impactsite::ImpactSite; # crater/impact parameters
                        plims::NamedTuple=(;), # Paramter prior distributions.
                        burnin::Int=0,      # Burn-in iterations
                        nsteps::Int,  # Post burn-in iterations
                        Δt::Number= 1.,    # Time-step (Ma)
                        tmax::Number=2000.,  # Max model duration (Ma, starts at CAIs)
                        Tmax::Number=1500.,  # maximum temperature (K, solidus after 1200C max solidus in Johnson+2016)
                        Tmin::Number=0.,     # minimum temperature (K)
                        nᵣ::Integer=100,    # Radial nodes
                        updateN::Integer=1_000, # Frequency of status updates (every `updateN` steps)
                        archiveN::Integer=0, # Save archive of output data every `archiveN` steps. Off (=0) by default.
                        archiveages::Bool=false, # Archive the thermal history at each timestep.
                        downscale::Integer=1, # Downscale high-res timesteps to `downscale`-times fewer bins
                        petrotypes::PetroTypes=PetroTypes(), # petrologic types, each with max Temp and rel. abundances in record
                        rng = Random.Xoshiro(), # Seed a specific random number generator
                        stepfactor = 2.9) # σ of the proposal function is stepfactor * last_σ; 2.9 -> ~50% acceptance...
# PREPARE OUTPUT DISTRIBUTIONS
    acceptanceDist = falses(nsteps)
    nᵥ = length(pvars)
    llDist = Array{float(eltype(mu))}(undef,nsteps) # Vector to track loglikelihood
    pDist = Array{float(eltype(mu))}(undef,nsteps,nᵥ) # Array to track proposal evolutions
    prt = similar(acceptanceDist,Symbol) # Vector to track proposed perturbations
# Prepare AsteroidHistory 

ah = AsteroidHistory(p.R, nnodes=nᵣ, Δt=Δt, tmax=tmax, downscale_factor=downscale)

# Preallocate cooling age archive if keeping track of it.
archiveages && (agearchive = fill(NaN,length(ah.agedist_downscaled),nsteps))

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

    strict_priors(p,:tss,plims[:tss]) || error("'Strict' priors are not met.  See docs on `ImpactChron.strict_priors` for requirements.")

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
        k = rand(rng,pvars)
        δpₖ = pσ[k] * randn(rng)
        pₚ = perturb(p,k,p[k]+δpₖ)

# Ensure "strict" priors are met.
        if strict_priors(pₚ,k,plims[k]) 
# Calculate thermochronologic history
            asteroid_agedist!(ah, pₚ, petrotypes, impactsite,nᵣ=nᵣ, Tmax=Tmax,Tmin=Tmin) 
# Calculate log-likelihood of thermochronologic history
            llₚ = ll_dist_params(ah,pₚ, plims,mu_sorted,sigma_sorted)
# Reject proposal if it fails `strict_priors` tests
        else
            llₚ = -Inf
        end

# Decide to accept or reject the proposal
        if log(rand(rng)) < (llₚ-ll)
# Record new step sigma
            pσ = perturb(pσ,k,abs(δpₖ)*stepfactor) #setproperty!(step_σ,k,abs(δpₖ)*stepfactor)
# Record new parameters
            p = pₚ  
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

# Adjust one parameter
        k = rand(rng,pvars)
        prt[i]=k
        δpₖ = pσ[k] * randn(rng) 
        pₚ = perturb(p,k,p[k]+δpₖ)

# Ensure "strict" priors are met.
        if strict_priors(pₚ,k,plims[k]) 
# Calculate thermochronologic history
            asteroid_agedist!(ah, pₚ, petrotypes, impactsite,nᵣ=nᵣ, Tmax=Tmax,Tmin=Tmin) 
# Calculate log-likelihood of thermochronologic history
            llₚ = ll_dist_params(ah,pₚ, plims,mu_sorted,sigma_sorted)
# Reject proposal if it fails `strict_priors` tests
        else
            llₚ = -Inf
        end

# Decide to accept or reject the proposal
        if log(rand(rng)) < (llₚ-ll)
# Record new step sigma
            pσ = perturb(pσ,k,abs(δpₖ)*stepfactor) 
# Record new parameters
            p = pₚ  
# Record new log likelihood
            ll = llₚ
            acceptanceDist[i]=true
            archiveages && (agearchive[:,i] .= ah.agedist_downscaled)
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
    archiveages && (MetOut[:ages] = agearchive)
    return MetOut
end
