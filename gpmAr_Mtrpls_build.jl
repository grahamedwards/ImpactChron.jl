## gpmAr Metropolis functions: construction site

"""
In production....
proposal =

        rAlo = p[2],     # initial solar ²⁶Al/²⁷Al
        tₐ = p[3],      # accretion time, My after CAIs
        #tₛₛ = p[1],       #solar system age, Ma
        #R = p[2]      # Body radius
        To = p[2],       # Disk temperature @ 2.5 au, K
        Al_conc = p[2],   # Fractional abundance of Al (g/g)
        Tc = p[2],       # Ar closure temperature, K
        ρ = p[2],       # rock density, kg/m³
        K = p[2],          # Thermal Conductivity
        Cₚ = p[2],     # Specific Heat Capacity

p=proposal[:,1]
pmax = p .+ proposal[:,2]
pmin = p .- proposal[:,2]
"""

## Log-likelihood calculation

"""
Log-likelihood notes...

Make sure everything is sorted.
    Use searchsortedfirst for speed.
Find corresponding date, interpolate, etc... and compare to each datum

"""
function ll_calc(   p_dist::Tuple,    # Proposed distribution (ages,proportion, radii)
                    mu::AbstractArray,      # 1D Array/Vector of observed μ's (sorted)
                    sigma::AbstractArray)   # 1D Array/Vector of observed σ's (sorted)

    # Define some frequently used variables
    ll = zero(float(eltype(p_dist[2])))
    nₘ = length(mu)           # n of "Measured" data
    nₚ = length(p_dist[1])         # N of "Proposed" distribution data

    # Sort relative to ages in x (this may be unnecessary)
    i_sort = sortperm(p_dist[1])
    x = p_dist[1][i_sort]
    dist = p_dist[2][i_sort]

    # Calculate integral of cooling age distribution with trapezoid rule
    ∫distdx = zero(float(eltype(p_dist[2])))
    for k = 1:length(x)-1
        ∫distdx += (x[k+1]-x[k]) * 0.5 * (dist[k+1] + dist[k])
    end
    # Convert values in discrete distribution to values
    dist ./= ∫distdx

    # Cycle through each datum in dataset
    # @inbounds
    for j = 1:nₘ
        # Find index of μ in the `dist` array
        iₓ = searchsortedfirst(x,mu[j]) # x[iₓ] ≥ mu[j]

        # If possible, prevent aliasing problems by interpolation
        if iₓ > 1 && (sigma[j] < abs(x[iₓ] - x[iₓ-1]))
            # Interpolate corresponding distribution value
            likelihood = dist[iₓ] - (x[iₓ]-mu[j]) * (dist[iₓ]-dist[iₓ-1])/(x[iₓ]-x[iₓ-1])

        # Otherwise, sum contributions from Gaussians at each point in distribution
        else
            likelihood = zero(float(eltype(dist)))
            # add @inbounds, then @turbo, then @tturbo.
            for i = 1:nₚ
                # Likelihood curve follows a Gaussian PDF.
                likelihood += dist[i] / (sigma[j] * sqrt(2*pi)) *
                        exp( - (x[i]-mu[j])*(x[i]-mu[j]) / (2*sigma[j]*sigma[j]) )
            end
        end
        ll += log(likelihood)
    end
    return ll
end



function MetropolisAr(  DistAr::Function,    # Proposal distribution calculator
                        x::AbstractArray,   # I think unneeded. timeseries from exp/gauss dist version
                        p::Proposal,   # Parameter proposal
                        pmin::Proposal,# Minimum parameter bounds
                        pmax::Proposal,# Maximum parameter bounds
                        pvars::Vector{Symbol}, # Variable parameters in proposal
                        mu::AbstractArray,  # Observed means
                        sigma::AbstractArray;# Observed 1σ's
                        burnin::Int=0,      # Burn-in iterations
                        nsteps::Int=10000,  # Post burn-in iterations
                        Δt::Number= 0.1,    # Time-step (Ma)
                        tmax::Number=2000,  # Max model duration (Ma, starts at CAIs)
                        nᵣ::Integer=100)    # Radial nodes
    # Prepare output Distributions
    acceptanceDist = falses(nsteps)

########### float(eltype(x))?????
    llDist = Array{float(eltype(x))}(undef,nsteps)
    pDist = Array{float(eltype(x))}(undef,nsteps,length(pvars))

    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize acceptance probability at 50%
    stepfactor = 2.9

    # Sort the dataset from youngest to oldest
    sI = sortperm(mu)
    mu_sorted = mu[sI] # Sort means
    sigma_sorted = sigma[sI] # Sort uncertainty

    # These quantities will be used more than once
    datarows = length(mu_sorted)
    p = copy(p)
    pₚ = copy(p)

    # Step sigma for Gaussian proposal distributions
    # Assumes a 1% σ step for the first Metropolis step.
    step_sigma = copy(p)./100

    # Log likelihood of initial proposal
    dist = similar(x)
    ll = llₚ = ll_calc(dist!(dist, x, p), mu_sorted, sigma_sorted, x[1], x[end])

############
# Change to appropriate
############
    PlntsmlAr(  tₛₛ = p.tss,       #solar system age, Ma
            rAlo = p.rAlo,     # initial solar ²⁶Al/²⁷Al
            tₐ = p.ta,      # accretion time, My after CAIs
            R = p.R,      # Body radius
            To = p.Tm,       # Disk temperature @ 2.5 au, K
            Al_conc = p.cAl,   # Fractional abundance of Al (g/g)
            Tc = p.Tc,       # Ar closure temperature, K
            ρ = p.ρ,       # rock density, kg/m³
            K = p.k,          # Thermal Conductivity
            Cₚ = p.Cp,     # Specific Heat Capacity
            Δt = Δt,      # absolute timestep, default 10 ka
            tmax = tmax,     # maximum time allowed to model
            nᵣ = nᵣ,        # radial nodes
            rmNaN=false) # allow NaNs to remain

    # Burnin
    for i=1:burnin
        # Start with fresh slate of parameters
        copyto!(pₚ, p)

        # Adjust one parameter
        k = rand(pvars)
        δpₖ = getproperty(step_sigma,k)*randn()
        setproperty!(pₚ,k,getproperty(pₚ,k)+δpₖ)
            #pₚ[k] += δpₖ

        # Calculate log likelihood for new proposal, ensuring bounds are not exceeded
        if getproperty(pmin,k) < getproperty(pₚ,k) < getproperty(pmax,k)
                #if  pmin[k] < pₚ[k] < pmax[k]
            llₚ = dist_ll(dist!(dist, x, pₚ), mu_sorted, sigma_sorted, x[1], x[end])
        else
            llₚ = -Inf # auto-reject proposal if bounds exceeded
        end

        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            # Record new step sigma
            setproperty!(step_sigma,k,abs(δpₖ)*stepfactor)
                #step_sigma[k] = abs(δpₖ)*stepfactor
            # Record new parameters
            copyto!(p, pₚ)
            # Record new log likelihood
            ll = llₚ
        end
    end
    # Step through each of the N steps in the Markov chain
    @inbounds for i=1:nsteps
        # Start with fresh slate of parameters
        copyto!(pₚ, p)

        # Adjust one parameter
        k = rand(1:length(p))
        δpₖ = step_sigma[k]*randn()
        pₚ[k] += δpₖ

        # Calculate log likelihood for new proposal, ensuring bounds are not exceeded
        if  pmin[k] < pₚ[k] < pmax[k]
            llₚ = dist_ll(dist!(dist, x, pₚ), mu_sorted, sigma_sorted, x[1], x[end])
        else
            llₚ = -Inf # auto-reject proposal if bounds exceeded
        end

        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            # Record new step sigma
            step_sigma[k] = abs(δpₖ)*stepfactor
            # Record new parameters
            copyto!(p, pₚ)
            # Record new log likelihood
            ll = llₚ
            acceptanceDist[i]=true
        end
        pDist[i,:] .= p
        llDist[i] = ll
    end
    return pDist, llDist, acceptanceDist
end
