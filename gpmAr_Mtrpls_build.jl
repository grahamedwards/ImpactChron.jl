## gpmAr Metropolis functions: construction site

## Log-likelihood calculation

"""
Log-likelihood notes...

"""
function ll_calc(   p_dist::Tuple{AbstractVector,AbstractVector,AbstractVector},    # Proposed distribution (ages,proportion, radii)
                    mu::AbstractVector,      # 1D Array/Vector of observed μ's (sorted)
                    sigma::AbstractVector)   # 1D Array/Vector of observed σ's (sorted)

# Define some frequently used variables
    ll = zero(float(eltype(p_dist[2])))
    nₘ = length(mu)            # n of "Measured" data
    nₚ = length(p_dist[1])     # N of "Proposed" distribution data

    #tkr = zeros(nₘ)
# Sort relative to ages in x this saves a lot of extra abs() tests
    i_sort = sortperm(p_dist[1])
    x = p_dist[1][i_sort]
    dist = p_dist[2][i_sort]

# Calculate area under each dist[x] and sum to estimate integral ∫distdx.
    distdx = Vector{float(eltype(dist))}(undef,nₚ)
    distdx[1] = 0.5 * (x[2] - x[1]) * dist[1]
    distdx[end] = 0.5 * (last(x) - x[end-1]) * dist[end]

    ∫distdx = distdx[1] + distdx[end]

    if nₚ > 2
        @inbounds for k ∈ 2:nₚ-1
            distdx[k] = 0.5 * (x[k+1] - x[k-1]) * dist[k]
            ∫distdx += distdx[k]
        end
    end

# Cycle through each datum in (mu,sigma)
    @inbounds for j ∈ 1:nₘ
# Find index of μ in the `dist` array
        iₓ = searchsortedfirst(x,mu[j]) # x[iₓ] ≥ mu[j]

# If possible, prevent aliasing problems by interpolation
        if (iₓ>1) && (iₓ<=nₚ) && ( (2sigma[j]) < (x[iₓ]-mu[j]) ) && ( (2sigma[j])<(mu[j]-x[iₓ-1]) )
            # && (sigma[j] < (x[iₓ]-x[iₓ-1]) ) # original threshold, see notes.
# Interpolate corresponding distribution value, note: (x[iₓ]-x[iₓ-1]) cancels in second term
            likelihood = dist[iₓ] * (x[iₓ]-x[iₓ-1]) - (x[iₓ]-mu[j]) * (dist[iₓ]-dist[iₓ-1])
            #likelihood = 6sigma[j] * (dist[iₓ] - (x[iₓ]-mu[j]) * (dist[iₓ]-dist[iₓ-1]) / (x[iₓ]-x[iₓ-1]) )
                    #alternate likelihood calculation that only integrates "width" of mu distribution...
# Otherwise, sum contributions from Gaussians at each point in distribution
        else
            likelihood = zero(float(eltype(dist)))
            @turbo for i ∈ 1:nₚ     # @turbo faster than @tturbo
# Likelihood curve follows a Gaussian PDF.
                likelihood += ( distdx[i] / (sigma[j] * sqrt(2*π)) ) *
                        exp( - (x[i]-mu[j])*(x[i]-mu[j]) / (2*sigma[j]*sigma[j]) )
            end
        end
        ll += log(likelihood/∫distdx)
# Normalize by total area under curve for intercomparability among proposals.
    end
    return ll
end



function MetropolisAr(  DistAr::Function,    # Proposal distribution calculator
                        p::Proposal,   # Parameter proposal
                        step_σ::Proposal, # σ for Gauss. proposal distributions
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
    nᵥ = length(pvars)
########### used to be `float(eltype(x))``
    llDist = Array{Float64}(undef,nsteps)
    pDist = Array{Float64}(undef,nsteps,nᵥ)

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


    # Calculate initial proposal distribution
    dist = DistAr(
                tₛₛ = p.tss,     #solar system age, Ma
                rAlo = p.rAlo,  # initial solar ²⁶Al/²⁷Al
                tₐ = p.ta,      # accretion time, My after CAIs
                R = p.R,        # Body radius
                To = p.Tm,      # Disk temperature @ 2.5 au, K
                Al_conc = p.cAl,# Fractional abundance of Al (g/g)
                Tc = p.Tc,      # Ar closure temperature, K
                ρ = p.ρ,        # rock density, kg/m³
                K = p.k,        # Thermal Conductivity
                Cₚ = p.Cp,      # Specific Heat Capacity
                Δt = Δt,        # absolute timestep, default 10 ka
                tmax = tmax,    # maximum time allowed to model
                nᵣ = nᵣ,        # radial nodes
                rmNaN=true)     # remove NaNs

    # Log likelihood of initial proposal
    ll = llₚ = ll_calc(dist, mu_sorted, sigma_sorted)

    # Burnin
    #@inbounds
    for i = 1:burnin
        # Start with fresh slate of parameters
        copyto!(pₚ, p)

        # Adjust one parameter
        k = rand(pvars)
        δpₖ = getproperty(step_σ,k)*randn()
        setproperty!(pₚ,k,getproperty(pₚ,k)+δpₖ)
            #pₚ[k] += δpₖ

        # Calculate log likelihood for new proposal, ensuring bounds are not exceeded
        if getproperty(pmin,k) < getproperty(pₚ,k) < getproperty(pmax,k)
        # if  pmin[k] < pₚ[k] < pmax[k]
            distₚ = dist = DistAr(
                        tₛₛ = pₚ.tss,     #solar system age, Ma
                        rAlo = pₚ.rAlo,  # initial solar ²⁶Al/²⁷Al
                        tₐ = pₚ.ta,      # accretion time, My after CAIs
                        R = pₚ.R,        # Body radius
                        To = pₚ.Tm,      # Disk temperature @ 2.5 au, K
                        Al_conc = pₚ.cAl,# Fractional abundance of Al (g/g)
                        Tc = pₚ.Tc,      # Ar closure temperature, K
                        ρ = pₚ.ρ,        # rock density, kg/m³
                        K = pₚ.k,        # Thermal Conductivity
                        Cₚ = pₚ.Cp,      # Specific Heat Capacity
                        Δt = Δt,        # absolute timestep, default 10 ka
                        tmax = tmax,    # maximum time allowed to model
                        nᵣ = nᵣ,        # radial nodes
                        rmNaN=true)     # remove NaNs

            llₚ = ll_calc(distₚ , mu_sorted, sigma_sorted)
        else
            llₚ = -Inf # auto-reject proposal if bounds exceeded
        end

        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            # Record new step sigma
            setproperty!(step_σ,k,abs(δpₖ)*stepfactor)
                #step_sigma[k] = abs(δpₖ)*stepfactor
            # Record new parameters
            copyto!(p, pₚ)
            # Record new log likelihood
            ll = llₚ
        end
    end

    # Step through each of the N steps in the Markov chain
    #@inbounds
    for i=1:nsteps

        # Start with fresh slate of parameters
        copyto!(pₚ, p)

        # Adjust one parameter
        k = rand(pvars)
        δpₖ = getproperty(step_σ,k)*randn()
        setproperty!(pₚ,k,getproperty(pₚ,k)+δpₖ)
            #pₚ[k] += δpₖ

        # Calculate log likelihood for new proposal, ensuring bounds are not exceeded
        if getproperty(pmin,k) < getproperty(pₚ,k) < getproperty(pmax,k)
        # if  pmin[k] < pₚ[k] < pmax[k]
            distₚ = dist = DistAr(
                        tₛₛ = pₚ.tss,     #solar system age, Ma
                        rAlo = pₚ.rAlo,  # initial solar ²⁶Al/²⁷Al
                        tₐ = pₚ.ta,      # accretion time, My after CAIs
                        R = pₚ.R,        # Body radius
                        To = pₚ.Tm,      # Disk temperature @ 2.5 au, K
                        Al_conc = pₚ.cAl,# Fractional abundance of Al (g/g)
                        Tc = pₚ.Tc,      # Ar closure temperature, K
                        ρ = pₚ.ρ,        # rock density, kg/m³
                        K = pₚ.k,        # Thermal Conductivity
                        Cₚ = pₚ.Cp,      # Specific Heat Capacity
                        Δt = Δt,        # absolute timestep, default 10 ka
                        tmax = tmax,    # maximum time allowed to model
                        nᵣ = nᵣ,        # radial nodes
                        rmNaN=true)     # remove NaNs

            llₚ = ll_calc(distₚ , mu_sorted, sigma_sorted)
        else
            llₚ = -Inf # auto-reject proposal if bounds exceeded
        end

        # Decide to accept or reject the proposal
        if log(rand()) < (llₚ-ll)
            # Record new step sigma
            setproperty!(step_σ,k,abs(δpₖ)*stepfactor)
                #step_sigma[k] = abs(δpₖ)*stepfactor
            # Record new parameters
            copyto!(p, pₚ)
            # Record new log likelihood
            ll = llₚ
            acceptanceDist[i]=true
        end

        for j = 1:nᵥ
            pDist[i,j] = getproperty(p,pvars[j])
        end

        llDist[i] = ll
    end
    return pDist, llDist, acceptanceDist
end
