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
function ll_calc(   p_dist::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}},    # Proposed distribution (ages,proportion, radii)
                    mu::AbstractArray,      # 1D Array/Vector of observed μ's (sorted)
                    sigma::AbstractArray)   # 1D Array/Vector of observed σ's (sorted)

    # Define some frequently used variables
    ll = zero(float(eltype(p_dist[2])))
    nₘ = length(mu)           # n of "Measured" data
    nₚ = length(p_dist[1])         # N of "Proposed" distribution data
    tkr = zeros(nₘ)
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
        if iₓ > 1 && iₓ <= nₚ && (sigma[j] < abs(x[iₓ] - x[iₓ-1]))
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




## Testing


# Build main, max, min, σ `p` structs
prpsl = Proposal(4567.4,5.11e-5,1.5e5,2.13,0.011,250,550,3210,900,4)
prpsl_max = Proposal(4567.60 + 0.36,    # Jacobsen2008 CAI AJEF Pb-Pb upperbound
                    5.23e-5 + 0.13e-5, #Jacobsen2008 Allende CAI whole rock max
                    2.1e5,
                    2.3,
                    .012,
                    600.,
                    550 + 20,
                    3360.,
                    950.,       # Cp ~ max calculated for LLs EdwardsBlackburn2020
                    5.)

prpsl_min = Proposal(4567.4 - 0.34, # Jacobsen2008 AJEF-A34 Pb-Pb lowerbound
                    5.11e-5 - 0.14e-5, # Jacobsen2008 AJEF lowerbound
                    1.1e5,
                    1.8,
                    0.010,
                    0.,
                    550 - 20,
                    3160.,
                    750.,
                    3.)

prpsl_σ = Proposal(0.5 * 0.34,      # tss
                    0.5 * 0.14e-5,  # rAlo
                    1e3,            # R
                    0.1,            # ta
                    0.0001,          # cAl
                    50.,            # Tmidplane
                    5.,             # Tc
                    300.,           # ρ
                    25.,            # Cp
                    0.5)           # k

vars =[:ta,:cAl,:Tm,:Tc,:ρ,:Cp,:k]
# Includes: [:tss,:rAlo,:R,:ta,:cAl,:Tm,:Tc,:ρ,:Cp,:k]

d = PlntsmlAr(
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
            Δt = .01,        # absolute timestep, default 10 ka
            tmax = 2000.,    # maximum time allowed to model
            nᵣ = 100,        # radial nodes
            rmNaN=true)     # remove NaNs

ages = [4548.0, 4544.0, 4541.0, 4538.0, 4533.0, 4533.0, 4532.0, 4530.0, 4530.0, 4522.0, 4520.0, 4520.0, 4520.0, 4517.0, 4514.0, 4514.0, 4511.0, 4505.0, 4505.0, 4503.0, 4500.0, 4500.0, 4500.0, 4497.0, 4495.0, 4494.0, 4490.0, 4490.0, 4490.0, 4483.0, 4480.0, 4480.0, 4480.0, 4480.0, 4480.0, 4477.0, 4470.0, 4470.0, 4469.0, 4461.0, 4460.0, 4460.0, 4454.0, 4452.0, 4450.0, 4450.0, 4450.0, 4450.0, 4450.0, 4444.0, 4440.0, 4440.0, 4435.0, 4433.0, 4430.0, 4430.0, 4430.0, 4430.0, 4420.0, 4420.0, 4411.0, 4400.0, 4400.0, 4383.0, 4380.0, 4370.0, 4360.0, 4351.0, 4350.0, 4350.0, 4340.0, 4330.0, 4313.0, 4300.0, 4249.0, 4240.0, 4230.0, 4200.0, 4180.0, 4090.0, 4005.0, 3942.0, 3939.0, 3939.0, 3910.0, 3790.0, 3720.0, 3704.0, 3700.0, 3630.0, 3620.0, 3051.0]
uncert = [30.0, 18.0, 41.0, 13.0, 6.0, 8.0, 16.0, 20.0, 20.0, 8.0, 80.0, 10.0, 30.0, 11.0, 48.0, 20.0, 11.0, 10.0, 15.0, 52.0, 30.0, 30.0, 2.0, 9.0, 11.0, 46.0, 70.0, 30.0, 20.0, 14.0, 30.0, 30.0, 30.0, 8.0, 14.0, 20.0, 30.0, 20.0, 6.0, 8.0, 20.0, 10.0, 6.0, 9.0, 50.0, 30.0, 30.0, 30.0, 30.0, 17.0, 30.0, 40.0, 5.0, 4.0, 30.0, 30.0, 10.0, 40.0, 30.0, 20.0, 5.0, 30.0, 30.0, 10.0, 20.0, 10.0, 120.0, 8.0, 13.0, 10.0, 20.0, 40.0, 14.0, 70.0, 13.0, 20.0, 30.0, 50.0, 60.0, 40.0, 80.0, 23.0, 62.0, 62.0, 70.0, 40.0, 10.0, 35.0, 0.0, 10.0, 10.0, 8.0]

ll_calc(d,ages[1:50],uncert[1:50])


out_pDist , out_llDist , out_accept =
    MetropolisAr(   PlntsmlAr,    # Proposal distribution calculator
                prpsl,   # Parameter proposal
                prpsl_σ, # σ for Gauss. proposal distributions
                prpsl_min,# Minimum parameter bounds
                prpsl_max,# Maximum parameter bounds
                vars, # Variable parameters in proposal
                ages[1:80],  # Observed means
                uncert[1:80],# Observed 1σ's
                burnin=0,      # Burn-in iterations
                nsteps=1000,  # Post burn-in iterations
                Δt= 0.01,    # Time-step (Ma)
                tmax=2000,  # Max model duration (Ma, starts at CAIs)
                nᵣ=200)    # Radial nodes
