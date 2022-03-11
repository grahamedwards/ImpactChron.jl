## gpmAr Metropolis functions: construction site


## Log-likelihood calculation

"""
Log-likelihood notes for building from Chron.jl code...

[x] Set input dist as tuple of age and abundnace vectors.
This will not use the xmin xmax frame.

Make sure everything is sorted.
    Use searchsortedfirst for speed.
Find corresponding date, interpolate, etc... and compare to each datum

Interpolation conditional:
KEY component is that sigma is < bin size.
I'm not sure why this means you DON'T have to compare to uncertainty.
But I get the idea that if uncertainty is small, it doesn't have to be fully evaluated

"""
function ll_calc(   p_dist::Tuple{Vector,Vector,Vector},    # Proposed distribution (ages,proportion, ?radii?)
                    mu::AbstractArray,      # 1D Array/Vector of observed μ's (sorted)
                    sigma::AbstractArray)   # 1D Array/Vector of observed σ's (sorted)

    # Define some frequently used variables
    ll = zero(float(eltype(p_dist[2])))
    nₘ = length(mu)           # n of "Measured" data
    nₚ = length(p_dist[1])         # N of "Proposed" distribution data
    #dist_yave = mean(dist)          # mean value of distribution (presumably for normalization)
    #nbins = distrows - 1            # bins (or steps) in dist
    #dx = abs(xmax-xmin)             # Timespan of model

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
        if iₓ > 1 && (sigma[j] < abs(x[iₓ] - x[iₓ-1])) && iₓ < nₚ
            # Interpolate corresponding distribution value
            likelihood = dist[iₓ] - (x[iₓ]-mu[j]) * (dist[iₓ]-dist[iₓ-1])/(x[iₓ]-x[iₓ-1])

        # Otherwise, sum contributions from Gaussians at each point in distribution
        else
            likelihood = zero(float(eltype(dist)))
            # add @inbounds, then @turbo, then @tturbo.
            for i = 1:nₚ
                # Likelihood curve follows a Gaussian PDF. Note: dx cancels
                likelihood += dist[i] / (sigma[j] * sqrt(2*pi)) *
                        exp( - (x[i]-mu[j])*(x[i]-mu[j]) / (2*sigma[j]*sigma[j]) )
            end
        end
        ll += log(likelihood)
    end
    return ll
end




function MetropolisAr(dist!::Function, x::AbstractArray, p::AbstractArray, pmin::AbstractArray, pmax::AbstractArray, mu::AbstractArray, sigma::AbstractArray; burnin::Integer=0, nsteps::Int=10000)
    acceptanceDist = falses(nsteps)
    llDist = Array{float(eltype(x))}(undef,nsteps)
    pDist = Array{float(eltype(x))}(undef,nsteps,length(p))

    # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
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
    step_sigma = copy(p)./100

    # Log likelihood of initial proposal
    dist = similar(x)
    ll = llₚ = dist_ll(dist!(dist, x, p), mu_sorted, sigma_sorted, x[1], x[end])

    # Burnin
    for i=1:nsteps
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
