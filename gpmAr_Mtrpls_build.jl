

"""
[x] Set input dist as tuple of age and abundnace vectors.
This will not use the xmin xmax frame.

Make sure everything is sorted.
Use searchsortedfirst for speed.
Find corresponding date, interpolate, etc... and compare to that datum


KEY part of condition is that sigma is < bin size.
I'm not sure why this means you DON'T have to compare to uncertainty.
But I get the idea that if uncertainty is small, it doesn't have to be fully flushed.

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
