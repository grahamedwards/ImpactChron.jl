

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
function ll_calc(    p_dist::Tuple{Vector,Vector,Vector},    # Proposed distribution (ages,proportion, ?radii?)
                    mu::AbstractArray,      # 1D Array/Vector of observed μ's (sorted)
                    sigma::AbstractArray)   # 1D Array/Vector of observed σ's (sorted)

    # Define some frequently used variables
    ll = zero(float(eltype(p_dist[2])))
    nₘ = length(mu)           # n of "Measured" data
    nₚ = length(p_dist[1])         # N of "Proposed" distribution data
    #dist_yave = mean(dist)          # mean value of distribution (presumably for normalization)
    #nbins = distrows - 1            # bins (or steps) in dist
    #dx = abs(xmax-xmin)             # Timespan of model

    # Remove NaN's from model output
    Xnan = .!isnan.(p_dist[1])
    x_unsorted = p_dist[1][Xnan]
    dist_unsorted = p_dist[2][Xnan]

    # Sort relative to ages in x (this may be unnecessary)
    i_sort = sortperm(x_unsorted)
    x = x_unsorted[i_sort]
    dist = dist_unsorted[i_sort]

    ∫distdx = diff(x) .*

    # Cycle through each datum in dataset
    # @inbounds
    for j = 1:nₘ
        # Find index of μ in the `dist` array
        iₓ = searchsortedfirst(x,mu[j]) # x[iₓ] ≥ mu[j]
        Uₓ = x[iₓ] # upper bound age in distribution

        # If possible, prevent aliasing problems by interpolation
        if iₓ > 1 && (sigma[j] < (Uₓ - x[iₓ-1])) && iₓ < nₚ
            # Interpolate corresponding distribution value
            Lₓ = x[iₓ-1] # lower bound age in distribution

            likelihood = ( dist[iₓ]*(mu[j]-Lₓ) + dist[iₓ-1]*(Uₓ-mu[j]) ) /
                ( (last(x)-x[1]) * mean(dist) )
                # is this valid though?
                # I think it should be divided by x-distance.

        # Otherwise, sum contributions from Gaussians at each point in distribution
        else
            likelihood = zero(float(eltype(dist)))
            # add @inbounds, then @turbo, then @tturbo.
            for i = 2:nₚ-1
                # Likelihood curve follows a Gaussian PDF. Note: dx cancels
                likelihood += dist[i] * 0.5 * (x[i+1] - x[i-1])
                / (mean(dist) * nₚ * sigma[j] * sqrt(2*pi)) *
                        exp( - (x[i]-mu[j])*(x[i]-mu[j]) / (2*sigma[j]*sigma[j]) )
            end
        end
        ll += log(likelihood)
    end
    return ll
end
