## gpmAr Metropolis functions: construction site
using LoopVectorization
## Plot Evolution of proposals

function plotproposals(d::Dict,plims::NamedTuple,cols::Integer;vars::Tuple=(),
                        ll::Bool=true,bounds::Bool=true)
    isempty(vars) ? v=keys(plims) : v=vars
    ll && (v=tuple(:ll,v...))
    nᵥ=length(v)
# Calculate number of rows needed to accomodate all variables in `cols` columns.
    rows = Int(ceil(nᵥ/cols,digits=0))

    panels = Vector{Any}(nothing,nᵥ)
    for i ∈ 1:nᵥ
        k = v[i]
        y = d[k]
        x = 1:length(y)
        if k == :ll
            panels[i] = plot(x,y,xticks=[],ylabel="ll",linecolor=:black) #use \scrl eventually
            r = 100 * sum(d[:accept])/length(d[:accept])
            annotate!(last(x), (y[end]+y[1])/2, text("acceptance = $r %", :black,:right,6))
        elseif isnan(last(y))
            panels[i] = plot([1,last(x)],fill(y[1],2),xticks=[],ylabel="$k",linecolor=:black)
        else
            panels[i] = plot(x,y,xticks=[],ylabel="$k",linecolor=:black)
        end

        if bounds && k != :ll
            B = plims[k]
            if isa(B,Unf)
                plot!([1,last(x)],fill(B.a,2),linecolor=:grey,linestyle=:solid)
                plot!([1,last(x)],fill(B.b,2),linecolor=:grey,linestyle=:solid)
            elseif isa(B,Nrm)
                plot!([1,last(x)],fill(B.μ+B.σ,2),linecolor=:grey,linestyle=:dash)
                plot!([1,last(x)],fill(B.μ-B.σ,2),linecolor=:grey,linestyle=:dash)
            elseif isa(B,lNrm)
                plot!([1,last(x)],fill(exp(B.μ+B.σ),2),linecolor=:grey,linestyle=:dashdot)
                plot!([1,last(x)],fill(exp(B.μ-B.σ),2),linecolor=:grey,linestyle=:dashdot)
            end
        end
    end
    sbplts=rows*cols
    Δplts = sbplts-length(panels)
    if Δplts > 0
        blnkplt = plot(legend=false,grid=false,foreground_color_subplot=:white)
        [ push!(panels,blnkplt) for j ∈ 1:Δplts]
    end

    plot(panels...,layout=grid(rows,cols),labels="")
end

## Log-likelihood calculation

"""
Log-likelihood notes...

"""
function ll_dist(   p_dist::Tuple{AbstractVector,AbstractVector},    # Proposed distribution (ages,proportions)
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

ll_param(x::Number,D::Nrm) = -(x-D.μ)*(x-D.μ)/(2*D.σ*D.σ)
ll_param(x::Number,D::lNrm) = -(x-D.μ)*(x-D.μ)/(2*D.σ*D.σ)
    #lnx = log(x); return -lnx-(lnx-D.μ)*(lnx-D.μ) / (2*D.σ*D.σ)
ll_param(x::T,D::Unf) where T<:Number = zero(T)

function ll_params(p::Proposal,d::NamedTuple)
    ll = 0.
    ll += ll_param(p.Cp, d.Cp)
    ll += ll_param(p.R, d.R)
    ll += ll_param(p.Tc, d.Tc)
    ll += ll_param(p.Tm, d.Tm)
    ll += ll_param(p.cAl, d.cAl)
    ll += ll_param(p.k, d.k)
    ll += ll_param(p.rAlo, d.rAlo)
    ll += ll_param(p.ta, d.ta)
    ll += ll_param(p.tss, d.tss)
    ll += ll_param(p.ρ, d.ρ)
    return ll
end

# Status function to keep user updated...
    function MetropolisStatus(p::Proposal,vars::Tuple,ll::Number,stepI::Integer,stepN::Integer,stage::String,t::Number;accpt::AbstractVector=[])
        println("---------------------------")
        stepI != 0 && println("Step $stepI of $stepN in $stage. \n")
        println("run time: ",round((time()-t)/60.,digits=2)," minutes \n")
        isempty(accpt) || println("acceptance rate =",reduce(+,accpt)/stepI)
        println("ll=$ll \n")
        for v ∈ vars
            println(v," → ",getproperty(p,v))
        end
        println("---------------------------")
    end

function MetropolisAr(  time_domain::AbstractRange,
                        p::Proposal,   # Parameter proposal
                        pσ::Proposal, # proposed σ for pertrubations.
                        pvars::Tuple, # Variable parameters in proposal
                        mu::AbstractArray,  # Observed means
                        sigma::AbstractArray;# Observed 1σ's
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
    bincenters = LinRange(first(time_domain)+0.5*Δd,last(time_domain) - 0.5Δd,length(time_domain)-1)

# If no plims given, set ranges to
    plims[1] == () && ( plims = (;zip(pvars,fill((-Inf,Inf,:U),length(pvars)))...) )

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
    step_σ=copy(pσ)

    # Calculate initial proposal distribution
    dates,Vfrxn = PlntsmlAr(
                tₛₛ = pₚ.tss,     #solar system age, Ma
                rAlo = pₚ.rAlo,  # initial solar ²⁶Al/²⁷Al
                tₐ = pₚ.ta,      # accretion time, My after CAIs
                R = pₚ.R,        # Body radius
                To = exp(pₚ.Tm), # (lNrm) Disk temperature @ 2.5 au, K
                Al_conc = exp(pₚ.cAl), # (lNrm) Fractional abundance of Al (g/g)
                Tc = pₚ.Tc,      # Ar closure temperature, K
                ρ = exp(pₚ.ρ),   # (lNrm) rock density, kg/m³
                K = exp(pₚ.k),   # (lNrm) Thermal Conductivity
                Cₚ = pₚ.Cp,      # Specific Heat Capacity
                Δt = Δt,tmax=tmax,nᵣ=nᵣ,Tmax=Tmax,Tmin=Tmin)
# Convert thermal code output into a binned histogram
    if iszero(pₚ.Fχ)
        distₚ = histogramify(time_domain,dates,Vfrxn,Δd=Δd)
    else
        Iages,Ivols = ImpactResetAr(dates,Vfrxn,pₚ,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
        distₚ = histogramify(time_domain,Iages,Ivols,Δd=Δd)
    end

# Log likelihood of initial proposal
    ll = llₚ = ll_dist((bincenters,distₚ), mu_sorted, sigma_sorted) + ll_params(p,plims)

# Start the clock
    start = time()

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
        if !isa(plims[k], Unf) || plims[k].a < getproperty(pₚ,k) < plims[k].b
# Calculate cooling history if  pₚ[k] ∈ ( plims[k][1] , plims[k][2] )
            PlntsmlAr!(dates,Vfrxn,
                tₛₛ = pₚ.tss,     #solar system age, Ma
                rAlo = pₚ.rAlo,  # initial solar ²⁶Al/²⁷Al
                tₐ = pₚ.ta,      # accretion time, My after CAIs
                R = pₚ.R,        # Body radius
                To = exp(pₚ.Tm), # (lNrm) Disk temperature @ 2.5 au, K
                Al_conc = exp(pₚ.cAl), # (lNrm) Fractional abundance of Al (g/g)
                Tc = pₚ.Tc,      # Ar closure temperature, K
                ρ = exp(pₚ.ρ),   # (lNrm) rock density, kg/m³
                K = exp(pₚ.k),   # (lNrm) Thermal Conductivity
                Cₚ = pₚ.Cp,      # Specific Heat Capacity
                Δt = Δt,tmax=tmax,nᵣ=nᵣ,Tmax=Tmax,Tmin=Tmin)
            #k == problem && println("problem"); flush(stdout)

# If >10% of interior radius melts, reject proposal
            if isnan(dates[div(nᵣ,10)])
                printstyled("meltdown rejected\n"; color=:light_magenta);flush(stdout)
                fill!(distₚ,zero(eltype(distₚ)))
# Only calculate Impact Resetting if flux is nonzero
            elseif iszero(pₚ.Fχ)
                histogramify!(distₚ,time_domain,Δd,dates,Vfrxn)
            else
                Iages,Ivols = ImpactResetAr(dates,Vfrxn,pₚ,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
                histogramify!(distₚ,time_domain,Δd,Iages,Ivols)
            end
# Ensure the returned distribution is nonzero
            if vreduce(+,distₚ) > 0 # actually faster than iszero() when there's lots of zeros
                llₚ = ll_dist((bincenters,distₚ) , mu_sorted, sigma_sorted) + ll_params(pₚ,plims)
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
            setproperty!(step_σ,k,abs(δpₖ)*stepfactor)
                #step_sigma[k] = abs(δpₖ)*stepfactor
            # Record new parameters
            copyto!(p, pₚ)
            # Record new log likelihood
            ll = llₚ
        end
        i%updateN == 0 && MetropolisStatus(p,pvars,ll,i,burnin,"Burn In",start); flush(stdout)
    end

# Hooray, we finished the burn-in, let's tell someone!
    println("===  BURN IN COMPLETE  ===\n\n")
    println("Post-Burn-In Status:")
    MetropolisStatus(p,pvars,ll,0,0,"",start)
    println("== == == == == == == == ==")
    println("Now beginning $nsteps Markov chain iterations...")
    flush(stdout)

    # Step through each of the N steps in the Markov chain
    #@inbounds
    for i=1:nsteps

        # Start with fresh slate of parameters
        copyto!(pₚ, p)

        # Adjust one parameter
        k = rand(pvars)
        prt[i]=k
        δpₖ = getproperty(step_σ,k)*randn()
        setproperty!(pₚ,k,getproperty(pₚ,k)+δpₖ)
            #pₚ[k] += δpₖ

# Calculate log likelihood for new proposal, ensuring bounds are not exceeded
        if !isa(plims[k], Unf) || plims[k].a < getproperty(pₚ,k) < plims[k].b
# Calculate cooling history if  pₚ[k] ∈ ( plims[k][1] , plims[k][2] )
            PlntsmlAr!(dates,Vfrxn,
                tₛₛ = pₚ.tss,     #solar system age, Ma
                rAlo = pₚ.rAlo,  # initial solar ²⁶Al/²⁷Al
                tₐ = pₚ.ta,      # accretion time, My after CAIs
                R = pₚ.R,        # Body radius
                To = exp(pₚ.Tm), # (lNrm) Disk temperature @ 2.5 au, K
                Al_conc = exp(pₚ.cAl), # (lNrm) Fractional abundance of Al (g/g)
                Tc = pₚ.Tc,      # Ar closure temperature, K
                ρ = exp(pₚ.ρ),   # (lNrm) rock density, kg/m³
                K = exp(pₚ.k),   # (lNrm) Thermal Conductivity
                Cₚ = pₚ.Cp,      # Specific Heat Capacity
                Δt = Δt,tmax=tmax,nᵣ=nᵣ,Tmax=Tmax,Tmin=Tmin)
            #k == problem && println("problem"); flush(stdout)

# If >10% of interior radius melts, reject proposal
            if isnan(dates[div(nᵣ,10)])
                printstyled("meltdown rejected\n"; color=:light_magenta); flush(stdout)
                fill!(distₚ,zero(eltype(distₚ)))
# Only calculate Impact Resetting if flux is nonzero
            elseif iszero(pₚ.Fχ)
                histogramify!(distₚ,time_domain,Δd,dates,Vfrxn)
            else
                Iages,Ivols = ImpactResetAr(dates,Vfrxn,pₚ,Δt=Δt,tmax=tmax,nᵣ=nᵣ)
                histogramify!(distₚ,time_domain,Δd,Iages,Ivols)
            end
# Ensure the returned distribution is nonzero
            if vreduce(+,distₚ) > 0
                llₚ = ll_dist((bincenters,distₚ) , mu_sorted, sigma_sorted) + ll_params(pₚ,plims)
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
        i%updateN == 0 && MetropolisStatus(p,pvars,ll,i,nsteps,"Main Chain",start,accpt=acceptanceDist); flush(stdout)
    end
    MetOut = Dict{Symbol,Any}((pvars[i],pDist[:,i]) for i ∈ 1:length(pvars))
    for x ∈ keys(plims)
# Record proposal values of unvaried parameters
        in(x,pvars) || (MetOut[x]= vcat(getproperty(p,x),fill(NaN,nsteps-1)))
# Calculate exponent of lognormally distributed variables
        isa(plims[x],lNrm) && vmapt!(exp,MetOut[x],MetOut[x])
    end
    MetOut[:ll] = llDist
    MetOut[:accept] = acceptanceDist
    MetOut[:prt] = prt
    return MetOut #(; MetOut...) # convert to NamedTuple
end
