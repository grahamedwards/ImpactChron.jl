## Functions to aid visualization of Metropolis code outputs.
using Plots;gr()
using DelimitedFiles

## Load data from a csv

function csv2nt(filename::String;symbol::Bool=true)
    A = readdlm(filename,',') #Input Array
    kz = Tuple(Symbol.(A[1,:]))
    pN = Vector{Union{Vector{Float64},Vector{Symbol},Vector{Bool},Number}}(undef,length(kz)) #proto-NamedTuple Vector{Vector}
    for i = 1:length(kz)
        T = typeof(A[2,i])
        if T <: Number
            if isnan(A[3,i])
                pN[i] = float(A[2,i])
            else
                pN[i] = convert(Vector{T},A[2:end,i])
            end
        elseif T <: AbstractString
            symbol ? pN[i] = convert(Vector{Symbol},Symbol.(A[2:end,i])) : pN[i] = A[2:end,i]
        else
            k = kz[i]
            println("key $k did not match any criteria"); flush(stdout)
        end
    end
    return NamedTuple{kz}(pN)
end

## Plot Evolution of proposals
function plotproposals(d::Dict,plims::NamedTuple,cols::Integer;vars::Tuple=(),
                        ll::Bool=true,bounds::Bool=true)
    isempty(vars) ? v=keys(plims) : v=vars
    ll && (v=tuple(:ll,v...))
    nᵥ=length(v)

#Convert to strings if necessary
    isequal(eltype(keys(d)),String) ? (v= String.(v); llₛ = "ll" ; acpt = "accept") : (llₛ = :ll; acpt = :accept)
# Calculate number of rows needed to accomodate all variables in `cols` columns.
    rows = Int(ceil(nᵥ/cols,digits=0))

    panels = Vector{Any}(nothing,nᵥ)
    for i ∈ 1:nᵥ
        k = v[i]
        y = d[k]
        x = 1:length(y)
        if k == llₛ
            panels[i] = plot(x,y,xticks=[],ylabel="ll",linecolor=:black) #use \scrl eventually
            r = 100 * round(sum(d[acpt])/length(d[acpt]),digits=3)
            xlabel!("acceptance = $r %",xguidefontsize=6)
            #annotate!(last(x), (y[end]+y[1])/2, text("acceptance = $r %", :black,:bottomleft,6))
        elseif isnan(last(y))
            panels[i] = plot([1,last(x)],fill(y[1],2),xticks=[],ylabel="$k",linecolor=:black)
        else
            panels[i] = plot(x,y,xticks=[],linecolor=:black)
        end

        if bounds && k != llₛ
            B = plims[Symbol(k)]
            if isa(B,Unf)
                plot!([1,last(x)],fill(B.a,2),ylabel="$k",linecolor=:grey,linestyle=:solid)
                plot!([1,last(x)],fill(B.b,2),linecolor=:grey,linestyle=:solid)
            elseif isa(B,Nrm)
                plot!([1,last(x)],fill(B.μ+B.σ,2),ylabel="$k",linecolor=:grey,linestyle=:dash)
                plot!([1,last(x)],fill(B.μ-B.σ,2),linecolor=:grey,linestyle=:dash)
            elseif isa(B,lNrm)
                plot!([1,last(x)],fill(B.μ+B.σ,2),ylabel="log[" * "$k" * "]",linecolor=:grey,linestyle=:dashdot)
                plot!([1,last(x)],fill(B.μ-B.σ,2),linecolor=:grey,linestyle=:dashdot)
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


using KernelDensity, StatGeochem
function BootstrapKDE_xD(data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}; cutoff::Number=-0.05)
# Array to hold stacked, scaled data
    allscaled = Array{float(eltype(data)),1}([])
    agemin = agemax = zero(float(eltype(ages)))
# For each row of data
    for i=1:size(data,2)
        μ, σ = data[:,i], sigma[:,i]

# Maximum extent of expected analytical tail (beyond eruption/deposition/cutoff)
        maxTailLength = nanmean(σ) .* norm_quantile(1 - 1/(1+countnotnans(μ)))
        included = (μ .- nanminimum(μ)) .>= maxTailLength
        included .|= μ .> nanmedian(μ) # Don't exclude more than half (could only happen in underdispersed datasets)
        included .&= .!isnan.(μ) # Exclude NaNs

# Include and scale only those data not within the expected analytical tail
        if sum(included)>0
            agemin = minimum(data[included,i])
            agemax = maximum(data[included,i])
            scaled = data[included,i] .- minimum(data[included,i])
            if maximum(scaled) > 0
                scaled = scaled ./ maximum(scaled)
            end
            append!(allscaled, scaled)
        end
    end

# Calculate kernel density estimate, truncated at 0
    kd = kde(allscaled,npoints=2^7)
    t = kd.x .> cutoff
    return range(agemin,agemax,sum(t)), kd.density[t]
end
