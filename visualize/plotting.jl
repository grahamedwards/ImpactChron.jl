## Plotting: Functions to produce visualizations of of ImpactChron outputs

if !isdefined(@__MODULE__,:interleave)
    @warn "Trying to load datamgmt.jl"
    try
        ImpactChron.
    catch
        @warn "Load data management tools for plotting functionality"
    end
end

# try 
#     using CairoMakie
# catch
#     @warn "All plotting functions rely on CairoMakie.\n\nTo make plots, please install it by typing into the REPL:  ]add CairoMakie\n\n"
# end


"""

```julia
proposalhists_priordists(   data_in::Dict, plims::NamedTuple, v::Tuple; 
                            cols=3, nbins=32, figsize=(800,600),darkmode=false)
```

Plots a grid of pretty (top-traced with translucent fill) histograms for parameters listed in `v` and Markov chains in `data_in` (output of `thermochron_metropolis`), overlain by prior distributions in `plims`.

| kwarg | description |
| :---- | :---------- |
|`cols` | number of columns in the histogram grid |
|`nbins`| number of bins/histogram |
|`figsize` | specifies figure size (px) |
|`darkmode=true`| transparent, dark background color scheme |

"""
function proposalhists_priordists(data_in::Dict,plims::NamedTuple,v::Tuple;cols::Int=3, nbins::Int=32,figsize=(800,600),darkmode::Bool=false)
    rows = ceil(Int,length(v)/cols)
    f=Figure(size=(figsize), backgroundcolor=ifelse(darkmode,:transparent,:white))
    for j in CartesianIndices((rows,cols))
        i = j[1] + rows * (j[2]-1) #calculate index in v
        if i <= length(v)
            proposalhist_priordist(v[i], data_in[v[i]], plims[v[i]], f=f[j[1],j[2]], nbins=nbins, darkmode=darkmode)
        end
    end
    f
end

function proposalhist_priordist(v::Symbol, data_in::Vector,B::ImpactChron.PriorDistribution; f=Figure(), nbins::Int=32,darkmode::Bool=false)

    names = ParamTitles
    x=deepcopy(data_in)
# Convert

    if v==:R 
        x .+= log(1e-3) # m -> km
        B = lNrm(B.μ +log(1e-3),B.σ)
    elseif v== :cAl
        x .+= log(100) # Multiply by 100 -> wt%, but in natural-log-space
        B = lNrm(B.μ +log(100),B.σ)
    elseif v== :rAlo 
        x .*= 1e5
        B = Nrm(B.μ *1e5,B.σ*1e5)
    end
    pltclr = fillcolor=ifelse(darkmode,:white,:black)

    isa(B,lNrm) && map!(exp,x,x)

    h = cleanhist(x, nbins=nbins,scooch_nbins=1)
    ax = Axis(f[1,1], xlabel=names[v],bottomspinecolor=pltclr,xtickcolor=pltclr,xticklabelcolor=pltclr, xlabelcolor=pltclr,backgroundcolor=ifelse(darkmode,:transparent,:white),
    xgridvisible=false,ygridvisible=false,yticklabelsvisible=false,yticksvisible=false,rightspinevisible=false,leftspinevisible=false,topspinevisible=false,)
    Makie.band!(ax,h.x,h.y,zero(h.y), color=(fillcolor,0.1))
    Makie.lines!(ax,h.x,h.y, color=pltclr,linewdith=2,)

# Plot prior distributions
prdst_x = LinRange(first(h.x),last(h.x),100)
prdst_y = if isa(B,Unf)
        fill(1/(ifelse(isinf(B.b),maximum(x),B.b)-B.a), length(prdst_x))
    elseif isa(B,Nrm)
         normdens.(prdst_x,B.μ,B.σ)
    elseif isa(B,lNrm)
        lognormdens.(prdst_x,B.μ,B.σ)
    end
    lines!(prdst_x, prdst_y, linewidth=2,linestyle=:dash,color=(fillcolor,0.6))
    f
end



## Plot Markov chain of one parameter

function plotproposal(v::Symbol,x::Vector;f=Figure(),B=(),linecolor=:black)
    lc = linecolor
    if v==:ll 
        yl = "ℓ"
    elseif B isa lNrm
        yl = "ln(" * ParamVars[v] * ")"
    else
        yl = ParamVars[v]
    end

    ax = Axis(f[1,1], ylabel=yl,xgridvisible=false, ygridvisible=false,xticksvisible=false,xticklabelsvisible=false)
    lines!(ax,x, color=(lc,0.6))
    #text!(ax,yl,position=(0,maximum(x)), align=(:left,:top),color=lc,fontsize=18)
    if B isa Unf
        b = (low=B.a,high=B.b)
    elseif B isa Nrm || B isa lNrm
        b = (low=B.μ-B.σ, high=B.μ+B.σ)
    else
        b=()
    end
    isempty(b) || hlines!(ax,[b.low,b.high], color=lc,linewidth=2,linestyle=:dash)
    f
end

# Plot Markov chains of several parameters
function plotproposals(data_in::Dict,plims::NamedTuple,v::Tuple;ll::Bool=true, cols::Int=3, figsize=(800,600),linecolor=:black, acceptance::Bool=false)
    llind=0
    
    lv=length(v) + ifelse(ll,1,0)
    rows = ceil(Int,lv/cols)
    
    f=Figure(size=(figsize))

    for j in CartesianIndices((rows,cols))
        i = j[1] + rows * (j[2]-1) #calculate index in v
        if i <= length(v)
            plotproposal(v[i], data_in[v[i]],B=plims[v[i]], f=f[j[1],j[2]],linecolor=linecolor)
        elseif i==lv
            llind = j
        end
    end

## Extras 
    ll && plotproposal(:ll, data_in[:ll], f=f[llind[1],llind[2]],linecolor=linecolor)
    
    ar= round(100ImpactChron.vmean(data_in[:accept]), digits=1)
    println("acceptance rate = $ar %")
    #text!("acceptance rate = $ar %",position=(length(data_in[:ll]), minimum(data_in[:ll])), align=(:right,:top))
    acceptance && Label(f[end+1,:],"acceptance rate = $ar %", halign = :center)
    f
end

"~Plotting Functions Loaded Successfully~"