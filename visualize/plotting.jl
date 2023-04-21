## Functions to produce visualizations of of ImpactChron outputs

if !isdefined(@__MODULE__,:interleave)
    @warn "Trying to load datamgmt.jl"
    try
        localdir= @__DIR__
        include("$localdir/datamgmt.jl")
    catch
        @warn "Load `datamgmt.jl` for functionality"
    end
end

using CairoMakie
import Distributions

"""

```julia
proposalhists_priordists(data_in::Dict,plims::NamedTuple,v::Tuple; cols=3, nbins=32,figsize=(800,600),darkmode=false)
```

Plots a grid of pretty (top-traced with translucent fill) histograms for parameters listed in `v` and Markov chains in `data_in`, overlain by prior distributions in `plims`.

`cols` specifies the number of columns in the histogram grid,
`nbins` specifies the number of bins to plot,
`figsize` specifies full figure size in pixels.

Set to a transparent background color scheme by setting `darkmode=true`

"""
function proposalhists_priordists(data_in::Dict,plims::NamedTuple,v::Tuple;cols::Int=3, nbins::Int=32,figsize=(800,600),darkmode::Bool=false)
    rows = ceil(Int,length(v)/cols)
    f=Figure(resolution=(figsize), backgroundcolor=ifelse(darkmode,:transparent,:white))
    for j in CartesianIndices((rows,cols))
        i = j[1] + rows * (j[2]-1) #calculate index in v
        if i <= length(v)
            proposalhist_priordist(v[i], data_in[v[i]], plims[v[i]], f=f[j[1],j[2]], nbins=nbins, darkmode=darkmode)
        end
    end
    f
end

function proposalhist_priordist(v::Symbol, data_in::Vector,B::ImpactChron.PriorDistribution; f=Figure(), nbins::Int=32,darkmode::Bool=false)

    names = Dict(   
        :tss => "Age of CAIs (Ma)", 
        :rAlo=>"Initial ²⁶Al/²⁷Al (×10⁻⁵)",
        :Tm=> "Midplane temperature (K, 2.5 AU)",
        :R => "Radius (km)",
        :ta=> "Accretion time (Ma after CAIs)", 
        :cAl=> "Al abundance (wt%)",
        :ρ=> "Bulk density (kg/m³)",
        :Cp=> "Specific heat Capacity (J/kg•K)",
        :k => "Thermal conductivity (W/m•K)",
        :Tc=> "Ar closure temperature (K)",
        :tχα=> "Primordial bombardment onset (Ma after CAIs)",
        :τχα=> "Primordial bombardment ℯ-folding time (Ma)",
        :Fχα=> "Primordial initial impactor flux (Ma⁻¹)",
        :tχβ=> "Post-accretion bombardment onset (Ma after CAIs)",
        :τχβ=> "Post-accretion bombardment ℯ-folding time (Ma)",
        :Fχβ=> "Post-accretion initial impactor flux (Ma⁻¹)",
        :tχγ=> "2ⁿᵈ Post-accretion bombardment onset (Ma after CAIs)",
        :τχγ=> "2ⁿᵈ Post-accretion bombardment ℯ-folding time (Ma)",
        :Fχγ=> "2ⁿᵈ Post-accretion initial impactor flux (Ma⁻¹)",            )

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
    if isa(B,Unf)
        prdst = Distributions.Uniform(B.a,ifelse(isinf(B.b),maximum(x),B.b))
    elseif isa(B,Nrm)
        prdst = Distributions.Normal(B.μ,B.σ)
    elseif isa(B,lNrm)
        prdst = Distributions.LogNormal(B.μ,B.σ)
    end
    prdst_x = LinRange(first(h.x),last(h.x),100)
    lines!(prdst_x,Distributions.pdf.(prdst,prdst_x),linewidth=2,linestyle=:dash,color=(fillcolor,0.6))
    f
end



## Plot Evolution of proposals

function plotproposal(v::Symbol,x::Vector;f=Figure(),B=(),linecolor=:black)
    lc = linecolor
    if v==:ll 
        yl = "ℓ"
    elseif B isa lNrm
        yl = "ln[" * "$v" * "]"
    else
        yl = "$v"
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

function plotproposals(data_in::Dict,plims::NamedTuple,v::Tuple;ll::Bool=true, cols::Int=3, figsize=(800,600),linecolor=:black)
    llind=0
    lv=length(v) + ifelse(ll,1,0)
    rows = ceil(Int,lv/cols)
    f=Figure(resolution=(figsize))
    for j in CartesianIndices((rows,cols))
        i = j[1] + rows * (j[2]-1) #calculate index in v
        if i <= length(v)
            plotproposal(v[i], data_in[v[i]],B=plims[v[i]], f=f[j[1],j[2]],linecolor=linecolor)
        elseif i==lv
            llind = j
        end
    end
    ll && plotproposal(:ll, data_in[:ll], f=f[llind[1],llind[2]],linecolor=linecolor)
    ar= round(100vmean(data_in[:accept]), digits=1)
    #text!("acceptance rate = $ar %",position=(length(data_in[:ll]), minimum(data_in[:ll])), align=(:right,:top))
    #Label(f[end+1,:],"acceptance rate = $ar %", halign = :center)
    f
end

"~Plotting Functions Loaded Successfully~"