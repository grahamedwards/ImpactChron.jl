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

import Plots; Plots.gr()
import Distributions

"""

```julia
proposalhists_priordists(data_in::Dict,plims::NamedTuple,v::Tuple;
    cols::Integer=1, nbins::Integer=20,figsize=(800,600),darkmode::Bool=false)
```

Plots a grid of pretty (top-traced with translucent fill) histograms for parameters listed in `v`.
Parameter `plims` tracks whether log-normal distributions need to be converted back to linear space.

`cols` specifies the number of columns in the histogram grid,
`nbins` specifies the number of bins to plot,
`figsize` specifies full figure size in pixels.

Set to a transparent background color scheme by setting `darkmode=true`

"""
function proposalhists_priordists(data_in::Dict,plims::NamedTuple,v::Tuple;
    cols::Integer=1, nbins::Integer=20,figsize=(800,600),darkmode::Bool=false)

    names = Dict(   
        :tss => "Age of CAIs (Ma)", 
        :rAlo=>"Initial ²⁶Al/²⁷Al",
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

    nᵥ=length(v)

    d=deepcopy(data_in)
# Convert
    d[:R] .+= log(1e-3) # m -> km
    d[:cAl] .+= log(100) # Multiply by 100 -> wt%, but in natural-log-space
# INCORPORATE ADJUSTMENTS INTO BOUNDS BELOW

#Convert to strings if necessary
    isequal(eltype(keys(d)),String) ? (v= String.(v); acpt = "accept") : ( acpt = :accept)
# Calculate number of rows needed to accomodate all variables in `cols` columns.
    rows = ceil(Int,nᵥ/cols)

    pltclr,bkgrnd = ifelse(darkmode,(:white,:transparent),(:black,:white))


    panels = Vector{Any}(nothing,nᵥ)
    for i ∈ 1:nᵥ
        k = v[i]
        x = d[k]

        isa(plims[k],lNrm) && map!(exp,x,x)

        x_scooch = (maximum(x)-minimum(x))/ (nbins-4)
        binedges = LinRange(minimum(x)-2*x_scooch,maximum(x)+2*x_scooch,nbins+1)
        y = histcounts(x,binedges) ./ (length(x)*step(binedges))
        panels[i] = Plots.plot(binweave(binedges),interleave(y),yaxis=false,yticks=[],grid=false,label="", xlabel=names[k],
            xlabelfontsize=16,xtickfontsize=12,linewidth=2,linecolor=pltclr, background=bkgrnd,fillcolor=pltclr,fillrange=0,fillalpha=0.1)

# Plot prior distributions
        B = plims[Symbol(k)]
        isequal(k,:cAl) && (B= lNrm(B.μ +log(100),B.σ))
        isequal(k,:R) && (B = lNrm(B.μ +log(1e-3),B.σ))

        if isa(B,Unf)
            prdst = Distributions.Uniform(B.a,ifelse(isinf(B.b),maximum(x),B.b))
        elseif isa(B,Nrm)
            prdst = Distributions.Normal(B.μ,B.σ)
        elseif isa(B,lNrm)
            prdst = Distributions.LogNormal(B.μ,B.σ)
        end
        prdst_x = LinRange(first(binedges),last(binedges),100)
        Plots.plot!(prdst_x,Distributions.pdf.(prdst,prdst_x),linewidth=2,linecolor=pltclr,linestyle=:dash,linealpha=.6)

    end

    sbplts=rows*cols
    Δplts = sbplts-length(panels)
    if Δplts > 0
        blnkplt = Plots.plot(legend=false,grid=false,foreground_color_subplot=bkgrnd,axis=:none,ticks=:none)
        [ push!(panels,blnkplt) for j ∈ 1:Δplts]
    end

    Plots.plot(panels...,layout=Plots.grid(rows,cols),labels="",size=figsize,bottom_margin=10Plots.mm)
end

"""

```julia
proposal_histograms(data_in::Dict,plims::NamedTuple,v::Tuple;
    cols::Integer, nbins::Integer,bounds::Bool,figsize::Tuple, dark=false)
```

Plots a grid of pretty (only top-traced) histograms for parameters listed in `v`.
Parameter `plims` tracks whether log-normal distributions need to be converted back to linear space.

`cols` specifies the number of columns in the histogram grid,
`nbins` specifies the number of bins to plot,
`bounds`=`true` plots the metropolis-enforced bounds in grey (uniform: solid, normal: dash, log-normal: dot-dash),
and `figsize` specifies full figure size in pixels.

Set to a transparent background color scheme by setting `dark=true`

"""
function proposal_histograms(data_in::Dict,plims::NamedTuple,v::Tuple;
    cols::Integer=1, nbins::Integer=20,bounds::Bool=false,cent::Symbol=:none,c_interval=:none,figsize=(800,600),dark::Bool=false)

    names = Dict(   
        :tss => "Age of CAIs (Ma)", 
        :rAlo=>"Initial ²⁶Al/²⁷Al",
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

    nᵥ=length(v)

    d=deepcopy(data_in)
# Convert
    d[:R] .+= log(1e-3) # m -> km
    d[:cAl] .+= log(100) # Multiply by 100 -> wt%, but in natural-log-space
# INCORPORATE ADJUSTMENTS INTO BOUNDS BELOW

#Convert to strings if necessary
    isequal(eltype(keys(d)),String) ? (v= String.(v); acpt = "accept") : ( acpt = :accept)
# Calculate number of rows needed to accomodate all variables in `cols` columns.
    rows = ceil(Int,nᵥ/cols)

    pltclr,bkgrnd = ifelse(dark,(:white,:transparent),(:black,:white))


    panels = Vector{Any}(nothing,nᵥ)
    for i ∈ 1:nᵥ
        k = v[i]
        x = d[k]

        isa(plims[k],lNrm) && map!(exp,x,x)

        x_scooch = (maximum(x)-minimum(x))/ (nbins-4)
        binedges = LinRange(minimum(x)-2*x_scooch,maximum(x)+2*x_scooch,nbins+1)
        y=histcounts(x,binedges)
        panels[i] = Plots.plot(binweave(binedges),interleave(y),yaxis=false,yticks=[],grid=false,label="", xlabel=names[k],
            linewidth=2,linecolor=pltclr, background=bkgrnd,fillcolor=pltclr,fillrange=0,fillalpha=0.1)

        if bounds
            B = plims[Symbol(k)]
            isequal(k,:cAl) && (B= lNrm(B.μ +log(100),B.σ))
            isequal(k,:R) && (B = lNrm(B.μ +log(1e-3),B.σ))

            if isa(B,Unf)
                linestylin= :solid
                bound_lo = fill(B.a,2)
                bound_hi = fill(ifelse(isinf(B.b),maximum(x),B.b),2)
            elseif isa(B,Nrm)
                linestylin= :dash
                bound_lo = fill(B.μ-B.σ,2)
                bound_hi = fill(B.μ+B.σ,2)
            elseif isa(B,lNrm)
                linestylin= :dash
                bound_lo = fill(exp(B.μ-B.σ),2)
                bound_hi = fill(exp(B.μ+B.σ),2)
            end

            Plots.plot!(bound_lo,[0,maximum(y)],linecolor=pltclr,fillalpha=0.6,linewidth=2,linestyle=linestylin)
            Plots.plot!(bound_hi,[0,maximum(y)],linecolor=pltclr,fillalpha=0.6,linewidth=2,linestyle=linestylin)

        end

        if cent==:none

        else
            cent==:mean && (m=mean(x);mname="mean")
            cent==:median  && (m=median(x);mname="median")

            ymax = Plots.ylims(panels[i])[2] # Keep this value fixed.
            Plots.plot!(fill(m,2),[0,ymax],linecolor=pltclr,linewidth=3,label="")
            #Plots.annotate!(m,minimum(y),text(" $mname",:left,:bottom,12))

            if c_interval == :none
            else
                if c_interval == :sigma || c_interval ==:σ
                    cih = cil = std(x)
                elseif isa(c_interval,Number)
                    α = 1-c_interval
                    cil = quantile(x,α/2)
                    cih = quantile(x,1-α/2)
                end
                Plots.plot!(fill(cih,2),[0,ymax],linecolor=pltclr,linewidth=3,linestyle=:dash,label="")
                Plots.plot!(fill(cil,2),[0,ymax],linecolor=pltclr,linewidth=3,linestyle=:dash,label="")
            end
        end

    end

    sbplts=rows*cols
    Δplts = sbplts-length(panels)
    if Δplts > 0
        blnkplt = Plots.plot(legend=false,grid=false,foreground_color_subplot=bkgrnd,axis=:none,ticks=:none)
        [ push!(panels,blnkplt) for j ∈ 1:Δplts]
    end

    Plots.plot(panels...,layout=Plots.grid(rows,cols),labels="",size=figsize)
end


## Plot Evolution of proposals

"""
"""
function plotproposals(d::Dict,plims::NamedTuple,cols::Integer;vars::Tuple=(),ll::Bool=true,bounds::Bool=true,figsize=(800,600))
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
            panels[i] = Plots.plot(x,y,xticks=[],ylabel="ll",linecolor=:black) #use \scrl eventually
            r =  round(100*sum(d[acpt])/length(d[acpt]),digits=1)
            Plots.xlabel!("acceptance = $r %",xguidefontsize=6)
            #Plots.annotate!(last(x), (y[end]+y[1])/2, text("acceptance = $r %", :black,:bottomleft,6))
        elseif isnan(last(y)) || isone(length(y))
            panels[i] = Plots.plot([1,last(x)],fill(y[1],2),xticks=[],ylabel="$k",linecolor=:black)
        else
            panels[i] = Plots.plot(x,y,xticks=[],linecolor=:black)
        end

        if bounds && k != llₛ
            B = plims[Symbol(k)]
            if isa(B,Unf)
                Plots.plot!([1,last(x)],fill(B.a,2),ylabel="$k",linecolor=:grey,linestyle=:solid)
                Plots.plot!([1,last(x)],fill(B.b,2),linecolor=:grey,linestyle=:solid)
            elseif isa(B,Nrm)
                Plots.plot!([1,last(x)],fill(B.μ+B.σ,2),ylabel="$k",linecolor=:grey,linestyle=:dash)
                Plots.plot!([1,last(x)],fill(B.μ-B.σ,2),linecolor=:grey,linestyle=:dash)
            elseif isa(B,lNrm)
                Plots.plot!([1,last(x)],fill(B.μ+B.σ,2),ylabel="log[" * "$k" * "]",linecolor=:grey,linestyle=:dashdot)
                Plots.plot!([1,last(x)],fill(B.μ-B.σ,2),linecolor=:grey,linestyle=:dashdot)
            end
        end
    end
    sbplts=rows*cols
    Δplts = sbplts-length(panels)
    if Δplts > 0
        blnkplt = Plots.plot(legend=false,grid=false,foreground_color_subplot=:white)
        [ push!(panels,blnkplt) for j ∈ 1:Δplts]
    end

    Plots.plot(panels...,layout=Plots.grid(rows,cols),labels="",size=figsize,left_margin=10Plots.mm)
end

"~Plotting Functions Loaded Successfully~"