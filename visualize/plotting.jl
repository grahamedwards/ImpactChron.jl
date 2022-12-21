## Functions to produce visualizations of of ImpactChron outputs

using Plots;gr()

isdefined(@__MODULE__,:interleave) || throw("Load `datamgmt.jl` first.")


"""

```julia
proposal_histograms(data_in::Dict,plims::NamedTuple,v::Tuple;
    cols::Integer, nbins::Integer,bounds::Bool,figsize::Tuple)
```

Plots a grid of pretty (only top-traced) histograms for parameters listed in `v`.
Parameter `plims` tracks whether log-normal distributions need to be converted back to linear space.

`cols` specifies the number of columns in the histogram grid,
`nbins` specifies the number of bins to plot,
`bounds`=`true` plots the metropolis-enforced bounds in grey (uniform: solid, normal: dash, log-normal: dot-dash),
and `figsize` specifies full figure size in pixels.

"""
function proposal_histograms(data_in::Dict,plims::NamedTuple,v::Tuple;
    cols::Integer=1, nbins::Integer=20,bounds::Bool=false,cent::Symbol=:none,c_interval=:none,figsize=(800,600))

    nᵥ=length(v)

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
        :tχα=> "Post-accretion bombardment onset (Ma after CAIs)",
        :τχα=> "Post-accretion bombardment ℯ-folding time (Ma)",
        :Fχα=> "Post-accretion initial impactor flux (Ma⁻¹)",
        :tχβ=> "Primordial bombardment onset (Ma after CAIs)",
        :τχβ=> "Primordial bombardment ℯ-folding time (Ma)",
        :Fχβ=> "Primordial initial impactor flux (Ma⁻¹)"             )

    d=deepcopy(data_in)
# Convert
    d[:R] .+= log(1e-3) # m -> km
    d[:cAl] .+= log(100) # Multiply by 100 -> wt%, but in natural-log-space
# INCORPORATE ADJUSTMENTS INTO BOUNDS BELOW

#Convert to strings if necessary
    isequal(eltype(keys(d)),String) ? (v= String.(v); acpt = "accept") : ( acpt = :accept)
# Calculate number of rows needed to accomodate all variables in `cols` columns.
    rows = ceil(Int,nᵥ/cols)

    panels = Vector{Any}(nothing,nᵥ)
    for i ∈ 1:nᵥ
        k = v[i]
        x = d[k]

        isa(plims[k],lNrm) && map!(exp,x,x)

        x_scooch = (maximum(x)-minimum(x))/ (nbins-4)
        binedges = LinRange(minimum(x)-2*x_scooch,maximum(x)+2*x_scooch,nbins+1)
        y=histcounts(x,binedges)
        panels[i] = plot(binweave(binedges),interleave(y),yaxis=false,yticks=[],grid=false,label="", xlabel=names[k],
            linewidth=2,linecolor=:black,fillcolor=:black,fillrange=0,fillalpha=0.1)

        if bounds
            B = plims[Symbol(k)]
            isequal(k,:cAl) && (B= lNrm(B.μ +log(100),B.σ))
            isequal(k,:R) && (B = lNrm(B.μ +log(1e-3),B.σ))

            if isa(B,Unf)
                linestylin= :solid
                bound_lo = fill(B.a,2)
                bound_hi = fill(B.b,2)
            elseif isa(B,Nrm)
                linestylin= :dash
                bound_lo = fill(B.μ-B.σ,2)
                bound_hi = fill(B.μ+B.σ,2)
            elseif isa(B,lNrm)
                linestylin= :dash
                bound_lo = fill(exp(B.μ-B.σ),2)
                bound_hi = fill(exp(B.μ+B.σ),2)
            end

            plot!(bound_lo,[0,maximum(y)],linecolor=:grey60,linewidth=2,linestyle=linestylin)
            plot!(bound_hi,[0,maximum(y)],linecolor=:grey60,linewidth=2,linestyle=linestylin)

        end

        if cent==:none

        else
            cent==:mean && (m=mean(x);mname="mean")
            cent==:median  && (m=median(x);mname="median")

            ymax = ylims(panels[i])[2] # Keep this value fixed.
            plot!(fill(m,2),[0,ymax],linecolor=:black,linewidth=3,label="")
            #annotate!(m,minimum(y),text(" $mname",:left,:bottom,12))

            if c_interval == :none
            else
                if c_interval == :sigma || c_interval ==:σ
                    cih = cil = std(x)
                elseif isa(c_interval,Number)
                    α = 1-c_interval
                    cil = quantile(x,α/2)
                    cih = quantile(x,1-α/2)
                end
                plot!(fill(cih,2),[0,ymax],linecolor=:black,linewidth=3,linestyle=:dash,label="")
                plot!(fill(cil,2),[0,ymax],linecolor=:black,linewidth=3,linestyle=:dash,label="")
            end
        end

    end

    sbplts=rows*cols
    Δplts = sbplts-length(panels)
    if Δplts > 0
        blnkplt = plot(legend=false,grid=false,foreground_color_subplot=:white,axis=:none,ticks=:none)
        [ push!(panels,blnkplt) for j ∈ 1:Δplts]
    end

    plot(panels...,layout=grid(rows,cols),labels="",size=figsize)
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
            panels[i] = plot(x,y,xticks=[],ylabel="ll",linecolor=:black) #use \scrl eventually
            r =  round(100*sum(d[acpt])/length(d[acpt]),digits=1)
            xlabel!("acceptance = $r %",xguidefontsize=6)
            #annotate!(last(x), (y[end]+y[1])/2, text("acceptance = $r %", :black,:bottomleft,6))
        elseif isnan(last(y)) || isone(length(y))
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

    plot(panels...,layout=grid(rows,cols),labels="",size=figsize)
end

"~Plotting Functions Loaded Successfully~"