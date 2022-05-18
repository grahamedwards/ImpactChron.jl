## Functions to aid visualization of Metropolis code outputs.
using Plots;gr()

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
