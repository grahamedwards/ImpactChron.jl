var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ImpactChron","category":"page"},{"location":"#ImpactChron","page":"Home","title":"ImpactChron","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ImpactChron.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ImpactChron]","category":"page"},{"location":"#ImpactChron.Nrm","page":"Home","title":"ImpactChron.Nrm","text":"Nrm(μ::Float64,σ::Float64)\n\nImmutable struct to describe normally distributed data,     reported as mean (μ) and 1σ (σ)\n\n\n\n\n\n","category":"type"},{"location":"#ImpactChron.Unf","page":"Home","title":"ImpactChron.Unf","text":"Unf(a::Float64,b::Float64)\n\nImmutable struct to describe uniformly distributed data,     reported as minimum (a) and maximum (b).\n\n\n\n\n\n","category":"type"},{"location":"#ImpactChron.lNrm","page":"Home","title":"ImpactChron.lNrm","text":"lNrm(μ::Float64,σ::Float64)\n\nImmutable struct to describe lognormally distributed data,     reported as log-space mean (μ) and 1σ (σ)\n\n\n\n\n\n","category":"type"},{"location":"#ImpactChron.ImpactResetQuench-Tuple{AbstractArray, AbstractArray, NamedTuple, NamedTuple}","page":"Home","title":"ImpactChron.ImpactResetQuench","text":"This is the old ImpactRestAr. Needs to be largely rewritten to a framework similar to ImpactResetArray. Most importantly, replace probabilistic maths with number of hits at each timestep maths.\n\nCURRENTLY NOT FUNCTIONAL\n\nFORMER DOCUMENTATION:\n\nImpactResetAr ~ reheat volumes for an exponential impact flux                 described by parameters in p {::Proposal}                 p.tχ ~ instability start time                 p.τχ ~ e-folding timescale of impact flux                 p.Fχ ~ initial impact flux\n\nImpactResetAr(  dates::AbsractArray,                 Vfrxn::AbstractArray,                 p::NamedTuple,                 c::NamedTuple                 Δt::Number,tmax::Number,nᵣ::Integer)\n\nSimulates an impact history from _χ parameters in p, and resets Ar-Ar cooling dates and fractional volumes (Vfraxn) based on impact/crater properties described in c. Impact site morphologies may be described by a conical (cone), parabolic (pbla), or hemispheric (hemi) approximation.\n\nΔt, tmax, and nᵣ respectively define the timestep, model duration, and radial nodes, as in PlntsmlAr function.\n\nNote that Vfrxn is overwritten.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.area_at_depth-Tuple{Number, Number, Cone}","page":"Home","title":"ImpactChron.area_at_depth","text":"area_at_depth(rᵢ::Number, R::Number, x::T) where T<:{Cone,Parabola,Hemisphere}\n\nCalculates the circular area at body radial depth rᵢ (for a body with radius, R) of an impact-affected region, where the volume of the region is approximated by a cone (x::Cone), paraboloid of rotation (x::Parabola), or hemisphere (x::Hemisphere) with a given maximum depth (x.z) and surface radius (x.r).\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.csv2dict-Tuple{String}","page":"Home","title":"ImpactChron.csv2dict","text":"csv2dict(filename::String;symbol::Bool=true)\n\nRead a .csv file into a NamedTuple. Assumes columns reflect individual entries, with the first row specifying that entry's key. For NaN-buffered columns (see nt2csv), the NaNs are removed, returning a single-element entry.\n\nA true value for symbol will convert data saved as arrays of Strings into arrays of Symbols.\n\nsee also: nt2csv, csv2dict\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.csv2nt-Tuple{String}","page":"Home","title":"ImpactChron.csv2nt","text":"csv2nt(filename::String;symbol::Bool=true)\n\nRead a .csv file into a NamedTuple. Assumes columns reflect individual elements, with the first row specifying that element's key. For NaN-buffered columns (see nt2csv), the NaNs are removed, returning a single-element entry.\n\nA true value for symbol will convert data saved as arrays of Strings into arrays of Symbols.\n\nsee also: nt2csv, csv2dict\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.data2csv-Tuple{String, Any}","page":"Home","title":"ImpactChron.data2csv","text":"data2csv(filename::String,data)\n\nSave data to a .csv file. data may be a composite type of Dict or NamedTuple or an Array. Accepts fields of any type element. For single-element entries, the rest of the array/table is filled with NaNs.\n\nsee also: nt2csv, dict2csv\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.dict2csv-Tuple{String, Dict}","page":"Home","title":"ImpactChron.dict2csv","text":"dict2csv(filename::String,D::Dict)\n\nSave a Dict to a .csv. Accounts for fields of any type. For single-element entries, the rest of the array/table is filled with NaNs.\n\nsee also: nt2csv, data2csv\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.histogramify!-Tuple{AbstractVector, AbstractRange, AbstractVector, AbstractVector}","page":"Home","title":"ImpactChron.histogramify!","text":"histogramify!(dist::AbstractVector, domain::AbstractRange, x::AbstractVector, y::AbstractVector)\n\n\nIn-place histogramify that overwites a pre-allocated vector dist.\n\nsee histogramify for details\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.histogramify-Tuple{AbstractRange, AbstractVector, AbstractVector}","page":"Home","title":"ImpactChron.histogramify","text":"histogramify(domain::AbstractRange,x::AbstractVector,y::AbstractVector)\n\nConstructs histogram over (linear) midpoints of domain from model outputs in x with corresponding abundances in y. Does not require a constant step in x.\n\nNormalizes the output, such that for output dist ∑ dist[dᵢ] * Δd = 1 (for each dᵢ in the bincenters of domain with step-size Δd), so long as all x ∈ domain. If any x ∉ domain, ∑ dist[dᵢ] * Δd = 1- (∑yₒᵤₜ / ∑yₐₗₗ ) where the corresponding xₒᵤₜ of each yₒᵤₜ is ∉ domain.\n\n\n\n\n\nReturns only the histogram masses, centers of time bins must be calculated externally. e.g.\n\nΔd = step(domain)\nbincenters= LinRange( first(domain)+Δd/2, last(domain)-Δd/2, length(domain)-1)\n\nor see rangemidpoints\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.impact_reset_array!-Tuple{AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, NamedTuple, NamedTuple}","page":"Home","title":"ImpactChron.impact_reset_array!","text":"impact_reset_array!(tₓr::AbstractArray, solartime::AbstractRange, impacts::AbstractArray, tcoolₒ::AbstractArray,\n                            dates::AbstractArray,Vfrxn::AbstractArray,\n                            p::NamedTuple,c::NamedTuple;\n                                nᵣ::Integer)\n\nSimulates an impact history from -χ parameters in p (below), and resets Ar-Ar primary planetesimal cooling dates and fractional volumes (Vfraxn) based on impact/crater properties described in c.\n\nDepth-cooling age (relative) abundances are tracked in array tₓr (time x radial depth), with dimensions (time,radius) = (length(solartime),nᵣ), where nᵣ describes the number of radial nodes, as in the PlntsmlAr function.\n\nimpacts and tcoolₒ are pre-allocated vectors that respectively record the number of impacts at each time step and the index of the primary cooling date in solartime.\n\nImpact site geometries are codified by types (Cone,Parabola,Hemisphere) defined in c and calculated in the area_at_depth.\n\nImpact flux follows an exponential decay described by parameters in p:\n\np.tχ ~ instability start time\n\np.τχ ~ e-folding timescale of impact flux\n\np.Fχ ~ initial impact flux\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.ll_dist-Tuple{AbstractRange, AbstractVector, AbstractVector, AbstractVector}","page":"Home","title":"ImpactChron.ll_dist","text":"ll_dist(x::AbstractVector,dist::AbstractVector,mu::AbstractVector,sigma::AbstractVector)\n\nCalculate loglikelihood that observations in mu and sigma are drawn from modeled distribution described by x and dist where x contains the bincenters of a normalized histogram dist and mu and sigma respectively contain the mean and 1σ of the observations.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.ll_param-Tuple{Number, Nrm}","page":"Home","title":"ImpactChron.ll_param","text":"```julia ll_param(x::Number,D::T) -> T ∈ {Nrm,lNrm,Unf}\n\nCalculate the log-likelihood that x is drawn from a distribution D, where D may be... Normal (Nrm), with mean D.μ and 1σ = D.σ\n\nLognormal (lNrm), with logspace mean D.μ and 1σ = D.σ. Assumes x is already in logspace.\n\nUniform (Unf) with lowerbound D.a and upperbound D.b (D::Unf always returns loglikelihhood of zero, bounds test is done earlier to speed up code.)\n\nsee ImCh_parameters.jl for construction of T-structs\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.ll_params-Tuple{NamedTuple, NamedTuple}","page":"Home","title":"ImpactChron.ll_params","text":"ll_params(p::NamedTuple,d::NamedTuple)\n\nCalculate log-likelihood for a number of proposals in p with corresponding distributions in d\n\nCurrently does not evaluate impact (χ) parameters.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.perturb-Tuple{NamedTuple, Symbol, Number}","page":"Home","title":"ImpactChron.perturb","text":"perturb(p::NamedTuple,k::Symbol,n::Number)\n\nReturn a NamedTuple identical to p, with one field (key k) changed to the value of n. Note that == identity is preserved only if the order of fields in p is as below\n\nFields: tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.rangemidbounds-Tuple{AbstractRange}","page":"Home","title":"ImpactChron.rangemidbounds","text":"rangebinbounds(x::AbstractRange)\n\nCalculate a LinRange of the linear bounds for each \"midpoint\" step in x.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.rangemidpoints-Tuple{AbstractRange}","page":"Home","title":"ImpactChron.rangemidpoints","text":"rangemidpoints(x::AbstractRange)\n\nCalculate a LinRange of the midpoints of each step in x(<:AbstractRange).\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.tturbosum-Tuple{AbstractArray}","page":"Home","title":"ImpactChron.tturbosum","text":"tturbosum(x::AbstractArray)\n\nFast summing of x with the power of LoopVectorization.jl's @tturbo (multithreaded turbo).\n\n\n\n\n\n","category":"method"}]
}
