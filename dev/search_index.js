var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ImpactChron","category":"page"},{"location":"#ImpactChron","page":"Home","title":"ImpactChron","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ImpactChron.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ImpactChron]","category":"page"},{"location":"#ImpactChron.AsteroidHistory","page":"Home","title":"ImpactChron.AsteroidHistory","text":"AsteroidHistory(typeseed<:Number; nnodes, Δt, tmax, downscale_factor)\n\nstruct containing Arrays that record the evolution of a (bombarded) asteroid.\n\nConstructor function takes a parameter from the proposal to seed type and requires the number of radial nodes nnodes (::Int), the timestep used Δt, the full time duration tmax, and the downscale_factor (::Int).\n\nFields:      Vfrxn: volume fractions of each radial shell,     peakT: peak temperature of each radial shell,      cooltime: indices in t of the primary cooling date,      impacts: number of impacts at each timestep,      txr: time x radius array of proportional cooling ages,      agedist: distribution of ages corresponding to t,      agedist_downscaled: distribution of ages corresponding to t_downscaled,      t: timesteps of full-scale model,      t_downscaled: timesteps of downscaled model output, \n\n\n\n\n\n","category":"type"},{"location":"#ImpactChron.ImpactSite-Union{Tuple{T}, Tuple{Type{T}, Number}} where T","page":"Home","title":"ImpactChron.ImpactSite","text":"ImpactSite(heat<:ImpactSiteShape,C<:Number)\n\nstruct describing the shape of the simulated impact site and zone of reheating.\n\nCONSTRUCTOR FUNCTION\n\nImpactSite(shape,impactor_diameter)\n\nProviding an impactor_diameter (<:Number) calculates impact parameters based on approximate heat distribution modeled in Davison+ 2012 (GCA, http://dx.doi.org/10.1016/j.gca.2012.08.001).\n\nImpactSite(shape ; r, C)\n\nIf values of r and C are provided,  prepares an ImpactSite that extends to the center of an asteroid with radius r (km) and has a surface diameter C times the asteroid circumference (C=0.01 by default). Note that C ∈ [0,1]. If no r is provided this seeds an ImpactSite with zeroed ImpactSiteShape and a C value. \n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.Nrm","page":"Home","title":"ImpactChron.Nrm","text":"Nrm(μ::Float64,σ::Float64)\n\nImmutable struct to describe normally distributed data,     reported as mean (μ) and 1σ (σ)\n\n\n\n\n\n","category":"type"},{"location":"#ImpactChron.PetroTypes","page":"Home","title":"ImpactChron.PetroTypes","text":"PetroTypes(temps::NamedTuple,samples::Vector{String})\n\nstruct containing fields of type3-type6, reflecting petrologic types, each with subfilds T (maximum temperature in K) and p (proportion among the chondrite record). Note that ps sum to unity.\n\nConstructor takes a NamedTuple containing maximum temperatures as fields T3-T6 and a Vector{String} containing the petrologic types corresponding to chondrites used as a prior.\n\nIf no argument given –  PetroTypes() – returns zeroed type_ fields and weight=false, which short-circuits weighting by petrologic type. \n\n\n\n\n\n","category":"type"},{"location":"#ImpactChron.Unf","page":"Home","title":"ImpactChron.Unf","text":"Unf(a::Float64,b::Float64)\n\nImmutable struct to describe uniformly distributed data,     reported as minimum (a) and maximum (b).\n\n\n\n\n\n","category":"type"},{"location":"#ImpactChron.lNrm","page":"Home","title":"ImpactChron.lNrm","text":"lNrm(μ::Float64,σ::Float64)\n\nImmutable struct to describe lognormally distributed data,     reported as log-space mean (μ) and 1σ (σ)\n\n\n\n\n\n","category":"type"},{"location":"#ImpactChron.agerecal-Tuple{Number, Number}","page":"Home","title":"ImpactChron.agerecal","text":"agerecal(x,sig;monitor_age,n, KAr=false)\n\nRecalculate the age and uncertainty of an Ar-Ar age x±sig (1σ) Ma with the decay constants of Steiger & Jager (1977). Quantifies uncertainty by resampling n times (default: 1M). Optionally accepts monitor_age. If not given, resamples from a range of absolute (not necessarily K-Ar constrained) monitor ages ∈ [0,3] Ga.\n\nIf KAr=true, this recalculates the K-Ar age with the new decay constants, without correcting for a monitor.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.csv2dict-Tuple{String}","page":"Home","title":"ImpactChron.csv2dict","text":"csv2dict(filename::String;symbol::Bool=true)\n\nRead a .csv file into a NamedTuple. Assumes columns reflect individual entries, with the first row specifying that entry's key. For NaN-buffered columns (see nt2csv), the NaNs are removed, returning a single-element entry.\n\nA true value for symbol will convert data saved as arrays of Strings into arrays of Symbols.\n\nsee also: nt2csv, csv2dict\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.csv2nt-Tuple{String}","page":"Home","title":"ImpactChron.csv2nt","text":"csv2nt(filename::String;symbol::Bool=true)\n\nRead a .csv file into a NamedTuple. Assumes columns reflect individual elements, with the first row specifying that element's key. For NaN-buffered columns (see nt2csv), the NaNs are removed, returning a single-element entry.\n\nA true value for symbol will convert data saved as arrays of Strings into arrays of Symbols.\n\nsee also: nt2csv, csv2dict\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.data2csv-Tuple{String, Any}","page":"Home","title":"ImpactChron.data2csv","text":"data2csv(filename::String,data)\n\nSave data to a .csv file. data may be a composite type of Dict or NamedTuple or an Array. Accepts fields of any type element. For single-element entries, the rest of the array/table is filled with NaNs.\n\nsee also: nt2csv, dict2csv\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.dict2csv-Tuple{String, Dict}","page":"Home","title":"ImpactChron.dict2csv","text":"dict2csv(filename::String,D::Dict)\n\nSave a Dict to a .csv. Accounts for fields of any type. For single-element entries, the rest of the array/table is filled with NaNs.\n\nsee also: nt2csv, data2csv\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.downscale!-Union{Tuple{T}, Tuple{AbstractArray{T}, Any}} where T","page":"Home","title":"ImpactChron.downscale!","text":"downscale!(B::AbstractArray, A::AbstractArray)\n\nDownscales elements of 1-D array A into smaller 1-D array B by summing. Scales the downscaled values in B by the ratio of length(A)÷length(B) to preserve any normalizations.\n\nRequires that length(A) % length(B) = 0.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.histogramify!-Tuple{AbstractVector, AbstractRange, AbstractVector, AbstractVector}","page":"Home","title":"ImpactChron.histogramify!","text":"histogramify!(dist::AbstractVector, domain::AbstractRange, x::AbstractVector, y::AbstractVector)\n\n\nIn-place histogramify that overwites a pre-allocated vector dist.\n\nsee histogramify for details\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.histogramify-Tuple{AbstractRange, AbstractVector, AbstractVector}","page":"Home","title":"ImpactChron.histogramify","text":"histogramify(domain::AbstractRange,x::AbstractVector,y::AbstractVector)\n\nConstructs histogram over (linear) midpoints of domain from model outputs in x with corresponding abundances in y. Does not require a constant step in x.\n\nNormalizes the output, such that for output dist ∑ dist[dᵢ] * Δd = 1 (for each dᵢ in the bincenters of domain with step-size Δd), so long as all x ∈ domain. If any x ∉ domain, ∑ dist[dᵢ] * Δd = 1- (∑yₒᵤₜ / ∑yₐₗₗ ) where the corresponding xₒᵤₜ of each yₒᵤₜ is ∉ domain.\n\n\n\n\n\nReturns only the histogram masses, centers of time bins must be calculated externally. e.g.\n\nΔd = step(domain)\nbincenters= LinRange( first(domain)+Δd/2, last(domain)-Δd/2, length(domain)-1)\n\nor see rangemidpoints\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.impact_reset_array!-Tuple{AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, NamedTuple, ImpactSite}","page":"Home","title":"ImpactChron.impact_reset_array!","text":"function impact_reset_array!(tₓr::AbstractArray,solartime::AbstractArray,tcoolₒ::AbstractArray,Vfrxn::AbstractArray,\n                                impacts::AbstractArray,\n                                p::NamedTuple,c::NamedTuple;\n                                nᵣ::Integer,Δt::Number)\n\nSimulates an impact history from -χ parameters in p (below), and resets primary planetesimal cooling dates (indices of dates in solartime in tcoolₒ) and fractional volumes (Vfraxn) based on impact/crater properties described in c.\n\nDepth-cooling age (relative) abundances are tracked in array tₓr (time x radial depth), with dimensions (time,radius) = (length(solartime),nᵣ), where nᵣ describes the number of radial nodes, as in the planetesimal_cooling_dates function.\n\nimpacts and tcoolₒ are pre-allocated vectors that respectively record the number of impacts at each time step and the index of the primary cooling date in solartime.\n\nImpact site geometries are codified by types (Cone,Parabola,Hemisphere) defined in c and calculated in the radius_at_depth.\n\nImpact flux follows an exponential decay described by parameters in p:\n\np.tχ ~ instability start time\n\np.τχ ~ e-folding timescale of impact flux\n\np.Fχ ~ initial impact flux\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.js_arx2dict-Tuple{String, Tuple}","page":"Home","title":"ImpactChron.js_arx2dict","text":"\njs_arx2dict(file::String, vars::Tuple; n, ll=true, accept=true, perturbation=false)\n\n\nLoad a serialized archive file from thermochron_metropolis into a Dictfor \"post-run\" analysis, only incorporating the variables invarsand coveringnsteps (all steps by default). Provide aBoolto incorporatell,accept, orperturbation`\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.ll_dist-Tuple{AbstractRange, AbstractVector, AbstractVector, AbstractVector}","page":"Home","title":"ImpactChron.ll_dist","text":"ll_dist(x::AbstractVector,dist::AbstractVector,mu::AbstractVector,sigma::AbstractVector)\n\nCalculate loglikelihood that observations in mu and sigma are drawn from modeled distribution described by x and dist where x contains the bincenters of a normalized histogram dist and mu and sigma respectively contain the mean and 1σ of the observations.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.ll_dist_params-Tuple{AsteroidHistory, NamedTuple, NamedTuple, AbstractArray, AbstractArray}","page":"Home","title":"ImpactChron.ll_dist_params","text":"function ll_dist_params(a::AsteroidHistory, p, plims, mu,sig)\n\nCalculate the combined log-likelihood that \n\nObservations with mean mu and corresponding 1σ uncertainties sig (both ::Array) are drawn from the downscaled age distribution and timesteps contained in a (::AsteroidHistory)\nThe proposal parameters p used to calculate this history are drawn from the prior distributions plims (both ::NamedTuple)\n\nIf the age distribution is all zero, quickly returns a log-liklihood of -∞\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.ll_param-Tuple{Number, Nrm}","page":"Home","title":"ImpactChron.ll_param","text":"```julia ll_param(x::Number,D::T) -> T ∈ {Nrm,lNrm,Unf}\n\nCalculate the log-likelihood that x is drawn from a distribution D, where D may be... Normal (Nrm), with mean D.μ and 1σ = D.σ\n\nLognormal (lNrm), with logspace mean D.μ and 1σ = D.σ. Assumes x is already in logspace.\n\nUniform (Unf) with lowerbound D.a and upperbound D.b (D::Unf always returns loglikelihhood of zero, bounds test is done earlier to speed up code.)\n\nsee ImCh_parameters.jl for construction of T-structs\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.ll_params-Tuple{NamedTuple, NamedTuple}","page":"Home","title":"ImpactChron.ll_params","text":"ll_params(p::NamedTuple,d::NamedTuple)\n\nCalculate log-likelihood for a number of proposals in p with corresponding distributions in d\n\nCurrently does not evaluate impact (χ) parameters.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.mcmean-Tuple{AbstractArray, AbstractArray}","page":"Home","title":"ImpactChron.mcmean","text":"mcmean(x,xsig,n=10_000_000,fullpost=false)\n\nCalculate mean of age(s) x with 1σ uncertainties xsig by Monte Carlo method. A 10M resample (default) returns consistent results at the 10ka level.\n\nReturns a tuple of (mean,1σ) by default. \n\nIf fullpost=true, returns posterior samples rather than summary statistics.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.nan_regolith!-Tuple{AbstractArray, AbstractArray, Number}","page":"Home","title":"ImpactChron.nan_regolith!","text":"nan_regolith!(d::AbstractArray,T::AbstractArray,Tmin::Number)\n\nReplace all dates in d with NaN where corresponding peak temperatures in T are less than Tmin Fast with @turbo\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.perturb-Tuple{NamedTuple, Symbol, Number}","page":"Home","title":"ImpactChron.perturb","text":"perturb(p::NamedTuple,k::Symbol,n::Number)\n\nReturn a NamedTuple identical to p, with one field (key k) changed to the value of n. Note that == identity is preserved only if the order of fields in p is as below\n\nFields: tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ,tχγ,τχγ,Fχγ\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.planetesimal_cooling_dates!-Tuple{AbstractArray, AbstractArray, AbstractArray, NamedTuple}","page":"Home","title":"ImpactChron.planetesimal_cooling_dates!","text":"function planetesimal_cooling_dates!(ages::AbstractArray, Vfrxn::AbstractArray,peakT::AbstractArray, p::NamedTuple;\n    nᵣ::Integer, Δt::Number, tmax::Number, Tmax::Number, Tmin::Number)\n\nIn-place planetesimal_cooling_dates that updates ages, Vfrxn, and peakT.\n\nsee also planetesimal_cooling_dates\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.planetesimal_cooling_dates-Tuple{}","page":"Home","title":"ImpactChron.planetesimal_cooling_dates","text":"planetesimal_cooling_dates(p::NamedTuple; nᵣ, Δt, tmax, Tmax, Tmin)\n\nReturns an NTuple{4,Vector} containing (in this order) the thermochronologic cooling dates (Ma after CAIs) and corresponding volume fractions, radial depths (km from center), and peak temperatures (K) of nᵣ evenly spaced nodes in a spherical body. \n\nPhysical and environmental parameters are described in p. Alternatively, these parameters may be individually listed in lieu of ϕ. These parameters are outlined in the table below. Note that several of these parameters need to be entered as the natural logarithm of the value for easy compatibility with the inversion function.\n\nΔt gives the timestep (in Ma), tmaxdescribes the duration of the model (Ma after CAIs), andTmaxandTmindefine the maximum and minimum temperatures (K) allowed for chondritic material in the body. Default values are only given fortmax(2000 Ma),Tmax(1500 K), andTmin` (0 K).\n\n| Parameter                 | log?  | `NmTpl` | `func`  |\n| ------------------------- | ----  | ------ | -------- |\n| closure temperature (K)   | no    | `Tc`   | `Tc`     |\n| solar system age (Ma)     | no    | `tss`  | `tₛₛ`    |\n| initial ²⁶Al/²⁷Al         | no    | `rAlo` | `rAlo`   |\n| body radius (m)           | yes   | `R`    | `R`      |\n| accretion date (Ma)       | yes   | `ta`   | `tₐ`     |\n| disk temperature (K)      | yes   | `Tm`   | `To`     |\n| [Al] (g/g)                | yes   | `cAl`  | `Al_conc`|\n| density (kg/m³)           | yes   | `k`    | `K`      |\n| thermal diffusivity       | yes   | `ρ`    | `ρ`      |\n| spec. heat capacity       | yes   | `Cp`   | `Cₚ`     |\n| ------------------------- | ----  | ------ | -------- |\n\nsee also: planetesimal_cooling_dates!\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.planetesimal_cooling_timestep!-Tuple{AbstractRange, AbstractVector, AbstractVector, AbstractVector, NamedTuple}","page":"Home","title":"ImpactChron.planetesimal_cooling_timestep!","text":"function planetesimal_cooling_timestep!(solartime::AbstractRange, time_i::Vector, Vfrxn::Vector, peakT::Vector, p::NamedTuple; nᵣ, Tmax, Tmin)\n\nReturns thermochronologic cooling dates in time_i as indices of solartime, along with corresponding volumetric fractions (Vfrxn) and peak temperatures in K (peakT) for nᵣ nodes in a body with planetesimal and environmental parameters given in p. Tmax and Tmin respectively describe the maximum and minimum temperatures allowed in the chondritic planetesimal. Failing to exceed Tmin gives the date of accretion, and exceeding Tmax sets the volumetric fraction to zero (achondritic).\n\nsee also: planetesimal_cooling_dates, planetesimal_cooling_dates!\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.planetesimal_temperature-Tuple{AbstractArray, AbstractArray}","page":"Home","title":"ImpactChron.planetesimal_temperature","text":"function planetesimal_temperature(time::AbstractArray, radii::AbstractArray; To::Float64, Ao::Float64, λ::Float64, K::Float64, κ::Float64 )\n\nCalculates the evolution of temperature at a range of depths defined for a conductively cooling sphere with thermal conductivity K and thermal diffusivity κ, given ambient temperature To, initial heat production Ao, and heat-producing-element decay constant 'λ'.\n\nAdapted from: Carlslaw & Jäger (1959) and Hevey & Sanders (2006)\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.prior_bounds-Tuple{Number, Unf}","page":"Home","title":"ImpactChron.prior_bounds","text":"ImpactChron.prior_bounds(x,p<:PriorDistribution)\n\nEvaluates whether a proposal x falls within the permissible bounds of its prior p <: PriorDistribution.  Always returns true if p is of type Nrm or lNrm, tests if x ∈ (p.a,p.b) for p::Unf.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.radius_at_depth-Tuple{Number, Number, Cone}","page":"Home","title":"ImpactChron.radius_at_depth","text":"radius_at_depth(rᵢ::Number, R::Number, x::T) where T<:{Cone,Parabola,Hemisphere}\n\nCalculates the radius of a circular area at a radial distance of rᵢ from the center of a body with radius R, where the volume of the region is approximated by a cone (x::Cone), paraboloid of rotation (x::Parabola), or hemisphere (x::Hemisphere). x includes a maximum depth (x.z) and surface radius (x.r). When x::Hemisphere, only r is used.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.rangemidbounds-Tuple{AbstractRange}","page":"Home","title":"ImpactChron.rangemidbounds","text":"rangebinbounds(x::AbstractRange)\n\nCalculate a LinRange of the linear bounds for each \"midpoint\" step in x.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.rangemidpoints-Tuple{AbstractRange}","page":"Home","title":"ImpactChron.rangemidpoints","text":"rangemidpoints(x::AbstractRange)\n\nCalculate a LinRange of the midpoints of each step in x(<:AbstractRange).\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.strict_priors-Tuple{NamedTuple, Symbol, ImpactChron.PriorDistribution}","page":"Home","title":"ImpactChron.strict_priors","text":"strict_priors(p::NamedTuple,k,p_prior<:PriorDistribution)\n\nEvaluates strict priors related to Uniform distributions and other rules for paramter proposal p and perturbed variable k::Symbol, where p_prior is the prior distribution of p[k].\n\nReturns true if all priors are satisfied. Returns false if any priors fail.\n\nCurrently includes:\n\n1. Ensure bounds of uniform priors are not exceeded.\n\n2. Ensure bombardment events α, β, γ are in sequential order.\n\n3. The ℯ-folding time of the primordial flux must be longer than that of post-accretion bombardments.\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.tturbosum-Tuple{AbstractArray}","page":"Home","title":"ImpactChron.tturbosum","text":"tturbosum(x::AbstractArray)\n\nFast summing of x with the power of LoopVectorization.jl's @tturbo (multithreaded turbo).\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.turbosum-Tuple{AbstractArray}","page":"Home","title":"ImpactChron.turbosum","text":"turbosum(x::AbstractArray)\n\nFast summing of x with the power of LoopVectorization.jl's @turbo. Since vreduce does not accept views, we turbosum!\n\nFaster than reduce for length(x)>≈200\n\n\n\n\n\n","category":"method"},{"location":"#ImpactChron.weight_petro_types!-Tuple{AbstractArray, AbstractArray, PetroTypes}","page":"Home","title":"ImpactChron.weight_petro_types!","text":"weight_petro_types!(v::AbstractArray,T::AbstractArray,petrotypes::NamedTuple)\n\nReweight volumetric fractions relative to the abundance of each petrologic type in the meteorite record. Takes Vectors of volumetric fraction (v), peak temperature (T), and cooling date (d), as output by planetesimal_cooling_dates. Accounts for melted layers (r where v[i]==0 due to melting in planetesimal_cooling_timestep!. Requires all petrologic types to occupy at least one radial node, otherwise returns zero in all v.\n\npetrotypes is a NamedTuple of NamedTuples, such that petrotypes = (type#=(T<:Number,p<:Number), ...) where T and p respectively represent maximum temperature and relative abundance in the meteorite record.\n\n\n\n\n\n","category":"method"}]
}
