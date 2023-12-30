module ImpactChron

# For speed:
using LoopVectorization: @turbo, @tturbo, vreduce
using Polyester: @batch
using VectorizedStatistics: vsum, vmean, vstd, vmedian, vquantile

# For data handling.
using DelimitedFiles
using Serialization

# For reproducibility
using Random

export  perturb, PriorDistribution, Nrm, lNrm, Unf, ImpactSiteShape, Cone, Parabola, Hemisphere, ImpactSite, AsteroidHistory, PetroTypes
include("parameters.jl")

export  rangemidpoints, rangemidbounds, ll_param, ll_params, ll_dist, ll_dist_params
include("statistics.jl")

export plntsml_Tz, planetesimal_cooling_dates, planetesimal_cooling_dates!, impact_reset_array!, asteroid_agedist!
include("thermal.jl")

export thermochron_metropolis
include("metropolis.jl")

export nt2csv, dict2csv, data2csv, csv2nt, csv2dict, loadArAr, loadKArAr, loadKAr
include("datasaveload.jl")

include("calcpriors.jl") # all functions are loaded explicitly in the scripts in the data/ directory.

end
