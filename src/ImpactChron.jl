module ImpactChron

# For speed:
using LoopVectorization: @turbo, @tturbo, vreduce
using Polyester: @batch
using VectorizedStatistics: vsum, vmean, vstd, vmedian, vquantile

# For data handling and plotting
using DelimitedFiles, Serialization
import Requires

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


export load_datamgmt
"""

    load_datamgmt()

Load data management tools for summary statistics and plotting of ImpactChron data. 

"""
load_datamgmt() = include(string(@__DIR__,"../visualize/datamgmt.jl"))

function __init__()
    vizpath = string(@__DIR__,"../visualize/plotting.jl")
    Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include(vizpath)
    Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include(vizpath) 
    Requires.@require WGLMakie="276b4fcb-3e11-5398-bf8b-a0c2d153d008" include(vizpath) 
end


end
