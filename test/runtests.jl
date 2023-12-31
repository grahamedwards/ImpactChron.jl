using ImpactChron
using Test
using StableRNGs
using Suppressor

@testset "perturb & Constructors" begin include("testparameters.jl") end
@testset "data management" begin include("testdata.jl") end
@testset "statistics" begin include("teststatistics.jl") end
@testset "prior calculations" begin include("testcalcpriors.jl") end
@testset "thermal" begin include("testthermal.jl") end
@testset "metropolis" begin include("testmetropolis.jl") end

