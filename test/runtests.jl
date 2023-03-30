using ImpactChron
using Test
using SafeTestsets # @safetestset, add "using ImpactChron" after each `begin`.

@testset "perturb & Constructors" begin include("testparameters.jl") end
@testset "statistics" begin include("teststatistics.jl") end
@testset "prior calculations" begin include("testcalcpriors.jl") end
@safetestset "thermal" begin include("testthermal.jl") end
@safetestset "metropolis" begin include("testmetropolis.jl") end
