using ImpactChron
using Test
using SafeTestsets # @safetestset, add "using ImpactChron" after each `begin`.


@testset "perturb" begin include("testparameters.jl") end
@testset "ImCh_statistics" begin include("teststatistics.jl") end
@safetestset "thermal" begin using ImpactChron; include("testthermal.jl") end
