using ImpactChron
using Test
#using SafeTestsets --> @safetestset, add "using ImpactChron" after each `begin`.


@testset "perturb" begin include("testImCh_parameters.jl") end
@testset "ImCh_statistics" begin include("testImCh_statistics.jl") end
#@testset
