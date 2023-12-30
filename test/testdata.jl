
maxage= 4000
trimtest = ImpactChron.trim_ages(maxage,ImpactChron.readdlm(string(@__DIR__,"/../data/KArArages.csv"),','))
# test trim_ages output    
@test typeof(trimtest) <:  Tuple{Vector{Float64}, Vector{Float64}, Vector{String}, Vector{String}, Vector{String}}
# test that trim_ages removes all ages younger than maxage
@test isempty(findall(x -> x<maxage, trimtest[1]))

x = ImpactChron.loadArAr(0)
    @test 136 == length(x[1]) == length(x[2]) == length(x[3]) == length(x[4]) == length(x[5])

x = ImpactChron.loadKAr(0)
    @test 67 == length(x[1]) == length(x[2]) == length(x[3]) == length(x[4]) == length(x[5])

x = ImpactChron.loadKArAr(0)
    @test 67+136 == length(x[1]) == length(x[2]) == length(x[3]) == length(x[4]) == length(x[5])

x = ImpactChron.loadECs(0)
    @test 27 == length(x[1]) == length(x[2]) == length(x[3]) == length(x[4]) == length(x[5])

x = ImpactChron.loadTrie2003(0)
    @test 9 == length(x[1]) == length(x[2]) == length(x[3]) == length(x[4]) == length(x[5])