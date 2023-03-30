@test ImpactChron.mcmean([5,3,2],[1,2,1])[1] ≈ 3.333 atol=0.01
@test ImpactChron.mcmean([5,3,2],[1,2,1])[2] ≈ 0.817 atol=0.01

#agerecal could have a test but doesn't for now.

@test ImpactChron.draw((1,2,3)) ∈ (1,2,3)
@test ImpactChron.draw(0.2) === 0.2
@test 1 <= ImpactChron.draw(Unf(1,5)) <= 5
# no test for the normal draw without an explicit rng.

@test ImpactChron.lognorm((x=5,y=7,z=4)).μ ≈ 1.647214140869768
@test ImpactChron.lognorm((x=5,y=7,z=4)).σ ≈ 0.28171393310016424

@test ImpactChron.lognormMC((;n=Nrm(8,0.2), u=Unf(2,10), x=7.3, t = (6,7,8))).μ ≈ 1.9277 atol=0.001
@test ImpactChron.lognormMC((;n=Nrm(8,0.2), u=Unf(2,10), x=7.3, t = (6,7,8))).σ ≈ 0.2651 atol=0.001