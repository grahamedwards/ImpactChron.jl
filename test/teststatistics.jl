## Test functions in ImCh_statistics.jl

    # rangemidpoints / rangemidbounds
    # histogramify(!)
    # ll_param
    # ll_params
    # ll_dist

## Test rangemidpoints and rangemidbounds

testrange = LinRange(1., 10., 9)
@test testrange == rangemidbounds(rangemidpoints(testrange))

## Test downscale!

downscale_big = collect(1:90)
no_downscale = similar(downscale_big)
downscale_small = zeros(Int,10)
ImpactChron.downscale!(downscale_small,downscale_big)
ImpactChron.downscale!(no_downscale,downscale_big)
@test downscale_small == [5, 14, 23, 32, 41, 50, 59, 68, 77, 86]
@test downscale_big == no_downscale

## Test histogramify and histogramify!

# Vector-based histogramify(!)
# Set up domain and synthetic vector data that include NaNs extend beyond the domain bounds
hstdmn = 3:3:12.
xhst = [NaN, NaN, 2., 4., 5.5, 6.7, 10.2, 11.8, 14.8, NaN]
yhst = (LinRange(10,100,10).^3 .- LinRange(0,90,10).^3) ./ 100^3
# Calculate the correct answer
hist_test_ans = [yhst[4]+yhst[5],yhst[6],yhst[7]+yhst[8]] / (sum(yhst[.!isnan.(xhst)])*step(hstdmn))
# Let histogramify and histogramify! try...
dist1 = histogramify(hstdmn,xhst,yhst)
dist2=zero(dist1)
histogramify!(dist2,hstdmn,xhst,yhst)
# See how they do!
@test dist1 == dist2 == hist_test_ans

## Test ll_param for Unf

@test iszero(ll_param(5.,Unf(1.,3.)))
@test iszero(ll_param(2.,Unf(1.,3.)))


## Test: ll_param for Nrm & lNrm

@test iszero( ll_param(3.,Nrm(3.,1.)) )
@test iszero( ll_param(3.,lNrm(3.,1.)) )
@test ll_param(2.,Nrm(3.,1.)) === ll_param(2.,Nrm(3.,1.)) === -0.5


## Test ll_params

# Define proposal values
tss=rAlo=R=ta=cAl=Tm=Tc=ρ=Cp=k=tχα=τχα=Fχα=tχβ=τχβ=Fχβ=1.
ϕ = (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ)
# Define proposal distributions
tss=rAlo=R=ta=cAl=Tm=Tc=ρ=Cp=k=tχα=τχα=Fχα=tχβ=τχβ=Fχβ=Nrm(2,1)
ϕdist = (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k,tχα,τχα,Fχα,tχβ,τχβ,Fχβ)

@test ll_params(ϕ,ϕdist) === -0.5 * 10

## Test ll_dist

ˡˡx = 0:.2:10
dNrm = Nrm(5,1)
ˡˡy =  exp.(-( ˡˡx .- dNrm.μ ).^2 ./ (2*dNrm.σ*dNrm.σ)) /(dNrm.σ*sqrt(2π))

μₙₐ=[3.7,4.81,6.5] ; σₙₐ = [.3,.8,.4]
no_aliasing = ll_dist(ˡˡx,ˡˡy,μₙₐ,σₙₐ)

μₐ=[3.53,4.69] ; σₐ = [.03,.04]
aliasing = ll_dist(ˡˡx,ˡˡy,μₐ,σₐ)

@test no_aliasing ≈ -4.877525612777772
@test aliasing ≈ -6.184401753353056
@test (aliasing + no_aliasing) ≈ ll_dist(ˡˡx,ˡˡy,vcat(μₙₐ,μₐ),vcat(σₙₐ,σₐ))

## Test ll_dist_params

ah_lldp = AsteroidHistory(ϕ.R,nnodes=3, Δt=.2,tmax=10.,downscale_factor=1)

# Ensure zeroed agedist returns -Inf
ah_lldp.agedist_downscaled .= zero(eltype(ˡˡy))
@test ll_dist_params(ah_lldp, ϕ,ϕdist, μₙₐ,σₙₐ) === -Inf

# test summing of two ll_ functions.
ah_lldp.agedist_downscaled .= ˡˡy
@test ll_dist_params(ah_lldp, ϕ,ϕdist, μₙₐ,σₙₐ) === ll_dist(ˡˡx,ˡˡy,μₙₐ,σₙₐ) + ll_params(ϕ,ϕdist)

## Test weight_petro_types!

types_wpt = (type3=(T=800., p=1/8), type4=(T=900,p=1/8), type5=(T=1000.,p=1/4), type6=(T=1200,p=1/2))
V_wpt = [0., 0., 0.1,0.12,0.13,0.14,0.15,0.16,0.2]
T_wpt = [1600.,1520, 1150, 1100, 980, 860, 840, 790, 745]

ImpactChron.weight_petro_types!(V_wpt,T_wpt,types_wpt)

# correct answer:
wpt_6 = [0.1,0.12] .* (1/2)/.22; wpt_5 = .13 * (1/4)/.13; wpt_4 = [.14,.15] .* (1/8)/.29; wpt_3 = [.16,.2] .* (1/8)/.36
@test V_wpt ≈ vcat([0.,0.], wpt_6,wpt_5,wpt_4,wpt_3)

V_wpt = [.3,.2,.5]
T_wpt = [ 1150, 1100, 745]

ImpactChron.weight_petro_types!(V_wpt,T_wpt,types_wpt)
@test iszero(V_wpt)