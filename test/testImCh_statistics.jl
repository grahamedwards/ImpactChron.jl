## Test functions in ImCh_statistics.jl
    # histogramify(!)
    # ll_param
    # ll_params
    # ll_dist

## Test histogramify and histogramify!
# Set up domain and synthetic data that include NaNs extend beyond the domain bounds
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
