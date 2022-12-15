## Calculate or recalculate mean ages for specific chondrites and chondrite groups (1σ uncertainties)

# The mean ages and uncertainties from this file are used in "KArArages_compilation.csv"

using ImpactChron: mcmean

Shaw = mcmean([4420.,4400],[30,30.]) 

Y75097 = mcmean([489.4,505.2],[3.,2.7]) 

Y790964 = mcmean([1241,1264,1276.],[37, 23, 12.]) 

PortalesValley= mcmean([4477.,4458.],[11.,16.]) 

ALHA_L3_group = mcmean([3610.,3750.,3870.,3760.],[20,20,20,30.]) 

Rumuruti = mcmean([4470,4450.],[20,10.])

Trebbin = mcmean([4350.,4450.].-30,[20.,30.]) 

Marion = mcmean([4.38,4.36,4.40].*1e3, 20*ones(3)) 

Wellman = mcmean([4301, 4244,4232,4273,4342,4415,4589,4295,4309,4143,4139.],[42, 42, 42, 43,42,42,42,62,43,42,45,55.]) 

Chico = mcmean([550,550,590.],[50,20,20.]) 

Orvinio = mcmean([4050,4250,4170,4180,4350,4325,4400.],[160, 190, 130, 160, 50, 40, 60.])

Jilin = mcmean([3850.,4020.,3960.], [50.,160.,60.])