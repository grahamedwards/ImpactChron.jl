## See description in DATA.md
using ImpactChron: Nrm, lNrm, Unf, mcmean, lognorm, lognormMC



##################################
##                              ##
##  Environmental Parameters    ##
##                              ##
##################################


## Age of oldest CAIs (Ma)
    # Pb-Pb age from Connelly+ 2012 (doi: 10.1126/science.1226919)
tss = Nrm(4567.3,.08)


## Initial solar system ²⁶Al/²⁷Al composition
    # From Jacobsen+ 2008 (doi: 10.1016/j.epsl.2008.05.003)
rAlo= Nrm(5.23e-5,0.065e-5)


## Midplane temperature at 2.5 au (K)
    # Compiled from Fig. 2 of Woolum & Casten 1999 (doi: 10.1111/j.1945-5100.1999.tb01408.x)
    # This data was only available as a histogram binned at 100 K steps from 0 to 700 K.
    # The following lines describe a cumulative counts function estimated from the histogram and calculates the lognormal distribution's parameters from it.
Tₘ = (  lnT = log.(50:100:650), # natural log of midpoints of midplane temperature bins
        n=(1,18,13,5,1,1,1) ) # counts in each bin
μlTₘ = sum(Tₘ.lnT .* Tₘ.n)/sum(Tₘ.n) # log-normal mean

Tm = lNrm( μlTₘ , sqrt(sum(Tₘ.n .* ( Tₘ.lnT .- μlTₘ).^2) / (sum(Tₘ.n)-1))) #log-normal mean and standard deviation.



##################################
##                              ##
##     Material Parameters      ##
##                              ##
##################################


## Specific Heat Capacity (J kg⁻¹ K⁻¹)
    # All data from Wach+ 2013 (ADS Bibcode: 2013M&PSA..76.5017W), a MetSoc abstract. 
    # This is the only published data on chondrite specific heat capacities for temperatures >600K.
    # We include measurements at T>300 K
cₚ = (; Sołtmany_373K=753., Sołtmany_473K=839., Sołtmany_573K=881., Sołtmany_673K=879., Sołtmany_773K=907., Sołtmany_823K=922., NWA4560_373K=821., NWA4560_473K=900., NWA4560_573K=894., NWA4560_673K=803., NWA4560_773K=747., NWA4560_823K=750. )

Cp = lognorm(cₚ)


## Thermal Conductivity (W m⁻¹ K⁻¹)
    # From compilation in Table 2 in Opeil+ 2012 (doi: 10.1111/j.1945-5100.2012.01331.x)
    # Includes previously published values of Opeil+ 2010 (doi: 10.1016/j.icarus.2010.01.021) and Yomogida & Matsui 1983 (doi: 10.1029/JB088iB11p09513)
    # Contains only H and L ordinary chondrites.
    # Excludes extreme outliers from Abee and Pillistfer following the advice of C. Opeil (2022, pers. comm.).

Opeil2012 = ( 1.25, 3.05, 0.82, 1.90, 0.45, 1.15, 3.15, 2.72, 2.26 )
Opeil2010 = ( 1.88, 1.47 )
YomogidaMatsui1983 = ( 3.53, 0.75, 3.60, 2.16, 2.35, 3.85, 1.54, 1.15, 0.55, 1.20, 0.73, 0.85, 2.31, 1.03, 2.14, 1.86, 0.4, 0.47, 1.54, 0.78, 1.24, 0.97 )

k = lognorm((Opeil2012...,Opeil2010...,YomogidaMatsui1983...))


## Bulk density (kg m⁻³)
    # All data compiled from measurements of Macke 2010 (PhD Thesis,  Identifier: CFE0003424, https://stars.library.ucf.edu/etd/1638)
    # Mean ordinary and enstatite chondrite densities calculated from falls only, as reported in Flynn+ 2018 (doi: 10.1016/j.chemer.2017.04.002)
    # No R (Rumuruti-type) falls were measured, so we calculated the mean of all RC finds with weathering grades < W3 (p. 274 of Macke 2010)

densities = ( H=Nrm(3350,10), L=Nrm(3300,10), LL=Nrm(3180,20), EH=Nrm(3580,50),EL=Nrm(3480,50), R=Nrm(mcmean([3320,3300,3280,3400,3150],[80,50,40,100,30])...) )

ρ = lognormMC(densities)


## Ar closure temperature (K)
    # Compiled from a number of Ar-Ar cosmochronology studies. 
    # We include bulk chondritic Tc's (calculated from Arrhenius plots of Ar degassing data) as well as an effective oligoclase Tc.
    # Typically these Tc's are reported as upper and lower bounds corresponding to upper and lower bound cooling rates.
    # If a study reports Tc's for multiple samples, we wrap those in a NamedTuple for that study.

Tc_olg=Nrm(550,20) # assumed Ar closure in oligoclase, used in Trieloff+ 2003 (doi: 10.1038/nature01499)
MIL05029=Nrm(313,20) # Weirich+ 2010 (doi: 10.1111/j.1945-5100.2010.01124.x)
Forrest= Unf(608,668) # calculated for 10¹-10³ ºC/Ma cooling rates for the K-bearing phase with higher-T Ar loss, from Bogard+ 2010 (doi: 10.1111/j.1945-5100.2010.01060.x)
BogardGarrison2009 = ( PV1=600., PV2=500. ) # Two aliquots of Portales Valley (doi: 10.1016/j.gca.2009.08.009)
Turner1978 = ( Tieschitz=Unf(313,333), Menow=Unf(563,593), Ochansk=Unf(423,443), Sena=Unf(353,373), ForestVale=Unf(593,623), Richardton=Unf(653,683), Butsura=Unf(613,633), Guarena=Unf(643,673), Kernouve=Unf(643,673), QueensMercy=Unf(363,383), MountBrowne=Unf(313,333), Saratov=Unf(413,433), Barwell=Unf(453,473), Shaw=Unf(543,573), Olivenza=Unf(423,443) ) # calculated for 1-10ºC/Ma cooling rates, from Turner+ 1978 (ADS Bibcode: 1978LPSC....9..989T)

Tc = lognormMC((Tc_olg,MIL05029,Forrest,BogardGarrison2009...,Turner1978...))



##################################
##                              ##
##    Asteroidal Parameters     ##
##                              ##
##################################


## Aluminium abundance 
    # Al mass fraction in Ordinary, Enstatite, and Rumuruti-type chondrites 
    # From Table 16.11 in Lodders & Fegley, 1998 (ISBN: 9780195116946)

Al_abundances = (H=0.0106, L=0.0116, LL=0.0118, EH=0.0082, EL=0.01, R=0.0106)

cAl=lognorm(Al_abundances)


## Accretion time & initial radius
    # Models that estimate either accretion time or initial parent body radius also report the other.
    # The exception to this is Sugiura & Fujiya 2014, which only reports accretion times for several parent bodies.
    # Notes after each study's data identify the corresponding meteorite group.
    # Individual numbers in `Tuple`s reflect discrete "accepted" simulations.
SugiuraFujiya2014 = (taEC = Nrm(1.83,0.1), taOC = Nrm(2.14,0.1), taRC=Nrm(2.1,0.1)) # EC, OC, RC parent bodies (doi: 10.1111/maps.12292)
Henke2013 = (ta = (2.026,1.843,1.789,1.836), R= (189.6e3,195.5e3,180e3,199.9e3) ) # H parent body (doi: 10.1016/j.icarus.2013.05.034)
Blackburn2017 = (R=125e3, ta=Unf(2.14,2.25)) # H,L parent bodies (doi: 10.1016/j.gca.2016.11.038)
GailTrieloff2019 = (ta=(1.89,1.835), R=(115e3,160e3)) # L parent body (doi: 10.1051/0004-6361/201936020)
EdwardsBlackburn2022 = (ta=Unf(1.995,2.265), R=150e3) # LL parent body (doi: 10.1126/sciadv.aay8641)
Trieloff2022 = (R=(147e3, 139e3, 139e3, 179e3, 232e3, 131e3), ta=(2.047,2.042,2.038,1.93,1.84,2.11)) # EL parent body (doi: 10.1016/j.icarus.2021.114762)
    
ta = lognormMC((GailTrieloff2019.ta,Trieloff2022.ta,Henke2013.ta,Blackburn2017.ta,EdwardsBlackburn2022.ta,SugiuraFujiya2014...))

R = lognormMC((GailTrieloff2019.R,Trieloff2022.R,Henke2013.R,Blackburn2017.R,EdwardsBlackburn2022.R))



############
############



compiled_priors = (; tss,rAlo,R,ta,cAl,Tm,Tc,ρ,Cp,k)
## Values:
# (tss = Nrm(4567.3, 0.08), rAlo = Nrm(5.23e-5, 6.5e-7), R = lNrm(11.92033413442048, 0.18588555315038457), ta = lNrm(0.6983039440737816, 0.07923560161086551), cAl = lNrm(-4.566499965801431, 0.1316240996751273), Tm = lNrm(5.351706355691695, 0.4690641196454104), Tc = lNrm(6.1914065556334315, 0.2570768077824802), ρ = lNrm(8.11989284050272, 0.0399592655014891), Cp = lNrm(6.732160240426698, 0.0789929082866023), k = lNrm(0.33192938784399095, 0.6335231873675615))