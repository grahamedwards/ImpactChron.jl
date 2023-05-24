# Compiled data
The directory `/data/` contains the various compiled data used as priors in our Markov chain Monte Carlo inversions, including a database of [chondrite K-Ar and Ar-Ar ages](#database-of-chondrite-k-ar-and-ar-ar-ages) and a compilation of [constraints on parent body properties](#thermochronologic-model-parameters) used in the thermochronologic model. 

## Thermochronologic model parameters
`parameterpriors.jl` contains extensively annotated calculations of the distributions of the following thermochronologic model parameters:

- Environmental parameters of early solar system
  - solar system age (oldest Ca-Al-rich inclusions, CAIs)
  - initial solar system ²⁶Al/²⁷Al
  - solar system midplane temperature
- Material parameters of chondritic parent bodies
  - specific heat capacity
  - thermal conductivity
  - bulk density
  - Ar closure temperature
- Asteroidal parameters of parent bodies
  - aluminum abundance
  - body radius
  - accretion age (relative to the oldest CAIs)

The corresponding literature sources are identified by doi, Astrophysical Data System ([ADS](https://ads.harvard.edu/)) Bibcode, or another permanent identifier as necessary.

## Database of chondrite K-Ar and Ar-Ar ages

`KArArages_comipilation.csv` contains our database of chondritic K-Ar and Ar-Ar ages from the literature. The following table summarizes the columns in the database.

| Column | Description |
| :----- | :---------- |
| `name`   | meteorite name 
| `group` | genetic group/family (H, L, LL, EH, EL, or R)
| `type` | petrologic type ("im" indicates "impact melt")
| `ageMa` | K-Ar or Ar-Ar age in Ma
| `sigma` | 1σ standard deviation of age
| `lambda` | decay constants used to calculate age*
| `source` | literature source (see [References](#references))
| `srclocation` | location in literature source where age is reported.
| `notes` | any notes or additional information (e.g. denoting K-Ar ages)
||

*Abbreviations for decay constant references: `T69`= Turner 1969 (doi:[10.1007/978-94-010-3411-1](https://doi.org/10.1007/978-94-010-3411-1)), `H74` = Husain 1974 (doi: [10.1029/JB079i017p02588](https://doi.org/10.1029/JB079i017p02588)), `SJ77` = Steiger+Jäeger 1977 (doi: [10.1016/0012-821X(77)90060-7](https://doi.org/10.1016/0012-821X(77)90060-7))


### Calculating and recalibrating K-Ar/Ar-Ar ages

In cases where multiple viable ages were reported or only ages of individual heating steps were reported, we calculated a mean age using a Monte Carlo Method (`ImpactChron.mcmean`). These calculations are shown in `recalculateages.jl` and the results are reported in `KArArages_comipilation.csv`.

Since the ages in `KArArages_comipilation.csv` were originally calculated with different decay constants, we use a variety of techniques to recalibrate all ages to the near-ubiquitously used decay constants of Steiger+Jäeger 1977 (doi: [10.1016/0012-821X(77)90060-7](https://doi.org/10.1016/0012-821X(77)90060-7)). These calculations are summarized in `recalibrateages.jl`. In this script, we also divide the recalibrated ages into separate files containing only Ar-Ar ages (`ArArages.csv), only K-Ar ages (`KArages.csv`), or the entire recalibrated database (`KArArages.csv`).


### References:
Within the database csv file, literature sources of K-Ar and Ar-Ar dates are given as abbreviated reference codes, typically beginning with the first four letters of the first-author's name and the year of publication. The following table relates these codes to permanent identifiers such as a Digital Object Identifier (DOI) or Astrophysical Data System (ADS) Bibcode. If neither is available for a source, we provide a full citation.

---

| Reference code | DOI / ADS Bibcode / Citation | 
| :--- | :--- |
| Bene2008 | [10.1016/j.gca.2008.02.010](https://doi.org/10.1016/j.gca.2008.02.010) | 
| Boga1976 | [10.1029/JB081i032p05664](https://doi.org/10.1029/JB081i032p05664) | 
| Boga1980 | [10.1016/0016-7037(80)90219-7](https://doi.org/10.1016/0016-7037(80)90219-7) | 
| Boga1983 | [10.1016/0012-821X(83)90077-8](https://doi.org/10.1016/0012-821X(83)90077-8) | 
| Boga1987 | [10.1016/0016-7037(87)90192-X](https://doi.org/10.1016/0016-7037(87)90192-X) | 
| Boga1990 | [1990LPI....21..103B](https://ui.adsabs.harvard.edu/abs/1990LPI....21..103B) | 
| Boga1995 | [10.1111/j.1945-5100.1995.tb01124.x](https://doi.org/10.1111/j.1945-5100.1995.tb01124.x) | 
| Boga1995abs | [1995LPI....26..141B](https://ui.adsabs.harvard.edu/abs/1995LPI....26..141B) | 
| Boga1995chico | [10.1016/0016-7037(95)00051-Z](https://doi.org/10.1016/0016-7037(95)00051-Z) | 
| Boga2001 | [10.1111/j.1945-5100.2001.tb01813.x](https://doi.org/10.1111/j.1945-5100.2001.tb01813.x) | 
| Boga2009 | [10.1016/j.gca.2009.08.009](https://doi.org/10.1016/j.gca.2009.08.009) | 
| Boga2010 | [10.1111/j.1945-5100.2010.01060.x](https://doi.org/10.1111/j.1945-5100.2010.01060.x) | 
| Crab1981 | [10.1016/0016-7037(81)90097-1](https://doi.org/10.1016/0016-7037(81)90097-1) | 
| Dixo2003 | [10.1111/j.1945-5100.2003.tb00270.x](https://doi.org/10.1111/j.1945-5100.2003.tb00270.x) | 
| Dixo2004 | [10.1016/j.gca.2004.02.023](https://doi.org/10.1016/j.gca.2004.02.023) | 
| Eugs1988 | [10.1111/j.1945-5100.1988.tb00892.x](https://doi.org/10.1111/j.1945-5100.1988.tb00892.x) | 
| Folc2004 | [10.1016/j.gca.2003.11.023](https://doi.org/10.1016/j.gca.2003.11.023) | 
| Grie2004 | [10.1111/j.1945-5100.2004.tb00123.x](https://doi.org/10.1111/j.1945-5100.2004.tb00123.x) | 
| Hopp2014 | [10.1111/maps.12243](https://doi.org/10.1111/maps.12243) | 
| Jour2010 | [10.1016/j.gca.2009.11.032](https://doi.org/10.1016/j.gca.2009.11.032) | 
| Kane1979 |  I. Kaneoka, M. Ozima, M. Yanagisawa. 1979. <sup>40</sup>Ar-<sup>39</sup>Ar age studies of four Yamato-74 meteorites. <i>Memoirs of National Institute of Polar Research</i> (12): 186–206. | 
| Kane1980 |  I. Kaneoka. 1980.  <sup>40</sup>Ar-<sup>39</sup>Ar ages of L and LL chondrites from Alan Hills, Antarctica: ALHA77015, 77214 and 77304. <i>Memoirs of National Institute of Polar Research</i> (17): 177–188. | 
| Kane1981 | [1981PolRe......250K](https://ui.adsabs.harvard.edu/abs/1981PolRe......250K) | 
| Kane1988 |  I. Kaneoka, N. Takaoka, K. Yanai. 1988. <sup>40</sup>Ar-<sup>39</sup>Ar analyses of Yamato-75097 (L6) chondrite from Antarctica. <i>Proceedings of National Institute of Polar Research Symposium on Antarctic Meteorites</i> (1): 206–214.| 
| Keil1980 | [10.1016/0012-821X(80)90207-1](https://doi.org/10.1016/0012-821X(80)90207-1) | 
| Koro2007 | [10.1111/j.1945-5100.2007.tb00221.x](https://doi.org/10.1111/j.1945-5100.2007.tb00221.x) | 
| Koro2008 | [10.1016/j.gca.2008.05.014](https://doi.org/10.1016/j.gca.2008.05.014) | 
| Krin1996 | [10.1029/96JE03139](https://doi.org/10.1029/96JE03139) | 
| Kunz1997 | [10.1111/j.1945-5100.1997.tb01550.x](https://doi.org/10.1111/j.1945-5100.1997.tb01550.x) | 
| Mcco1988 | [10.1016/0016-7037(88)90307-9](https://doi.org/10.1016/0016-7037(88)90307-9) | 
| Mcco1995 | [10.1016/0016-7037(94)00231-A](https://doi.org/10.1016/0016-7037(94)00231-A) | 
| Metz2011 | [10.1111/j.1945-5100.2011.01181.x](https://doi.org/10.1111/j.1945-5100.2011.01181.x) | 
| Minh1984 | [1984LPI....15..552M](https://ui.adsabs.harvard.edu/abs/1984LPI....15..552M) | 
| Mull1985 | [1985LPI....16..595M](https://ui.adsabs.harvard.edu/abs/1985LPI....16..595M) | 
| Naka1994 | [1994AMR.....7..125N](https://ui.adsabs.harvard.edu/abs/1994AMR.....7..125N) | 
| Okan1990 | [10.1016/0016-7037(90)90301-Z](https://doi.org/10.1016/0016-7037(90)90301-Z) | 
| Pali2001 | [10.1111/j.1945-5100.2001.tb01958.x](https://doi.org/10.1111/j.1945-5100.2001.tb01958.x) | 
| Park2016 | [2016LPICo1921.6440P](https://ui.adsabs.harvard.edu/abs/2016LPICo1921.6440P) | 
| Podo1971 | [10.1016/0016-7037(71)90055-X](https://doi.org/10.1016/0016-7037(71)90055-X) | 
| Recc1986 | [10.1111/j.1945-5100.1986.tb01243.x](https://doi.org/10.1111/j.1945-5100.1986.tb01243.x) | 
| Rowe1965 | [10.1016/0016-7037(65)90001-3](https://doi.org/10.1016/0016-7037(65)90001-3) | 
| Rubi1981 | [10.1016/0016-7037(81)90073-9](https://doi.org/10.1016/0016-7037(81)90073-9) | 
| Rubi2011 | [10.1111/j.1945-5100.2011.01176.x](https://doi.org/10.1111/j.1945-5100.2011.01176.x) | 
| Ruzi2015 | [10.1016/j.gca.2015.04.030](https://doi.org/10.1016/j.gca.2015.04.030) | 
| Sche1998 | [10.1111/j.1945-5100.1998.tb01631.x](https://doi.org/10.1111/j.1945-5100.1998.tb01631.x) | 
| Schu1972 | [10.1016/0012-821X(72)90039-8](https://doi.org/10.1016/0012-821X(72)90039-8) | 
| Schu1977 | [10.1016/0012-821X(77)90061-9](https://doi.org/10.1016/0012-821X(77)90061-9) | 
| Srin1977 | [10.1016/0016-7037(77)90157-0](https://doi.org/10.1016/0016-7037(77)90157-0) | 
| Step1988 | [10.1111/j.1945-5100.1988.tb00926.x](https://doi.org/10.1111/j.1945-5100.1988.tb00926.x) | 
| Swin2009 | [10.1111/j.1945-5100.2009.tb00766.x](https://doi.org/10.1111/j.1945-5100.2009.tb00766.x) | 
| Swin2011 | [2011LPI....42.1897S](https://ui.adsabs.harvard.edu/abs/2011LPI....42.1897S) | 
| Swin2014 | [10.1144/SP378.6](https://doi.org/10.1144/SP378.6) | 
| Taka1989 | [10.1515/zna-1989-1006](https://doi.org/10.1515/zna-1989-1006) | 
| Taki1987 |  Y. Takigami & I. Kaneoka. 1987. Investigation of the effect of shock on the Antarctic meteorites by the <sup>40</sup>Ar-<sup>39</sup>Ar method. <i>Memoirs of National Institute of Polar Research</i> (46): 133–143. | 
| Trie1989 | [1989LPICo.712..243T](https://ui.adsabs.harvard.edu/abs/1989LPICo.712..243T) | 
| Trie1994 | [1994Metic..29Q.541T](https://ui.adsabs.harvard.edu/abs/1994Metic..29Q.541T) | 
| Trie2003 | [10.1038/nature01499](https://doi.org/10.1038/nature01499) | 
| Trie2007 | [2007M&PSA..42.5059T](https://ui.adsabs.harvard.edu/abs/2007M&PSA..42.5059T) | 
| Turn1969 | [10.1007/978-94-010-3411-1](https://doi.org/10.1007/978-94-010-3411-1) | 
| Turn1973 | [1973Metic...8..447T](https://ui.adsabs.harvard.edu/abs/1973Metic...8..447T) | 
| Turn1978 | [1978LPSC....9..989T](https://ui.adsabs.harvard.edu/abs/1978LPSC....9..989T) | 
| Turr2023 | [10.1111/maps.13953](https://doi.org/10.1111/maps.13953) | 
| Udry2019 | [10.1111/maps.13252](https://doi.org/10.1111/maps.13252) | 
| Weir2009 | [2009M&PSA..72.5368W](https://ui.adsabs.harvard.edu/abs/2009M&PSA..72.5368W) | 
| Weir2010 | [10.1111/j.1945-5100.2010.01124.x](https://doi.org/10.1111/j.1945-5100.2010.01124.x) | 
| Weir2012 | [10.1111/j.1945-5100.2012.01397.x](https://doi.org/10.1111/j.1945-5100.2012.01397.x) | 
| Welt2003 | [10.1111/j.1945-5100.2003.tb01052.x](https://doi.org/10.1111/j.1945-5100.2003.tb01052.x) | 
| Whit2000 | [10.1126/science.288.5472.1819](https://doi.org/10.1126/science.288.5472.1819) | 
| Witt2011 | [10.1016/j.gca.2011.07.037](https://doi.org/10.1016/j.gca.2011.07.037) | 
---


