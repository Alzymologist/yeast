# Yeast uniformity test

Version 0.2.0

## Overview

This test is performed on stored slants with yeast to increase certainty that population is indeed monoclonal and samplings from this storage unit are indeed genetically identical.

## Prior knowledge

Prior probability that culture is monoclonal: Pmp = 99%

Mutation rate: km = 1e-7/s (P(mutation>1 | t) = 1 - exp(-t * km))

Probability of sterile transfer violation: Pf = 0.1%

## Equipment and materials

- Single cell process transfer tools
- Sterilized Petri dish with wort agar culture
- Constant temperature incubator
- Sterilized test tubes with identical amount (~15 ml) of same (~1.040) wort capped with foil x4
- Sterilized pasteur pipettes
- Cell counting set
- Sterile graduated pipettes x3
- Pipette driver
- Sterilized 50 ml vial with 35 ml of wort x3 - all from single stock
- Sterilized 150 ml vial with 100 ml wort capped with foil x3 - all from single stock
- Refractometer
- Organoleptic test vials set

## Procedure

1. Inoculate 1 test tube with slant culture, leave to incubate at room temperature until clear signs of activity are observed
2. Streak tube 1 on Petri dish; leave to incubate at room temperature until colonies could be picked up
3. Inoculate 3 test tubes with colonies from dish 2; leave to incubate at 21C in incubator until clear signs of activity are observed
4. At some points in time perform sterile sampling for cell counting from test tubes 3; count cells; calculate growth rates; at reliable step within exponential phase (25-50 MCFU/ml) calculate pitch rate to inject identical amounts of yeast into next stage, approximately 0.5 GCFU. Identify outliers - if there are any, stop experiment, analyse data.
5. Inoculate 50ml vials with measured volumes of contents 3; incubate vials in identical conditions until concentrations around 100MCFU/ml are observed. Calculate pitch rates to introduce 1GCFU to next stage.
5. Inoculate 150ml vials with measured volumes of contents 5; incubate vials in identical conditions until fermentation stops completely and then further for duration of diacetyl rest.
6. At arbitrary intervals, perform density measurements and cell counting for the vials 6 simultaneously. Organoleptically analyze the samples taken before disposal.
7. After total completion of fermentation, perform organoleptic analysis of samples
8. Later, analyze data.

## Analysis

TODO based on amount of first test results

There are several related but not interchangeable parameters here:
- probability that slant contains single colony: Pm
- probability that yeast sampled from slant is clone of yeast sampled at time of this test: Pr
- probability that two samplings result in cloned cultures: Pc

Parameter directly measured in this test is Pc. This parameter itself is not usable, as it degrades with time. Desired parameter is Pr. Pm could be used to transition from one into another?

## Resulting knowledge

Inferred probability of sampling reproducible colony from slant is probability that from any two samplings from 

Acceptable expectation of sampling reproducible colony from slant should be above 99%.

Although we already have this certainty of purity, thus this test only gives peace of mind in a sense and identifies obviously erroneous cultures. On the other hand, it might be giving roughly one more "9" in certainty level, depending on priors on how good we can identify divergence.

## Cleanup and safety

All general safety measures required for methods involved are mandatory.

## Notes

As the tests are don in controlled conditions after deprecation of 0.1.0, results could be used in other analysis datasets.

How to ensure identical aeration in vials 5?

We should not use density meter for density measurements of active inoculated material

We should pitch at sqrt(10) ratio because growing cell wall is another variable we do not want to introduce

Se should count cells on 15->50 ml transfer to ensure that 50->150 ml transfer adds equal volume of diluted (by fermentation) wort

Results of yeast counts are expected to yield 10-25% errors; this should not be used to adjust transfer volumes, but instead to reject obvious outliers. Thus, cell counting could be performed after sterile transfer.

With typical doubling time in fastest stage being 90 min, we would get 10% error within 12 min; error of 25% within 30 min. Lag times could be this long. Otherwise, pitching simultaneously would be more important that cell count.

Maybe 2-3 counts should be performed on first tubes to measure multiplication rate at exponential stage.

Cell count points without dilution should be discarded in analysis, as clumping makes measurements too noisy; those should be used to identify when dilution could be used. If yeast is not clumping, they might be used - we should note that in data files in processable form.