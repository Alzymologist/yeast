# Yeast uniformity test

Version 0.1.0 - deprecated

## Deprecation note

After running this test 4 times, we've found that due to variation in plate-to-tube pitch rates, cultures diverge up to 2 orders of magnitude. With daily variations of temperature and other temporal deviations, this diverges results too much. To counter this, we need to create conditions that are constant over time.

## Overview

This test is performed on stored slants with yeast to increase certainty that population is indeed monoclonal and samplings from this storage unit are indeed genetically identical.

## Prior knowledge

Prior probability that culture is monoclonal: Pmp = 99%

Mutation rate: km = 1e-7/s (P(mutation>1 | t) = 1 - exp(-t * km))

Probability of sterile transfer violation: Pf = 0.1%

## Equipment and materials

- Single cell process transfer tools
- Sterilized Petri dish with wort agar culture
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

1. Inoculate 1 test tube with slant culture, leave to incubate at room temperature until clear signs of activity are observed.
2. Streak tube 1 on Petri dish; leave to incubate at room temperature until colonies could be picked up.
3. Inoculate 3 test tubes with colonies from dish 2; leave to incubate at identical room temperature until clear signs of activity are observed.
4. At some points in time perform sterile sampling for cell counting from test tubes 3; count cells; calculate growth rates; at reliable step within exponential phase (25-50 MCFU/ml) calculate pitch rate to inject identical amounts of yeast into next stage, approximately 0.5 GCFU. Identify outliers - if there are any, stop experiment, analyse data.
5. Inoculate 50ml vials with measured volumes of contents 3; incubate vials in identical conditions until concentrations around 100MCFU/ml are observed. Calculate pitch rates to introduce 1GCFU to next stage.
5. Inoculate 150ml vials with measured volumes of contents 5; incubate vials in identical conditions until fermentation stops completely and then further for duration of diacetyl rest.
6. At arbitrary intervals, perform density measurements and cell counting for the vials 6 simultaneously. Organoleptically analyze the samples taken before disposal.
7. After total completion of fermentation, perform organoleptic analysis of samples.
8. Later, analyze data.

## Analysis

TODO based on amount of first test results.

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

Exact parameters of stages are not important, as we are not comparing different samples or experiments, but instead focus on comparing processes within this test. We should be focused on that only.

How to ensure identical aeration in vials 5?

We should not use density meter for density measurements of active inoculated material

We should pitch at sqrt(10) ratio because growing cell wall is another variable we do not want to introduce

We should stop at 150ml

Se should count cells on 15->50 ml transfer to ensure that 50->150 ml transfer adds equal volume of diluted (by fermentation) wort

[white, zeinasheff]
cells to pitch = 1MCFU * ml of wort * Plato = 1MCFU * 150 * 10 = 1.5 GCFU

at 50 ml, this is 30MCFU/ml

This is clearly significantly less than what we'll have.

Instead, we should just pitch equal number of CFU on each stage, choosing safely lowest.



Results of yeast counts are expected to yeld 10-25% errors; this should not be used to adjust transfer volumes, but instead to reject obvious outliers. Thus, cell counting could be performed after sterile transfer.

With typical doubling time in fastest stage being 90 min, we would get 10% error within 12 min; error of 25% within 30 min. Lag times could be this long. Otherwise, pitching simultaneously would be more important that cell count.

Maybe 2-3 counts should be performed on first tubes to measure multiplication rate at exponential stage.

