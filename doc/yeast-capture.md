# Yeast capture statistical estimations

This is not a standard procedure; just general thoughts on priors estimation for a first slant.

## Starting media

Estimating yeast correspondence to particular expected strain on capture makes little sense; this would require us to trust manufacturer, to make assumptions on bottling process quality, some beers may have lagering yeast that is different from strain expected in style.

Laboratory yeast, of course, states some strain identification, but first, we capture commercial pure cultures from living beer fermentation, and then, whatever is listed by other laboratories is based on their less then open standards of processing and is thus below data quality we aspire to achieve here.

What we focus on here is yeast purity after capture, no matter how captured strain relates to origin.

## Initial capture

The yeast in starting media is often in poor condition and might be stuck to the walls. To pick it into solution, we've found the following to be effective:

1. Prepare approx. 50-200 ml of sterile approx. 1.040 pale wort.
2. Empty the container almost to the bottom, so that 10-50 ml remain.
3. Flame the opening and add sterile wort into container; mix thoroughly.
4. Wait until fermentation is visible; yeast could be sampled at high krauzen.
5. To sample, flame the opening of the bottle and pour some wort into sterile test tube. This wort could be used to streak the plate.

## Plate streak

The colonies on plate are observed as circles, often isolated, indicating effective single-cell distribution.

TODO: analyze plate photographs to obtain some information of colonies distribution

TODO: estimate how well we could identify overlapping circles

TODO: literature reference on this

PRIOR
Worst case estimate encountered in [Yeast: White, Zainacheff](ISBN-13: 978-0937381960) seems to be around 20% probability of pickup of colony from more than one origin cell. We will use worst-case for now.

Thus, plating from culture with probability of population being monoclonal Ps, we produce population with probability of being monoclonal Pf = Ps + (1-Ps) * 0.8

PRIOR
With worst-case scenario, Ps = 0 (capture from multiclonal source), Pf0 = 0.8

Just repeated streaking operations result in:

Pf1 = 0.8 + 0.2 * 0.8 = 0.96

Pf2 = 0.96 + 0.04 * 0.8 = 0.992

Having 99% probability of pure culture is our initial goal; process with 3 plating stages promises to result in this value if no other significant sources of contamination are introduced.

If divergence is observed on some stage, we should identify that stage as first streaking stage, as Ps = 1.

## Grow stages

Colonies picked up from plates are grown in tubes before new platings. This stage places additional selective pressure on non-monoclonal populations, but also intoduces probability of mutations.

PRIOR
TODO: Probability to contaminate test tube could be estimated as negligible for now (P = 0), as no events were identified yet.

PRIOR
TODO: Probability of mutation should be calculated from fission rate and mutation rate, that should be yet found. This happens occasionally, something that looks like petite mutants was observed in numerous experiments, especially where yeast that participated in long fermentations or storage was propagated.

## Population-based approach

Instead of computing probability of having pure culture, which makes error assessment difficult if possible at all, we should probaby use fraction of pedigree yeast in carrier as parameter.

It does not really matter what the distribution is in starting material on catch stage; any yeast can be chosen as "target" and we do not really care which at this point.

PRIOR
On plating stage, we can assume zero-truncated (since we definitely observe a colony) Poisson distribution of number of strains in a colony, with probability of 1 colony at 80%. This results in zero-truncated Poisson parameters:

P(1) = 0.8
0.8 * (e^l - 1) = l
l = 0.4

This is approximate value giving 81% probability of 1-colony pickup; let's re-define it as actual prior, since underlying value was also given with less than 1 significant figure.

Thus, after first pick-up, fraction of good cells in a tube is 1 with 81% probability, if it's 2 population, then fraction of good cells is uniformly distributed between 0.5 and 1 (since if it's less than 0.5, we name the other one "good"). And so forth.

