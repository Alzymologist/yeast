# Density analysis

Density is measured with U-tube for the starting wort and for the fermented wort
during the organoleptics testing.

## Starting wort

Starting wort must have density (when recalculated to 20C) in range 1.040 to
1.060, measured and recorded. It is assumed that the starting wort retains
density during final sterilization step prior to use.

In standard testing (15 mL -> 40 mL -> 100 mL sequence), all three vessels must
have identical wort.

## Fermented wort

Fermented wort density is measured during the organoleptics testing, in all
final vessels.

## Existing calculators

Online temperature correction calculators (such as
<https://www.brewersfriend.com/hydrometer-temp/>) use formula (see
<https://homebrewacademy.com/hydrometer-temperature-correction-calculator/> 
or <https://brewery.org/library/HydromCorr0992.html>) fitting the well-known
data for pure water (CRC, "Standard density of water").

Source: HBD #963, 9/7/92, by Christopher Lyons.

Note: available printed equations apparently have some typos, data is not
matching dramatically.

This data should not be used for fermented wort. It appears to be quite a
stretch even to use it for sugar solutions.

Temperature difference between the measured and the standard temperature is not
that drastic for typical fermented wort measurements (for example, 23.0 C vs
20.0 C), most likely the best course of action would be to increase the
error of an individual measurement.

Difference in measurements for pure water between 23.0 C and 20.0 C is
`6.651*10^(-4) g*cm^(-3)`, i.e. 0.067% of the density. Results spread between
different samples typically exceeds that greatly.

## Alcohol content

Simplest equation used for alcohol content is `ABV = (OG – FG) * 131.25`
(Balling equation). Balling equation assumes that (1) 1.00 g of ethanol is
produced consuming 0.11 g of carbohydrates and (2) only monosaccharides are
dissolved and thus affect the measured density (Cutaia, J. Am. Soc. Brew. Chem.
65(3):166-171, 2007).

Another source set: <https://www.meadmakr.com/how-to-determine-alcohol-by-volume/>
Apparently, coefficient has a range, 125..136.

All this is only good enough for an estimate.

Do we even need the alcohol content estimate? Possible reasons include real
attenuation calculation, and alcohol tolarance tests.

Hall estimate (Hall, 1995, Zymurgy magazine) could be better:

`ABW = 76.08 * (OG – FG) / (1.775 – OG)`

`ABV = ABW / 0.794, where 0.794 is the SG of ethanol`

## Attenuation

Everybody reports typical attenuation for yeast strains.

Apparent attenuation is simple to calculate, but is rather nonsensical
parameter, especially given that it can and regularly does go above 100%.

`AA = (FG - 1.0) / (OG - 1.0)`

Real attenuation seems more useful.

`E = –668.962 + 1262.45 SG – 776.43 SG^2 + 182.94 SG^3` extract equation, from
Hall (Hall, 1995, Zymurgy magazine).

Real extract calculation:
`q = 0.22 + 0.001 OE` attenuation coefficient
`RE = (q OE + AE) / (1 + q)`
(from Hall, after De Clerck)

Here OE is original extract (degrees Plato), calculated from original gravity,
AE is apparent extract, calculated from final gravity.

Equation for real attenuation:

`RA = 100% * (OE - RE) / OE`

Estimating uncertainties here would be a total pain in the neck.

## Protocols

All protocols must be recorded and versioned. When processing files must check
that a valid version is recorded.
Add protocol into organoleptic as well.
