# LifeTableFertility

**LifeTableFertility** is an R package that provides a Shiny application to construct
age-specific life tables and fertility schedules from individual female daily egg
records.
The package computes survival (`lx`), fertility (`mx`), reproduction (`lx·mx`),
and key demographic parameters including the net reproductive rate (R0),
mean generation time (T), intrinsic rate of increase (r), finite rate of increase (λ),
and doubling time (DT).
Confidence intervals for demographic parameters can be optionally estimated using
percentile bootstrap or delete-1 jackknife resampling at the female level.
## Installation
This package is intended to be installed from CRAN once released:
R
install.packages("LifeTableFertility")
Shiny application for life table and fertility analysis