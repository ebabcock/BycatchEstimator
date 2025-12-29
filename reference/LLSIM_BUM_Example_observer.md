# LLSIM-based example observer program data set

LLSIm BUM by trip, with 5% observer coverage including observed catch

## Usage

``` r
LLSIM_BUM_Example_observer
```

## Format

A tibble with 22 columns. Each row is a fishing trip.

- trip:

  Trip idenifier

- Year:

  Year

- month:

  month of the year

- gear:

  Gear code

- light:

  Light code

- fleet:

  Fleet number

- bait:

  Bait code

- hook:

  Hook code

- hooks:

  Effort variable, hooks

- sets:

  Number of sets in the trip

- SWO:

  Swordfish catch in numbers

- BUM:

  Blue marlin catch in numbers

- lat5:

  Latitude assigned to 5 degree grid

- lon5:

  Longitude assigned to 5 degree grid

- lat:

  Trip latitude

- lon:

  Trip longitude

- hbf:

  Hooks between floats

- habSWO:

  SDM habitat variable for swordfish

- habBUM:

  SDM habitat variables for blue marlin

- season:

  Categorical season variable

- area:

  Categorical season variable

## Source

Simulated data
