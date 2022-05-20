#' Observer simple example data set
#'
#' Example data set for logbook. Can be used to run bycatch estimation procedure along with observer example data set.
#'
#' @format A tibble with 5 columns.
#' \describe{
#' \item{EW}{East/West spatial grouping variable}
#' \item{season}{season temporal grouping variable}
#' \item{Year}{Year variable}
#' \item{sampled.sets}{The effort variable e.g., e.g. 1000 hook hours}
#' \item{Catch}{Bycatch species in numbers caught}
#' }
#' @source Simulated data
#' @keywords Simple simulated data set
"obsdatExample"

#' Logbook simple example data set
#'
#' Example data set for logbook. Can be used to run bycatch estimation procedure along with observer example data set.
#'
#' @format A tibble with 5 columns.
#' \describe{
#' \item{EW}{East/West spatial grouping variable}
#' \item{season}{season temporal grouping variable}
#' \item{Year}{Year variable}
#' \item{sets}{The effort variable e.g., e.g. 1000 hook hours}
#' \item{trips}{Number of sample units in the record, trips in this example. In this example logbook data is not aggregated, thus each row is 1 trip}
#' }
#' @source Simulated data
#' @keywords Simple simulated data set
"logdatExample"

#' LLSIM-based example observer program data set
#'
#' LLSIm BUM by trip, with 5% observer coverage including observed catch in totals April 17 2022
#'
#' @format A tibble with 22 columns. Each row is a fishing trip.
#' \describe{
#' \item{trip}{Trip idenifier}
#' \item{Year}{Year}
#' \item{month}{month of the year}
#' \item{gear}{Gear code}
#' \item{light}{Light code}
#' \item{fleet}{Fleet number}
#' \item{bait}{Bait code}
#' \item{hook}{Hook code}
#' \item{hooks}{Effort variable, hooks}
#' \item{sets}{Number of sets in the trip}
#' \item{SWO}{Swordfish catch in numbers}
#' \item{BUM}{Blue marlin catch in numbers}
#' \item{lat5}{Latitude assigned to 5 degree grid}
#' \item{lon5}{Longitude assigned to 5 degree grid}
#' \item{lat}{Trip latitude}
#' \item{lon}{Trip longitude}
#' \item{hbf}{Hooks between floats}
#' \item{habSWO}{SDM habitat variable for swordfish}
#' \item{habBUM}{SDM habitat variables for blue marlin}
#' \item{season}{Categorical season variable}
#' \item{area}{Categorical season variable}
#' }
#' @source Simulated data
#' @keywords Simulated data set from LLSIM
"LLSIM_BUM_Example_observer"

#' LLSIM-based example logbook program data set
#'
#' LLSIm BUM by trip, with 5% observer coverage including observed catch in totals April 17 2022
#'
#' @format A tibble with 25 columns. Each row is a fishing trip.
#' \describe{
#' \item{trip}{Trip idenifier}
#' \item{Year}{Year}
#' \item{month}{month of the year}
#' \item{gear}{Gear code}
#' \item{light}{Light code}
#' \item{fleet}{Fleet number}
#' \item{bait}{Bait code}
#' \item{hook}{Hook code}
#' \item{hooks}{Effort variable, hooks}
#' \item{sets}{Number of sets in the trip}
#' \item{SWO}{Swordfish catch in numbers}
#' \item{BUM}{Blue marlin catch in numbers}
#' \item{lat5}{Latitude assigned to 5 degree grid}
#' \item{lon5}{Longitude assigned to 5 degree grid}
#' \item{lat}{Trip latitude}
#' \item{lon}{Trip longitude}
#' \item{hbf}{Hooks between floats}
#' \item{habSWO}{SDM habitat variable for swordfish}
#' \item{habBUM}{SDM habitat variables for blue marlin}
#' \item{season}{Categorical season variable}
#' \item{area}{Categorical season variable}
#' \item{unsampledEffort}{Unsampled effort}
#' }
#' @source Simulated data
#' @keywords Simulated data set from LLSIM
"LLSIM_BUM_Example_logbook"
