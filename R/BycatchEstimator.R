
#' Bycatch Estimator
#'
#' Generic Model-Based Bycatch Estimation Procedure
#'
#' The code can estimate both total bycatch, calculated by expanding a sample, such as an observer database, to total effort from logbooks or landings records, and an annual index of abundance, calculated only from the observer data.
#'
#' @docType package
#' @name BycatchEstimator
"_PACKAGE"

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
"logdatExample"

#' LLSIM-based example observer program data set
#'
#' LLSIm BUM by trip, with 5% observer coverage including observed catch in totals April 17 2022
#'
#' @format A tibble with 22 columns. Each row is a fishing trip.
#' \describe{
#' \item{Year}{Year}
#' \item{month}{month of the year}
#' \item{gear}{Gear code}
#' \item{light}{Light code}
#' \item{fleet}{Fleet number}
#' \item{bait}{Bait code}
#' \item{hook}
#' \item{hooks}
#' \item{sets}
#' \item{SWO}
#' \item{BUM}
#' \item{lat5}
#' \item{lon5}
#' \item{lat}
#' \item{lon}
#' \item{hbf}
#' \item{habSWO}
#' \item{habBUM}
#' \item{season}
#' \item{area}
#' }
#' @source Simulated data
"LLSIM_BUM_Example_observer"

#' LLSIM-based example logbook program data set
#'
#' LLSIm BUM by trip, with 5% observer coverage including observed catch in totals April 17 2022
#'
#' @format A tibble with 25 columns. Each row is a fishing trip.
#' \describe{
#' \item{trip}{}
#' \item{Year}{Year}
#' \item{month}{month of the year}
#' \item{gear}{Gear code}
#' \item{light}{}
#' \item{fleet}{Fleet number}
#' \item{bait}{}
#' \item{hook}
#' \item{hooks}
#' \item{sets}
#' \item{SWO}
#' \item{BUM}
#' \item{lat5}
#' \item{lon5}
#' \item{lat}
#' \item{lon}
#' \item{hbf}
#' \item{habSWO}
#' \item{habBUM}
#' \item{trip.05}
#' \item{trips}
#' \item{season}
#' \item{area}
#' \item{unsampledEffort}
#' }
#' @source Simulated data
"LLSIM_BUM_Example_logbook"
