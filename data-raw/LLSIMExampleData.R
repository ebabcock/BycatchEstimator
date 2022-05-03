


#------------------------------------------------------------------------------------------------
# LLSIm BUM by trip, with 5% observer coverage including observed catch in totals April 17 2022
#------------------------------------------------------------------------------------------------

LLSIM_BUM_Example_observer<-readr::read_csv("data-raw/obstrip05.csv", col_names = TRUE)
LLSIM_BUM_Example_logbook<-readr::read_csv("data-raw/logtrip05.csv", col_names = TRUE)

LLSIM_BUM_Example_observer<-LLSIM_BUM_Example_observer[LLSIM_BUM_Example_observer$Year>=1990,]
LLSIM_BUM_Example_logbook<-LLSIM_BUM_Example_logbook[LLSIM_BUM_Example_logbook$Year>=1990,]

usethis::use_data(LLSIM_BUM_Example_observer, overwrite = TRUE)
usethis::use_data(LLSIM_BUM_Example_logbook, overwrite = TRUE)
