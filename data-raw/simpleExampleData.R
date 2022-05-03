


#------------------------------------------------------
#Simple example from tab 1. from Babcock original code
#------------------------------------------------------

obsdatExample<-readr::read_csv("data-raw/exampleobs.csv", col_names = TRUE)
logdatExample<-readr::read_csv("data-raw/examplelog.csv", col_names = TRUE)

usethis::use_data(obsdatExample, overwrite = TRUE)
usethis::use_data(logdatExample, overwrite = TRUE)
