## code to prepare `DATASET` dataset goes here
source("data-raw/gettingStarted-fn.R")

## ML1
ml1 <- haven::read_sav("data-raw/Data_ml1/CleanedDataset.sav")
ml1 <- get_data_ready_ml(ml1)

## ML5
# Alb 5
dtml5 <-  readxl::read_excel("data-raw/Data_alb5/ML5 Alb 5 Revised Protocol.xlsx")
dtmturk <- readxl::read_excel("data-raw/Data_alb5/ML5 Alb 5 RPP MTurk Protocol.xlsx")

dtml5 <- get_data_ready(dtml5, name = "alb5")
dtmturk <- get_data_ready(dtmturk, name = "alb5")


# Payne
path.temp <- "data-raw/Data_payne/"
file.temp <- "payne.csv"
if(!file.temp %in% list.files(path.temp)){
  # run script if datafile does not exist already.
  source(paste(path.temp, "ML5_Payne_Script.R", sep = ""), chdir = TRUE, local = new.env())
}
payne <- readr::read_csv(file = paste(path.temp, file.temp, sep = ""))
payne <- get_data_ready(payne, name = "payne")

# LoBue
lobue.dat <- haven::read_sav("data-raw/Data_LoBue/dataset_LoBue_raw.sav")
lobue.dat <- get_data_ready(dat = lobue.dat, name = "lobue")



usethis::use_data(ml1, dtml5, dtmturk, payne, lobue.dat, overwrite = TRUE)
