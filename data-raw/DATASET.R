set.seed(23532)



source("data-raw/gettingStarted-fn.R")
library(tidyverse)
## ML1
ml1 <- haven::read_sav("data-raw/Data_ml1/CleanedDataset.sav")
ml1 <- get_data_ready_ml(ml1)
scales <- ml1 %>% dplyr::select(ResponseId, Location, scales, scalesgroup)
scales.selection <- sample(unique(scales$Location), 10)
scales <- scales %>% filter(Location %in% scales.selection)


## ML5
# Alb 5
dtml5 <-  readxl::read_excel("data-raw/Data_alb5/ML5 Alb 5 Revised Protocol.xlsx")
dtmturk <- readxl::read_excel("data-raw/Data_alb5/ML5 Alb 5 RPP MTurk Protocol.xlsx")

dtml5 <- get_data_ready(dtml5, name = "alb5")
dtmturk <- get_data_ready(dtmturk, name = "alb5")

# remove unneccessary variables:
alb5_rev <- dtml5 %>% dplyr::select(ResponseId, Location, SATTotal, Condition)
alb5_rep <- dtmturk %>% dplyr::select(ResponseId, Location, SATTotal, Condition)


dtmturk <- dtmturk[sample(nrow(dtmturk), 150), ]


# Payne
path.temp <- "data-raw/Data_payne/"
file.temp <- "payne.csv"
if(!file.temp %in% list.files(path.temp)){
  # run script if datafile does not exist already.
  source(paste(path.temp, "ML5_Payne_Script.R", sep = ""), chdir = TRUE, local = new.env())
}
payne <- readr::read_csv(file = paste(path.temp, file.temp, sep = ""))
payne <- get_data_ready(payne, name = "payne")

payne <- payne %>% dplyr::select(ResponseId, Location, Indirect, Direct, Condition)


# LoBue
lobue <- haven::read_sav("data-raw/Data_LoBue/dataset_LoBue_raw.sav")
lobue <- get_data_ready(dat = lobue, name = "lobue")

lobue <- lobue.dat %>% dplyr::select(ID, Location, protocol, target_stimulus, child_parent, snake_experience, RT.correct)
length(lobue$Location)
usethis::use_data(scales, alb5_rev, alb5_rep, payne, lobue, overwrite = TRUE)
