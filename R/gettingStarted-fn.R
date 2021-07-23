###################################
###   Getting started   ###
###################################

##########################################################################
##  This function replaces values in a data.frame with new values
##  It allows for an overlap between the old and and the new valuses
##' @param table data.frame containing the data
##' @param column name of the column in which values need to be replaced
##' @param old_values old values that need to be replaced by new values
##' @param new_values new values that will replace the old values
##'        the ith entry of new:value replaces the ith entry of old_values
##' @return table data.frame where the values have been replaced
##' @author Lorenz Herger
##' 
replace_values <- function(table, column, old_values, new_values) {
  if (length(old_values) != length(new_values)) {
    stop("Vectors \"old_values\" and \"new_values\" are not of equal length.")
  }
  table$helper <- NA
  for (i in 1:length(old_values)) {
    table[table[[column]] == old_values[i] &!is.na(table[[column]]),]$helper <- new_values[i]
  }
  table[[column]] <- table$helper
  table$helper <- NULL
  table
}

##############################################################################
##  This function looks up the unique values of a binary column and swaps them
##' @param tb data.frame containing the data
##' @param col name of the column in which values need to be swapped
##' @return tb data.frame where the values have been replaced
##' @author Lorenz Herger
##' 
swap_binary_values <- function(tb, col) {
  trt <- unique(tb[[col]])
  trt <- trt[!is.na(trt)]
  if (! length(trt == 2)) stop("Treatment variable is not binary.")
  rev_trt <- rev(trt)
  tb <- replace_values(tb, col, trt, rev_trt)
}



####################################################################################
##  This function takes a data frame and its corresponding study name and makes adjustments necessary for our analysis
##' @param dat data frame of out study
##' @param name name of our study
##' @return data.frame containing data that is ready for analysis
##' @author Stefan Thoma (adapted from Lorenz Herger)

get_data_ready <- function(dat, name){
#dat <- dtml5
#name <- "alb5"
  # Remove SPSS-stlye attributes and classes

  dat[] <- lapply(dat, as.vector)
  dat <- as.data.frame(dat)

  
  # preprocess data
  if(name == "alb5"){
    dat <- dat %>% select(-starts_with("..."))
    dat <- try(replace_values(table = dat, column = "Condition", old_values = c("action", "inaction"), new_values = c(1, 0)))
    dat["Condition"] <- as_factor(dat["Condition"])
  }
  
  if(name=="payne"){
    dat <- dat %>% select(-ends_with("F")) %>% 
      rename(Location = Site,
             ResponseId = subject2)
    dat <- try(replace_values(table = dat, column = "Condition", old_values = c(-1, 1), new_values = c(0, 1)))
  }
  if(name == "lobue"){
    dat <- subset(dat, RT.correct <= 8.3 & number_errors < 8 &
                                      is.na(dat$snake_experience) == FALSE,
                                    select=Site:number_errors)
    dat["Location"] <- plyr::revalue(x=as.factor(dat$Site), c("1"="BG", "2"="NI", "3" = "LY", "4" = "NS"))
    dat["protocol"] <- plyr::revalue(x=as.factor(dat$protocol), c("1"="NP", "2"="RP"))
  }
  
  return(dat)
}



####################################################################################
##  This function takes the raw data and makes adjustments necessary for our analysis
## This function is specifically for the ManyLabs Project
##' @param dat the data.frame containing the raw data
##' @return dat_clean data.frame containing data that is ready for analysis
##' @author Lorenz Herger

get_data_ready_ml <- function(dat) {
  
  # Remove SPSS-stlye attributes and classes
  dat[] <- lapply(dat, as.vector)
  
  # The response variable allowedforbidden was incorrenctly encodend. Let's fix it:
  group_0 <- dat[dat$allowedforbiddenGroup == 0, ]
  dat[dat$allowedforbiddenGroup == 0, ] <- swap_binary_values(group_0, "allowedforbidden")
  
  # All effects presented in the original paper were positive. For some questions we have to swap
  # the tretments (usually given as 0 and 1) in order to get the postive sign
  swap_treatments <- c("sunkgroup", "gainlossgroup", "reciprocitygroup", "allowedforbiddenGroup")
  for (i in swap_treatments) dat <- swap_binary_values(dat, i)
  
  # The variable "sex" is encoded as "m" for male and "f" for female. We rencode "sex with "1" and "0"
  dat <- replace_values(dat, "sex", old_values = c("m", "f"), new_values = c(0, 1))
  
  # The response variable gainlossDV is binary with levels 1 and 2. We need to replace this by
  # 0 and 1 in order to be able to fit a logistic regression model
  dat <- replace_values(dat, "gainlossDV", old = c(1, 2), new = c(0, 1))
  
  # Change to factors:
  tofactor <- c("scales", "scalesgroup")
  dat[tofactor] <- as_factor(dat[tofactor])
  
  
  dat <- dat %>%
    rename(Location = referrer,
           ResponseId = session_id)
}
