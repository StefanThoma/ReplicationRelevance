###################################
###  Functions for replication  ###
###################################

# Either SMD or OR
use_vars_list <- list(alb5 = c("SATTotal", "Condition", "SMD")

                      )

###########################################################################
## Creates a table showing the standardized mean difference obtained by
## each of the 36 universities that tried to replicate a specific study. 
## Note that the SMD calculated by this function contains a bias correction
## and therefore differs slightly from the values calculated by the authors
## of the original paper and by our function "smd_cohen".
## This function is used in the function "create_overview_tables"
##' @param dat data.frame containing the full data from the manyLabs experiments
##' @param treatment_var Name of the column in dat that contains the treatment 
##' assignments.
##' @param response_var Name of the column in dat that contains the response.
##' @param orig Statistic of original study
##' @return output data.frame with one row for each of the 36 replication attempts.
##' Each row of output contains the sample sizes and means for both treatment groups
##' as well as the the standardized mean difference.
##' @author Lorenz Herger & Stefan Thoma
escalc_inputs_smd <- function(dat, treatment_var, response_var, orig){
  #treatment_var <- x_var 
  #response_var <- y_var 
  dat <- dat[!is.na(dat[[treatment_var]]) &!is.na(dat[[response_var]]), ]
  levels_treatment <- sort(unique(dat[[treatment_var]]), decreasing = TRUE) # sort since escalc interprets 1
                                                                            # as treatmenr 2 as control
  labs <- unique(dat$Location)
  #input <- data.frame(n1 = orig$n.1, n2 = orig$n.2, m1 = numeric(), m2 = numeric(),
  #                    sd1 = numeric(), sd2 = numeric())
  input <- as.data.frame(orig)
  rownames(input) <- "original Study"

  for (i in labs) {
    #i <- labs[1]
    # Treatment 1 (treatment)
    repl_trt <- dat[dat$Location == i & dat[[treatment_var]] == levels_treatment[1], ][[response_var]]
    input[i,"n.1"] <- length(repl_trt)
    input[i, "m.1"] <- mean(repl_trt)
    input[i, "sd.1"] <- sd(repl_trt)
    # Treatment 2 (control)
    repl_trt <- dat[dat$Location == i & dat[[treatment_var]] == levels_treatment[2], ][[response_var]]
    input[i,"n.2"] <- length(repl_trt)
    input[i,"m.2"] <- mean(repl_trt)
    input[i,"sd.2"] <- sd(repl_trt)
  }
  output <- escalc(measure="SMD", m1i=m.1, sd1i=sd.1, n1i=n.1, m2i=m.2, sd2i=sd.2, n2i=n.2, data=input)
  output
}

########################################################################
## Creates a table showing the odds ratios obtained by each of the 36 
## universities that tried to replicate a specific study. 
## This function is used in the function "create_overview_tables"
##' @param dat data.frame containing the full data from the manyLabs experiments
##' @param treatment_var Name of the column in dat that contains the treatment 
##' assignments.
##' @param response_var Name of the column in dat that contains the response.
##' @return output data.frame with one row for each of the 36 replication attempts.
##' Each row of output contains the numbers of participants in each treatment group
##' that responded with 1 and 0 as well as the odds ratio..
##' @author Lorenz Herger
escalc_inputs_or <- function(dat, treatment_var, response_var){
  dat <- dat[!is.na(dat[[treatment_var]]) &!is.na(dat[[response_var]]), ]
  levels_treatment <- sort(unique(dat[[treatment_var]]), decreasing = TRUE) # sort since escalc interprets 1
                                                                            # as treatmetn 2 as control
  labs <- unique(dat$Location)
  input <- data.frame(n1p = numeric(), n2p = numeric(), n1n = numeric(), n2n = numeric())
  for (i in labs) {
    # Treatment 0 (control)
    repl_trt <- dat[dat$Location == i & dat[[treatment_var]] == levels_treatment[1], ][[response_var]]
    input[i, "n1p"] <- length(repl_trt[repl_trt == 1])
    input[i, "n1n"] <- length(repl_trt[repl_trt == 0])
    # Treatment 1 (treatment)
    repl_trt <- dat[dat$Location == i & dat[[treatment_var]] == levels_treatment[2], ][[response_var]]
    input[i, "n2p"] <- length(repl_trt[repl_trt == 1])
    input[i, "n2n"] <- length(repl_trt[repl_trt == 0])
  }
  output <- escalc("OR", ai=n1p, bi=n1n, ci=n2p, di=n2n, data=input)
}


##############################################################################
# This function calls the functions escalc_inputs_smd and escalc_inputs_or to 
# create a table with the effectsizes of all the 36 universites for a specific
# research question. In addition, the function calculates one-sided CIs for 
# the effect found by each university. Whether the effect size for a specific
# question is calculated as the standardized mean difference or as the odds
# ratio is determined in the argument "vars_for_question"
##' @param dat.list list containing the full data and other infos from the manyLabs experiment
##' @param coverage_probability Double. Coverage probability of the CIs. Defaults
##' to .95 
# @returns after_esclalc data.frame. Contains the university-level effect sizes
# and CIs for one research question.
##' @author Lorenz Herger
create_overview_tables <- function(dat.list, coverage_prob = .95, fix.ef = TRUE){
 #dat.list <- alb5
 #vars_for_question <- use_vars_list$alb5
 #
  if(!"heterogen" %in% names(dat.list)){fix.ef <- FALSE}
  
    study <- names(dat.list)[1]
    dat <- dat.list[[study]]
    vars_for_question <- variables_used_in_manyLabs[[study]]
    
    x_var <- vars_for_question[2]
    y_var <- vars_for_question[1]
    measure <- vars_for_question[3]
    # Get effect sizes and variance
    if (measure == "SMD") {
      after_escalc <- escalc_inputs_smd(dat, x_var, y_var, orig =  dat.list$orig)
    } else if (measure == "OR") { 
      after_escalc <- escalc_inputs_or(dat, x_var, y_var)
    }
    # Get one sided CIs for effect sizes
    effect <- after_escalc["yi"]
    sd <- sqrt(after_escalc["vi"])
    
    
    if(onesided){
      lower_bounds <- effect - qnorm(coverage_prob)*sd
      upper_bounds <- Inf
    } else{
      lower_bounds <- effect - qnorm(1-(1-coverage_prob)/2)*sd
      upper_bounds <- effect + qnorm(1-(1-coverage_prob)/2)*sd
    }
    
    after_escalc[, "CI_lower_bound"] <- lower_bounds
    after_escalc[, "CI_upper_bound"] <- upper_bounds
    
    if(fix.ef){
      stnd.fix <- stnd.beta.lmer(dat.list$heterogen$lrt[[3]][[1]])
      after_escalc[x_var, names(stnd.fix)] <-  stnd.fix
    }
    
    # Add resulting data.frame to the output list
    return(after_escalc)
}

####################################################################################r 
# This function takes a table containing university-level effectsizes and constructs
# a confidence interval for the difference in effectsizes found by a specific pair of 
# univerities. This function is used in the function "CI_pairwise_all"
##' @param pair Character verctor. Contains the names of the two universities whose
##' findings are to be compared.
##' @param input_table data.frame as created by the functions "escalc_inputs_SMD"
##' and "escalc_inputs_or"
##' @param coverage_prob Double. Coverage probability of the CIs. Defaults to .95
##' @returns CI. numeric vector containing the difference in effect sizes and the 
##' lower and upper boundaries of a confidence interval for this difference.
##' @author Lorenz Herger
CI_pairwise_one_pair <- function(pair, input_table, coverage_prob = .95) {
  pair <- as.vector(pair)
  quantile <- 1 - (1 - coverage_prob) / 2
  pair_data <- input_table[pair, ]
  difference <- pair_data[1, "yi"] - pair_data[2, "yi"]
  joint_sd <- sqrt(sum(pair_data[, "vi"]))
  lower <- difference - qnorm(quantile) * joint_sd
  upper <- difference + qnorm(quantile) * joint_sd
  CI <- c(diff = difference, CI_lower_bound = lower, CI_upper_bound = upper)
}


###############################################################################
# For each research question this function uses the function "CI_pairwise_one_pair" to 
# calculate confidence intervals for the pariwise differences in effect sizes for
# all pairs of universities..
##' @param effect_table data.frame as produced by the function 
##' "create_overview_tables". Contains the university level effect sizes for one
##' research question
##' @combinations data.frame. Optional. Contains all pairs of universities in the
##' data.  The data.frame Will be created within the funciton if it is not already
##' passed as an argument. Note that the universities are the same for all research
##' questions
##' @param coverage_prob Double. Coverage probability of the CIs that. will be 
##' constructed. Defaults to.95
##' @returns pairwise_diffs data.frame. Contains  CIs for all pairwise differences
##'  in effect sizes between universities for one specific research question
##' @author Lorenz Herger
CI_pairwise <- function(effect_table, combinations, cover_prob = .95) {
  labs <- row.names(effect_table)
  combinations <- combn(labs, 2)
  CIs <- lapply(as.data.frame(combinations), CI_pairwise_one_pair,effect_table, coverage_prob = cover_prob)
  CIs <- t(as.data.frame(CIs))
  # Add row.names to all tables so that we can see which pair of labs each row deals with
  pairnames <- paste(as.vector(combinations[1, ]), as.vector(combinations[2, ]), sep = "-")
  row.names(CIs) <- pairnames
  CIs <- as.data.frame(CIs)
}

#########################################################################
# This function creates a forestplot showing the confidence intervals for 
# all university-level effect sizes for one research question. A vertical
# line in each plot marks the benchmark for practical relevance.
##' @param table data.frame in the style of the entries of the list
##' produced by the function "create_overview_tables". Contains  all the 
##' university-level confidence intervals for a research question. 
##' @param measure character. Either takes value "SMD" for standardized
##' mean differnece or "OR" for odds ratio. The appropriate value for 
##' each research question can be found in use_vars_ist
##' @returns NULL
##' @author Lorenz Herger
plot_CI_vs_benchmark <- function(table, measure, main) {
 #table <- alb_table
 #measure <- "SMD"
  
    effects <- as.vector(table$yi)
    lower_bounds <- as.vector(table$CI_lower_bound)
    upper_bounds <- as.vector(table$CI_upper_bound)
    upper_bounds[upper_bounds>1e6] <- 1e6
    
    benchmark <- switch(measure, OR = log(1.25), SMD = .2)
    xlimit <- ceiling(max(effects)/.5) * .5
    
    
    
    # Upper bound will be given as 1e6. This is a hack we use because the forestplot 
    # function can not deal with Inf as an upper bound for a confidence interval
    forestplot(mean=effects, upper = upper_bounds, lower = lower_bounds, 
               labeltext=as.matrix(row.names(table)), clip = c(-1, xlimit),
               title = main, boxsize=.5, xlab=paste0("Effect Size (", measure, ")"),
               grid=benchmark, txt_gp=fpTxtGp(label =gpar(cex=.9), xlab = gpar(cex=.9), 
               ticks = gpar(cex=.8), title=gpar(cex=1.2)),
               xticks = c(-.5, 0, .5, 1, xlimit, round(benchmark, 2)))
}

##########################################################################################
# This function creates a forestplot showing confidence intervals for pairwise 
# differences in effectsizes. Since there are 630 pairs of univeritiesfor each research 
# question we restrict ourselfs to plotting the differences between the effect found by 
# one university and those found by all other universities.
##' @param table of data.frame as produced by the function "CI_pairwise". Contains the
##' confidence intevals for all pairwise university-level differences for one 
##' research question
##' @param measure character. Either takes value "SMD" for standardized mean differnece
##' or "OR" for odds ratio. The appropriate value for each research question can be found 
##' in use_vars_ist.
##' @param lab Name of the university for which the comparisons are to be plotted.
##' @param main Character. Title of the plot that will be produced
##' @returns NULL
##' @author Lorenz Herger
plot_CI_effect_diff <- function(table, measure, lab, main) {
  lab_entries <- table[grep(lab, row.names(table), value = TRUE), ]
  effects <- as.vector(lab_entries[, "diff"])
  lower_bounds <- as.vector(lab_entries[, "CI_lower_bound"])
  upper_bounds <- as.vector(lab_entries[, "CI_upper_bound"])
  lbl_txt <- gsub(paste0(lab, "-"), "", row.names(lab_entries))
  forestplot(mean=effects, upper = upper_bounds, lower = lower_bounds, 
             labeltext = lbl_txt,
             title = main, boxsize=.5, xlab=paste0("Effect Size (", measure, ")"),
             txt_gp=fpTxtGp(label =gpar(cex=.9), xlab = gpar(cex=.9), 
             ticks = gpar(cex=.8), title=gpar(cex=1.2)))
}

########################################################################################
## This function checks how many of the university-level results for a research question
## violate the three replication criteria significance, practical relevance and data
## compatibility.
##' @param eff_CI_table data.frame as contained in the list prodcued by the function
##' "create_overview_tables" Contains confidence intervals for all university-level effect
##' sizes
##' @param diff_CI_table data.frame as contained in the list prodcued by the function
##' "CI_pairwise_all" The data.frame contains confidence intevals for all pairwise university-
##' level differences for one research question.
##' @param benchmark double. Benchmark for practically relevant effect size.
##' @returns replicable data.frame with one row. Shows percentage of effects/pairwise comparisons
##' which violate the three replicability criteria.
##' @author Lorenz Herger
replication_one_question <- function(eff_CI_table, diff_CI_table, benchmark) {
  cover_0 <- mean(eff_CI_table$CI_lower_bound < 0)
  cover_benchmark <-  mean(eff_CI_table$CI_lower_bound < benchmark)
  diff_not_0 <- mean(sign(diff_CI_table$CI_lower_bound) == sign(diff_CI_table$CI_upper_bound))
  replicable <- round(data.frame(cover_0, cover_benchmark, diff_not_0), 2)
}

#####################################################################################
## This function uses the function "replication_one_question" to create a table which
## shows how many of the university-level results for each research question violate
## the three replication criteria significance, practical relevance and data
## compatibility.
##' @param eff_CI_list list of data.frames as produced by the function 
##' "create_overview_tables". Each table in the list contains the university-level
##' effect sizes for one reasearch question
##' @param diff_CI_list list of data.frames as produced by the function 
##' "CI_pairwise_all". Each data.frame contains the confidence intervals for
##'  all pairwise university-level differences for one research question. 
##' @param var_list list in the style of use_vars_list. Needed because it contains
##' information on the effect size measure used for each research question.
##' @returns rep_tb data.frame showing percentages of university-level results
##' violating significance, practical relevance and data compatibility.
##' @author Lorenz Herger
replication_table <- function(eff_CI_list, diff_CI_list, var_list) {
  rep_tb <- data.frame(numeric(), numeric(), numeric())
  for(i in names(eff_CI_list)) {
  measure <- var_list[[i]][3]
  bm <- switch(measure, OR = log(1.25), SMD = .2)
  new <- replication_one_question(eff_CI_table = eff_CI_list[[i]],
                                     diff_CI_table = diff_CI_list[[i]],
                                     benchmark = bm)
  rep_tb <- rbind(rep_tb, new)
  }
  names(rep_tb) <- c("non_sig", "non_sig_benchmark", "non_dat_comp")
  row.names(rep_tb) <- names(eff_CI_list)
  rep_tb <- rep_tb[order(rep_tb$non_sig), ]
}



###########################################################################
## Creates a table showing the standardized mean difference obtained by
## each of the 36 universities that tried to replicate a specific study. 
## This is an alternative function to escalc_inputs_smd and is based on the 
## relevance package by Stahel. 
## 
## 
##' @param dat data.frame containing the full data from the manyLabs experiments
##' @param treatment_var Name of the column in dat that contains the treatment 
##' assignments.
##' @param response_var Name of the column in dat that contains the response.
##' @param orig Statistic of original study
##' @return output data.frame with one row for each of the 36 replication attempts.
##' Each row of output contains the sample sizes and means for both treatment groups
##' as well as the the standardized mean difference.
##' @author  Stefan Thoma


#dat <- alb5@data.revised
#treatment_var <- alb5@variables$iv
#response_var <- alb5@variables$dv
#var.of.interest <- treatment_var
#orig <- alb5@original
relevance_smd <- function(dat, treatment_var, response_var, orig, var.of.interest){

  if(length(treatment_var)==1){
    dat <- dat[!is.na(dat[[treatment_var]]) &!is.na(dat[[response_var]]), ]
  }
  # as treatmenr 2 as control
  labs <- unique(dat$Location)
  if(length(orig)==6){
    
    
    #orig <-   list(m.1 = 12.83, m.2 = 10.78, sd.1 = 1.86, sd.2 = 3.15, n.1 = 18, n.2 = 18)
    
    # Get effect sizes of original study
    orig.est <- summary(escalc(measure="MD", m1i=orig$m.1, sd1i=orig$sd.1, n1i=orig$n.1, m2i=orig$m.2, sd2i=orig$sd.2, n2i=orig$n.2))
    # Get standardised effect sizes of original study
    orig.st.est <- summary(escalc(measure="SMD", m1i=orig$m.1, sd1i=orig$sd.1, n1i=orig$n.1, m2i=orig$m.2, sd2i=orig$sd.2, n2i=orig$n.2))
    
    estimate <- orig.est$yi[1]
    orig.vec <- round(with(orig.est, inference(estimate = yi[1], se = sei, df = orig$n.1+orig$n.2-1, 
                             stcoef = orig.st.est$yi[1])), digitsForRounding)
    
  }
  # Create result data frame
  output <- data.frame("estimate" = numeric(), "se" = numeric(), "statistic" = numeric(),
                       "p.value" = numeric(), "Sig0" = numeric(), "ciLow" = numeric(), 
                       "ciUp" = numeric(), "stcoef" = numeric(), "stciLow" = numeric(), 
                       "stciUp" = numeric(), "Rle" = numeric(), "Rls" = numeric(), "Rlp" = numeric()) 
  
  if(length(orig)==6){
    output["Original", names(orig.vec)] <- orig.vec
  } else{
    output["Original", names(orig)] <- orig
  }
  
  # Create formula based on input
  formula.temp <- as.formula(paste0(response_var, " ~ ", paste(treatment_var, collapse = " + "),  "| Location"))
  
  # Fit all models with this function from lmer
  model.list <- lmList(formula.temp, dat)
  # Get relevance inference results
  model.inference <- lapply(model.list, FUN = function(x)round(inference(x)[var.of.interest,], digitsForRounding))
  
  # Calculate power for each lab given the estimated effect size of the original study and the x vector of the replication study
  power.temp <- t(sapply(labs, FUN = function(l) power.f(x = as.vector(subset(dat, subset = Location == l)[[treatment_var]]), std.effect = output$stcoef[1], threshold = .1)))
  
  
  
  # Bind power-results and results to the output df
  output[labs, ] <- do.call("rbind", model.inference)
  output[labs, c("p.power", "r.power")] <- power.temp

  # Return output vector
  output
}






##############################################################################
# This function calls the functions relevance_smd to 
# create a table with the effectsizes of all the 36 universites for a specific
# research question. In addition, the function calculates one-sided CIs for 
# the effect found by each university. Whether the effect size for a specific
# question is calculated as the standardized mean difference or as the odds
# ratio is determined in the argument "vars_for_question"
##' @param dat.list list containing the full data and other infos from the manyLabs experiment
##' @param coverage_probability Double. Coverage probability of the CIs. Defaults
##' to .95 
# @returns after_esclalc data.frame. Contains the university-level effect sizes
# and CIs for one research question.
##' @author Lorenz Herger
relevance_overview_tables <- function(dat.list, coverage_prob = .95, fix.ef = FALSE){
#dat.list <- alb5
#vars_for_question <- use_vars_list$alb5
  
  if(!"heterogen" %in% names(dat.list)){fix.ef <- FALSE}
  
  study <- names(dat.list)[1]
  dat <- dat.list[[study]]
  vars_for_question <- variables_used_in_manyLabs[[study]]
  
  x_var <- vars_for_question[2]
  y_var <- vars_for_question[1]
  measure <- vars_for_question[3]
  # Get effect sizes and variance
  if (measure == "SMD") {
    after_escalc <- relevance_smd(dat, x_var, y_var, orig =  dat.list$orig)
  } else if (measure == "OR") { 
    stop("OR not defined yet")
    #after_escalc <- escalc_inputs_or(dat, x_var, y_var)
  }
  # Get one sided CIs for effect sizes
  if(fix.ef){
    stnd.fix <- stnd.beta.lmer(dat.list$heterogen$lrt[[3]][[1]])
    after_escalc["All", c("stcoef", "stciLow", "stciUp")] <-  stnd.fix
  }
  
  # Add resulting data.frame to the output list
  return(after_escalc)
}




##############################################################################
# This function plots the relevance_overview_tables
#
#
#
#
#
##' @param dat.list list containing the full data and other infos from the manyLabs experiment
##' @param coverage_probability Double. Coverage probability of the CIs. Defaults
##' to .95 
##' 
##' 
##' TODO: change measure
##' 
# @returns after_esclalc data.frame. Contains the university-level effect sizes
# and CIs for one research question.
##' @author Lorenz Herger
plot_relevance_table <- function(table, old = FALSE){

#table <- alb_table2
#dat_list <- alb5
  lab <- rownames(table)
  lab_entries <- table
  effects <- as.vector(lab_entries[, "stcoef"])
  lower_bounds <- as.vector(lab_entries[, "stciLow"])
  upper_bounds <- as.vector(lab_entries[, "stciUp"])
  
  # Get category of effects
  category <- success(table)
  
  
  if(old){
    
  forestplot(mean=effects, upper = upper_bounds, lower = lower_bounds,
             labeltext = matrix(paste(lab, category, sep = ", ")),
             title = main,
             #boxsize=.5,
             xlab=paste0("Effect Size (", "measure", ")"),
             txt_gp=fpTxtGp(label =gpar(cex=.9), xlab = gpar(cex=.9), 
                            ticks = gpar(cex=.8), title=gpar(cex=1.2)))
  
  }else{
    
  
  ## With ggplot (adapted from https://www.selfmindsociety.com/post/a-forest-plot-in-ggplot2)
  
  ## Have to add the name automatically
  xname <- "Alb5"
  
  ggtable <- table %>% rownames_to_column("Location") %>% 
    mutate(Location = factor(x = Location, levels = lab),
           category = success(table))
  
  #setting up the basic plot
ggplot(data=ggtable, aes(y=as.numeric(Location), x=stcoef, xmin=stciLow, xmax=stciUp, col = category))+ 
    
    #this adds the effect sizes to the plot
    geom_point() + 
    
    #this changes the features of the overall effects
    #one could resize, reshape, or recolor these points if desired
    geom_point(data=subset(ggtable, as.numeric(Location)==length(lab)), aes(y = as.numeric(Location), x=stcoef), color="Black", size=2) + 
    
    #adds the CIs
    geom_errorbarh(height=.1) +
    
    #sets the scales
    #note that I reverse the y axis to correctly order the effect #sizes based on my index variable
  #  scale_x_continuous(limits=c(-2.5,1), breaks = c(-2.5:1), name=xname)+
    scale_y_continuous(name = "", breaks=1:nrow(table), labels = paste(rownames(table)," (",category, ")",  sep = ""), trans="reverse")+
    
    #adding a vertical line at the effect = 0 mark
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)+
    #adding a vertical line at relevance threshold
    geom_vline(xintercept = threshold, alpha = .3, col = "red") +
    
    #faceting based on my subgroups
    #facet_grid(Location~., scales= "free", space="free")+
    
    #thematic stuff
    ggtitle("Target Effects")+
    theme_minimal()+
    theme(text=element_text(family="Times",size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines")) 
    
    
  
  }
  
  
}

##############################################################################
# This function takes as input the output of relevance_overview_tables
# It creates a new column that indicates whether a replication has been successful or not. 
#
##' @param table list containing the full data and other infos from the manyLabs experiment
##' to .95 
##' @param threshold relevance threshold. Needed if Rlp & Rls are not in table. 
##' @returns vector indicating whether repl. has been successful or not.
##' @author Stefan Thoma
##' 
success <- function(table, threshold = NULL){
  
  if(is.null(threshold)){
    
    vec <- with(table, ifelse(Rls>1, "Rlv",
                    ifelse(stciLow>0 & Rlp>1, "Amb.Sig",
                           ifelse(stciLow<0 & Rlp>1, "Amb",
                                  ifelse(stciLow<0 & Rlp<1, "Ngl",
                                         ifelse(stciUp<0, "Ctr", NA)
                                         )
                                  )
                           )
                    )
                )
    
    } else{
      
  # get threshold from options
    lo <- table$stciLow
    hi <- table$stciUp
    
    vec <- ifelse(lo>threshold, "Rlv",
                  ifelse(lo>0 & hi>threshold, "Amb.Sig",
                         ifelse(lo<0 & hi>threshold, "Amb",
                                ifelse(lo<0 & hi<threshold, "Ngl",
                                       ifelse(hi<0, "Ctr", NA)
                                )
                         )
                  )
    )
  }

}  




##############################################################################
# Calculates Power given a relvance table object. 
# 
# 
# 
# 
# 
##' @param dat.list list containing the full data and other infos from the manyLabs experiment
##' @param coverage_probability Double. Coverage probability of the CIs. Defaults
##' to .95 
# @returns power estimate based on simulation
##' @author Stefan Thoma

#rel.table <- try(scale@table)
#
#x <- rep(c(0,1), each = 40)
#std.effect <- .3
power.f <- function(x, std.effect, threshold, family ){

  if(family=="binomial"){
    z <- std.effect*x
    p <- 1/(1+exp(-z))
  }
  
  yf <- function(){
    if(family=="gaussian"){
      rnorm(mean = std.effect * x, n = length(x))
    } else if(family=="binomial"){
      rbinom(prob = p, size = 1, length(p))
    } else stop("family not defined")
  }
  
  
  sim <- function(){
    y <- yf()
    cl <- confint(profile(glm(y~x, family = family)), parm = "x")[1] # call profile first, so there is no message.
    return(c(cl > 0, cl > threshold))
#    return(coef(glm(y~x, family = family))) This line to get parameter estimates
  }
  
  res <- (do.call("rbind", replicate(sim(), n = 1000, simplify = FALSE)))
  power <- apply(res,MARGIN = 2, mean)
  names(power) <- c("p.power", "r.power")
  return(power)
  
}
#
#l <- labs[1]
#power.f(x = as.vector(subset(dat, subset = Location == l)[[treatment_var]]),
#        std.effect = output$stcoef[1], threshold = .1, family = family)
#t(sapply(labs, FUN = 
#           function(l) power.f(x = as.vector(subset(dat, subset = Location == l)[[treatment_var]]),
#                               std.effect = output$stcoef[1], threshold = .1, family = family)))
#