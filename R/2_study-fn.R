###################################-
###    Study Class Functions    ###-
###################################-

# define test study object:
#try(object <- alb5)
#
##object <- try(sunk)
#try(object <- payne)
#try(object <- scales)
#
# get one function from web
source("http://pcwww.liv.ac.uk/~william/R/crosstab.r")


# Is slot empty? ------------------------
##' Function throws error if specified slot is empty
##' @param object Object of class study
##' @param slot Name of slot to check.
##' @return no return value.
##' @author Stefan Thoma

check_slot.study <- function(object, slot){
  if(rlang::is_empty(slot(object, slot))){
    stop(paste("The slot ", slot, "is not defined in the study object" ))
  }
}
check_slot.study <- Vectorize(check_slot.study, vectorize.args = "slot", SIMPLIFY = TRUE)

# Function which summarises the data of a study object. ------------------------
##' summary.study, function to summarise the data frames of the study objects.
##' result depends on the measure OR or SMD defined in the study object.
##' @param object Object of class study
##' @return some summary measures of the data, mainly about the data
##' @author Stefan Thoma
##' @export

summary.study <- function(object){
  if(!isClass(object, Class = "study")){
    stop("input object is not of class study")
  }

  dv <- object@variables$dv
  iv <- object@variables$iv
  return.list <- list()

  # general function
  if(object@family=="gaussian"){

    sum.f <- function(df){
      df %>% dplyr::group_by(Location) %>%
        dplyr::summarise("N" = length(unique(ResponseId)),
                  dv_mean = mean(get(dv), na.rm = TRUE),
                  dv_sd =  sd(get(dv), na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
  } else if(object@family == "binomial"){
    # define the summary function for the binomial case
    sum.f <- function(df){
      if(!exists(x = "crosstab")){
          aggregate(data = df, get(dv)~get(iv)+Location, FUN = table)
      # the crosstab function is imported from: source("http://pcwww.liv.ac.uk/~william/R/crosstab.r")
      } else {crosstab(df, row.vars = c( "Location", dv), col.vars = iv, type = "f", subtotals=FALSE)}

    }
    } else{stop("family is not specified as either binomial or gaussian")}



  if(object@manyLabs == "ml5"){
    if(object@family == "binomial"){
      return(list("revised" = sum.f(df1), "replication" = sum.f(df2)))
    }
    # extract first both data frames
    df1 <- object@data.revised
    df2 <- object@data.replication
    summ.temp1 <- sum.f(df1)
    summ.temp2 <- sum.f(df2)

    cols.temp <- colnames(summ.temp1)
    summ.temp1["type"] <- "revised"
    summ.temp2["type"] <- "replication"
    summ <- rbind(summ.temp1, summ.temp2)
    # summarise both dataframes in one table
    summ <- summ %>% dplyr::select(type, dplyr::all_of(cols.temp))

    return.list[["summary"]] <- summ
  }
  if(object@manyLabs == "ml1"){
    if(object@family == "binomial"){
      return(sum.f(df2))
    }

    df2 <- object@data.replication
    summ.temp2 <- sum.f(df2)

    cols.temp <- colnames(summ.temp2)
    summ.temp2["type"] <- "replication"
    summ <- rbind(summ.temp2)

    summ <- summ %>% dplyr::select(type, dplyr::all_of(cols.temp))

    return.list[["summary"]] <- summ
  }

    return(return.list)

}

# Set method
setMethod("summary", signature = "study", summary.study)



# MixedModels.study --------------------------------------------------------------------------------
##'  Function which fits two mixed models and a fixed effect model and returns objects containing them.
##' @param object takes an object of class study as input
##' @param methodFamily which family does the error distribution belong to?
##' not required, can be used to overwrite the default defined in the object.
##' @param transformation name of the transformation that should be applied to the response (e.g. log).
##' @return list with three objects of fitted models. The first is the "empty model", i.e. a model without random effects
##' but still with the fixed treatment effect. The second is the model with the random interaction. The third is the model
##' with an interaction between treatment and (random) university effect.
##' @author Stefan Thoma, adapted from Federico Rogai
##' @export

MixedModels.study <-function(object, methodFamily=NULL, transformation=NULL){
  if(!isClass(object, Class = "study")){
    stop("input object is not of class study")
  }

  # Get name of x and y variables
  yVar<-object@variables$dv
  xVar<-object@variables$iv
  base_iv <- paste(xVar, collapse = " + ")

  if(object@manyLabs=="ml5"){
    dat <- object@data.revised
  } else{dat <- object@data.replication}


  if(is.null(methodFamily)){
    methodFamily <- object@family
    }

  # Get transformation function according to input
  if(!is.null(transformation)) {
    FUN <- match.fun(transformation)
    dat[, yVar]<-FUN(dat[,yVar])
  }


  message <- ""
  sing <- FALSE

  if(tolower(methodFamily)!="gaussian"){
    fit.empty<-glm(formula=as.formula(paste0(yVar," ~", base_iv)), dat, family=methodFamily, na.action = na.omit)
    fit0 <- lme4::glmer(formula=as.formula(paste0(yVar," ~", base_iv, "+ (1|Location)")), dat, family=methodFamily)
    fit1 <- lme4::glmer(formula=as.formula(paste0(yVar," ~", base_iv, "+ (1+", xVar[1] ,"|Location)")), dat, family=methodFamily)
    sing <- lme4::isSingular(fit1, tol = 1e-5)
    if(sing){
      fit1 <- lme4::glmer(formula= as.formula(form <- paste0(yVar," ~", base_iv, " + (1 | Location) + ( - 1 +" , xVar[1] ," | Location)")), dat, family = methodFamily)
      sing <- lme4::isSingular(fit1, tol = 1e-5)
      message <- paste(("complex model is singular. Alternative model was fitted: "), form)
      if(sing){
        fit1 <- lme4::glmer(formula=as.formula(form2 <- paste0(yVar," ~", base_iv, " + ( - 1 +" , xVar[1] ," | Location)")), dat, family = methodFamily)
        sing <- lme4::isSingular(fit1, tol = 1e-5)
        message <- paste("model ", form, " is singular. Alternative model was fitted: ", form2)
      }
    }
  } else{
    fit.empty<-aov(formula=as.formula(paste0(yVar," ~", paste(base_iv, collapse = " + "))), dat, na.action=na.omit)
    fit0 <- lme4::lmer(formula=as.formula(paste0(yVar," ~", base_iv, "+ (1|Location)")), dat)
    fit1 <- lme4::lmer(formula=as.formula(form <- paste0(yVar," ~", base_iv, "+ (1+", xVar[1] ,"|Location)")), dat)
    sing <- lme4::isSingular(fit1, tol = 1e-5)
    # delete intercept-slope covariance if model was singular
    if(sing){
      fit1 <- lme4::lmer(formula=as.formula(form <- paste0(yVar," ~", base_iv, " + (1 | Location) + ( - 1 +" , xVar[1] ," | Location)")), dat)
      sing <- lme4::isSingular(fit1, tol = 1e-5)
      message <- paste("Model ", form, " is singular. Alternative model was fitted: ", form, ". Singularity:", sing)

    }

  }


  # get intraclass correlation (icc)
  icc.names <- c(names(lme4::ranef(fit1, condVar = TRUE)$Location))
  icc <- as.data.frame(specr::icc_specs(fit1))

  # Name the grp variable a bit nicer
  # Works for most cases I think
  if(nrow(icc) == length(icc.names)+1){
    icc$grp <- c(icc.names, "Residual")

  } else if(nrow(icc) == length(icc.names)+2){
    icc$grp <- c(icc.names, "Covariance", "Residual")
  }

  # Get Rle based on a 10% Rel. Threshold
  icc$Rle <- icc$percent/10

    # Round values
  icc[-1] <- round(icc[-1], digits = digitsForRounding)



  models <- list("NullModel" = fit.empty,
                 "RandomIntercept" = fit0,
                 "RandomCoefficients" = fit1)
  return.list <- list("Models" = models,
                      "RandomCoefficients.singular" = sing,
                      "icc" = icc)

  print(message)
  return(return.list)
}



#   LRT.study ------------------------------------------------------------

##'  Function computes a likelihood ratio test (LRT) between the model with and without random intercept. If the p-value for this test is
##'  smaller than 0.05, we compute the LRT between the model with the random slope and the one without this term (but with the
##'  random intercept). If the larger model delivers a superior fit to the data, we print the summary of that model. Otherwise we
##'  report the summary of the model we had chosen before.
##'  Note: to compute the LRT we need to have the ML estimate. anova() automatically refits the models,
##'  whose parameters were estimated using REML
##' @param object is an object of class study, which contains an output of object mixedModels
##' @return It returns a list with five elements.
##' The first element is the model we have chosen (model with 1-merely fixed effects, 2- random intercept only,
##' 3- with random intercept and random slope)
##' The second element reports the p-value of the hypothesis we have tested in order to choose the model
##' The third element the orginal fit for the "winning model" as from MixedModels
##' The fourth returned argument is the summary of the chosen model
##' The fifth element are the confidence interval of the winning model
##' @author Stefan Thoma, adapted from Federico Rogai
##' @export


LRT.study <-function(object){

  # check if object is of correct class
  if(!isClass(object, Class = "study")){
    stop("input object is not of class study")
  }


  if(rlang::is_empty(object@mixedModels)){
    warning("the object does not contain mixed models. The function mixedModels.study is called and results are used")
    object@mixedModels <- MixedModels.study(object)
  }


  mod.temp <- object@mixedModels$Models

  pval0<-anova(mod.temp$RandomIntercept, mod.temp$NullModel, test = "LRT")[["Pr(>Chisq)"]][[2]]
  pval1<-anova(mod.temp$RandomIntercept, mod.temp$RandomCoefficients, test = "LRT")[["Pr(>Chisq)"]][[2]]
  if(pval0<0.05 & pval1<0.05){
    testResult<- cat("Random interaction needed! \n Larger model fits data better than
                       purely fixed effect \n and random intercept only models. \n
                       p-values of LRT",round(pval0,digitsForRounding), " \n
                       respectively",round(pval1,digitsForRounding), ".")
    finalModel<-mod.temp$RandomCoefficients
  }
  else{
    pval2<-anova(modelsToCompare[[2]], modelsToCompare[[1]], test = "LRT")[["Pr(>Chisq)"]][[2]]
    if(pval1 > 0.05 & pval2 < 0.05){
      testResult<- cat("Random intercept only is the preferred model. \n
                       p-value of LRT",round(pval1,digitsForRounding), " comparing it to the random slope model \n
                       respectively",round(pval2,digitsForRounding)," to the fixed effects only model")
      finalModel<-mod.temp$RandomIntercept
    }
    if(pval2 > 0.05 & pval0 > 0.05){
      testResult <- cat("Pure fixed effect model is the preferred model. \n
                       p-value of LRT",round(pval0,digitsForRounding), " comparing it to the random slope model \n
                       respectively",round(pval2,digitsForRounding)," to the random intercept only model")
      finalModel <- mod.temp$NullModel
    } else{ cat("WARNING: contraddictory evidence from the LRTs")}
  }
  winningModel <- summary(finalModel)

  list("finalModel" = finalModel, "testResult" = testResult, modelsToCompare[finalModel], list("Summary of best model", winningModel,confint(modelsToCompare[[finalModel]], method="Wald")))

}



#  relevance_table.study ------------------------------------------------------------
##'  creates a table with the effectsizes of all locations for a specific
##'  research question. In addition, the function calculates one-sided CIs for
##'  the effect found by each university. Whether the effect size for a specific
##'  question is calculated as the standardized mean difference or as the odds
##'  ratio is determined in the argument "vars_for_question"
##'
##' @param object study object
##' @param coverage_probability Double. Coverage probability of the CIs. Defaults
##' to .95
##' @param mem boolean, whether fixed effect from mixed effects model should be included
##' @param use.both  boolean, should both data frames be used?
##' @returns after_esclalc data.frame. Contains the university-level effect sizes
##' and CIs for one research question.
##' @author Stefan Thoma, adapted from Lorenz Herger
##' @export
##'
#object <- scales.study
relevance_table.study <- function(object,
                                        coverage_prob = .95,
                                        mem = TRUE,
                                        use.both = TRUE){


  # check if we need drop effect, then an entirely other function is called.
  if(object@variables$measure=="drop"){
    return(r_squared.study(object))
  }



  if(rlang::is_empty(object@original)){
    original <- NULL
  } else(original <- object@original)
  if(rlang::is_empty(object@mixedModels)){
    warning("the object does not contain mixed models. The function mixedModels.study is called and results are used")
    object@mixedModels <- MixedModels.study(object)
  }


  # extract information
  study <-object@name
  y_var <- object@variables$dv
  x_var <- object@variables$iv
  measure <- object@variables$measure


  # depending on which many labs project we look at, we save either one or two df in a list
  if(object@manyLabs=="ml5"){
    if(use.both){
      dat <- list(object@data.revised, object@data.replication)
    } else{dat <- list(object@data.revised)}
  } else{
    dat <- list(object@data.replication)
  }


  #relevance_f.study(dat = dat[[1]], x_var, y_var, orig =  original, var.of.interest = object@var.of.interest, family = object@family)
  # Get effect sizes and variance
  after_escalc <- lapply(dat, FUN = function(x){
    relevance_f.study(dat = x, x_var, y_var, orig =  original, var.of.interest = object@var.of.interest, family = object@family, relevance.threshold = object@relevance.threshold)})

  #if(object@manyLabs=="ml5" & use.both){}
  # bind the list output from previous function together
  after_escalc <- do.call("rbind", after_escalc)

  # do the mixed effects modeling effect size extraction and standardisation
  if(mem){
    stnd.fix <- stnd.beta.glmer(mod = object@mixedModels$Models$RandomCoefficients)[object@var.of.interest,]
    after_escalc["All", c(names(stnd.fix),
                          "n.total",
                          "Rle",
                          "Rls",
                          "Rlp")] <-  c(stnd.fix,                                           # estimate and CI's
                                        nobs(object@mixedModels$Models$RandomCoefficients), # n.total
                                        stnd.fix["stcoef"]/object@relevance.threshold,      # Rle
                                        stnd.fix["stciLow"]/object@relevance.threshold,     # Rls
                                        stnd.fix["stciUp"]/object@relevance.threshold)      # Rlp
    }



  #
  after_escalc["type"] <- NA
  if(object@manyLabs=="ml5" & use.both){

    rn.revision <- unique(object@data.revised$Location)
    rn.replication <- unique(object@data.replication$Location)
    if(anyDuplicated(c(rn.revision, rn.replication))){
      stop("some locations in both df's have the same name")
    }

    after_escalc[ rn.revision,  "type"] <- "revision"
    after_escalc[ rn.replication,  "type"] <- "replication"

    after_escalc["Original1",] <- NA
  } else if(object@manyLabs=="ml5" & !use.both){
    rn.revision <- unique(object@data.revised$Location)
    after_escalc[ rn.revision,  "type"] <- "revision"
  } else{
    rn.replication <- unique(object@data.replication$Location)
    after_escalc[ rn.replication,  "type"] <- "replication"
  }



  after_escalc["type"] <- ifelse(rownames(after_escalc)=="Original", "original",
                                 ifelse(rownames(after_escalc)=="All", "fixef.rand.coef",
                                        after_escalc$type))
  # Add resulting data.frame to the output list
  after_escalc <- as.data.frame(after_escalc)
  return(after_escalc)

  }














# Plot diagnostics ----------------------------------
##'  Wrapper function for performance::check_model().
##'  Used to be a custom funciton, but the check_model function just seems so much better.
##'  Why reinvent the wheel?
##' @param object study object
##' @param ... forwarded to performance::check_model()
##' @return output of performance::check_model()
##' @author Stefan Thoma, adapted from Federico Rogai
##' @export
diagnosticPlot.study<-function(object, ...){

  if(rlang::is_empty(object@mixedModels)){
    warning("the object does not contain mixed models. The function mixedModels.study is called and results are used")
    object@mixedModels <- MixedModels.study(object)
  }

  model<-object@mixedModels$Models$RandomCoefficients
  performance::check_model(model, ...)

}


##  plotting function ------------------------------------------------------------
##' This function plots the relevance_overview_tables
##' @param object a study object
##' @param coverage_probability Double. Coverage probability of the CIs. Defaults
##' to .95
##'
##'
##' TODO: change measure
##'
# @returns after_esclalc data.frame. Contains the university-level effect sizes
# and CIs for one research question.
##' @author Stefan Thoma
##' @export

plot_estimates.study <- function(object, cutoff = NULL, standardise = TRUE){
  if(!is.null(cutoff) & length(cutoff)==1){cutoff <- c(-cutoff, cutoff)}
  if(rlang::is_empty(object@table)){
    warning("the object does not contain the table. The function mixedModels.study is called and results are used")
      object@table <- relevance_table.study(object)}

  if(object@variables$measure %in% c("OR", "drop")){
    standardise <- FALSE
  }


  table <- object@table

  table <- table[!is.na(table$type), ]

  lab <- rownames(table)
# effects <- as.vector(table[, "stcoef"])
# lower_bounds <- as.vector(table[, "stciLow"])
# upper_bounds <- as.vector(table[, "stciUp"])


  # if we want non-standardised effect sizes, we just replace the standardised ones.
  if(!standardise){
    table <- table %>% dplyr::mutate(
      stciLow = as.numeric(ciLow),
      stciUp = as.numeric(ciUp),
      stcoef = as.numeric(estimate)
    )
  }

  # Get category of effects

  category <- success(table, threshold = object@relevance.threshold)

  xname <- object@name

  ggtable <- table %>% tibble::rownames_to_column("Location") %>%
    dplyr::mutate(Location = factor(x = Location, levels = lab),
           category = success(table, threshold = object@relevance.threshold),
           alp = ifelse(type =="replication" & !rlang::is_empty(object@data.revised), .5, 1))

  # helping function from: https://stackoverflow.com/questions/45857787/adding-column-if-it-does-not-exist
  fncols <- function(data, cname) {
    add <-cname[!cname%in%names(data)]

    if(length(add)!=0) data[add] <- as.numeric(NA)
    data
  }
# add the two columns needed (if they don't exist yet) so the plotting function does not throw an error.
  ggtable <- fncols(ggtable, c("stpredUp", "stpredLow"))
  #setting up the basic plot
  if(is.null(cutoff)){
    plot.table.study(ggtable, category = category, threshold = object@relevance.threshold)+
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=stpredLow, xmax=stpredUp), height=.1, size = .4)+
      ggplot2::xlab(expression(theta))
  #  ggtitle(paste("Target Effects for study ", xname, " of project ", object@manyLabs, sep = ""))
  } else {
    plot.table.study(ggtable, category = category, threshold = object@relevance.threshold) +
     ggplot2::coord_cartesian(xlim =cutoff) +
     ggplot2::geom_errorbarh(ggplot2::aes(xmin=stpredLow, xmax=stpredUp), height=.1, size = .4)+
     ggplot2::xlab(expression(theta))
   # ggtitle(paste("Target Effects for study ", xname, " of project ", object@manyLabs, sep = ""))
  }

}

#  difference of estimates plotting function ------------------------------------------------------------
##' This function plots the difference table as a forest plot
##' Adapted from https://www.selfmindsociety.com/post/a-forest-plot-in-ggplot2
##' @param object a study object
##' @param coverage_probability Coverage probability of the CIs. Defaults
##' @param ci.type either "wald" or "newcombe"
##' to .95
##' @returns plot
##' @author Stefan Thoma
##' @export

plot_difference.study <- function(object, cutoff = NULL, ci.type = "wald", standardise = TRUE){

  # Deal with the cutoff for plotting: If only one is supplied, apply it to both sides.
  if(!is.null(cutoff) & length(cutoff) == 1){cutoff <- c(-cutoff, cutoff)}
  if(rlang::is_empty(object@difference.table)){
    warning("the object does not contain the table. The function mixedModels.study is called and results are used")
    object@difference.table <- effect_differences.study(object)
    }


 # if(object@variables$measure=="OR"){
 #   standardise <- FALSE
 # }


# Extract information
  table <- object@difference.table[[ci.type]]

# Filter out unwanted rows
  table <- table[!startsWith(rownames(table), prefix = "Original1"),]
  table <- table[!startsWith(rownames(table), prefix = "Diff.: NA"),]
# Get names of locations
  lab <- rownames(table)

  # if no stcoef in names, then just use non-std-values
  # The same if we do not standardise!
  if(!"stcoef" %in% names(table) | !standardise){
    table <- table %>% dplyr::mutate(
      stciLow = as.numeric(ciLow),
      stciUp = as.numeric(ciUp),
      stcoef = as.numeric(estimate)
    )
  }



  # Get category of effects
  category <- success(table, threshold = object@relevance.threshold)



  ## With ggplot (adapted from https://www.selfmindsociety.com/post/a-forest-plot-in-ggplot2)

  xname <- object@name

  ggtable <- table %>% tibble::rownames_to_column("Location") %>%
    dplyr::mutate(Location = factor(x = Location, levels = lab),
           category = success(table, threshold = object@relevance.threshold),
           alp = 1,
           threshold = object@relevance.threshold)




  #setting up the basic plot
  if(is.null(cutoff)){
    plot.table.study(ggtable, category = category, threshold = object@relevance.threshold) +
    ggplot2::geom_vline(xintercept = -object@relevance.threshold, alpha = .3, col = "red") +
    # ggtitle(paste("Effect differences of study ", xname, " of project ", object@manyLabs, sep = "")) +
    ggplot2::xlab(expression(Delta))
  } else {
    plot.table.study(ggtable, category = category, threshold = object@relevance.threshold) +
    ggplot2::coord_cartesian(xlim =cutoff) +
    ggplot2::geom_vline(xintercept = -object@relevance.threshold, alpha = .3, col = "red")+
    ggplot2::xlab(expression(Delta))
     # ggtitle(paste("Effect differences of study ", xname, " of project ", object@manyLabs, sep = ""))
  }

}

# plot ggtable  ---------------------------------------------------
##' Function used to plot estimates and differences in estimates
##' @param ggtable either from difference or estimate table
##' @param category the categorisation of the results for color choice
##' @return plot of studies
##' @author Stefan Thoma, adapted from Federico Rogai

plot.table.study <- function(ggtable, category, threshold){

  ## define colors for all possible success possibilities:
  cls <- RColorBrewer::brewer.pal(6, "Set1")
  sccs <- c( "Rlv" =  cls[1],  "Amb.Sig" =  cls[2],  "Amb" =  cls[3], "Ngl.Sig" = cls[4], "Ngl" =  cls[5], "Ctr" =  cls[6])

  #as.numeric(ggtable$stciLow)<as.numeric(ggtable$stcoef)
  ggplot2::ggplot(data=ggtable, ggplot2::aes(y=as.numeric(Location), x=stcoef, xmin=stciLow, xmax=stciUp, col = category))+

    #this adds the effect sizes to the plot
    ggplot2::geom_point(alpha = ggtable$alp, size = ggtable$alp) +
    #adds the CIs
    ggplot2::geom_errorbarh(height=.1, alpha = ggtable$alp, size = ggtable$alp) +

    #this changes the features of the overall effects
    #one could resize, reshape, or recolor these points if desired
    ggplot2::geom_point(data=subset(ggtable, as.numeric(Location)==length(Location)),ggplot2::aes(y = as.numeric(Location), x=stcoef)) +


    #sets the scales
    #note that I reverse the y axis to correctly order the effect #sizes based on my index variable

    #ggplot2::scale_y_continuous(name = "", breaks=1:nrow(ggtable), labels = paste(ggtable$Location," (",category, ")",  sep = ""), trans="reverse")+
    ggplot2::scale_y_continuous(name = "", breaks=1:nrow(ggtable), labels = paste(ggtable$Location,  sep = ""), trans="reverse")+

    #adding a vertical line at the effect = 0 mark
    ggplot2::geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)+
    #adding a vertical line at relevance threshold
    ggplot2::geom_vline(xintercept = threshold, alpha = .3, col = "red") +

    #faceting based on my subgroups
    #facet_grid(type~., space="free")+

    #thematic stuff

    ggplot2::theme_minimal()+
    ggplot2::theme(text=ggplot2::element_text(family="Times",size=18, color="black"))+
    ggplot2::theme(panel.spacing = ggplot2::unit(1, "lines")) +
    ggplot2::scale_color_manual(values = sccs)
}




# create relevance table ---------------------------------------------------
##' Creates a table showing the standardized mean difference obtained by
##' each of the 36 universities that tried to replicate a specific study.
##' This is an alternative function to escalc_inputs_smd and is based on the
##' relevance package by Stahel. The function is called by the relevance_table.study function
##'
##'
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
##relevance_smd.study(dat, treatment_var = x_var, response_var = y_var, orig =  object@original, var.of.interest = object@var.of.interest, family = object@family)
#
#var.of.interest = object@var.of.interest
#family = object@family
#orig =  object@original
#response_var = object@variables$dv
#treatment_var = object@variables$iv
#dat <- dat[[1]]
relevance_f.study <- function(dat, treatment_var, response_var, orig, var.of.interest, family = "gaussian", relevance.threshold = .1){

  # remove rows with NA
  if(length(treatment_var)==1){
    dat <- dat[!is.na(dat[[treatment_var]]) &!is.na(dat[[response_var]]), ]
  }
  labs <- unique(dat$Location)
  # Do stuff with the original results, if orig vector is supplied.
  # Create result data frame
  output <- data.frame("estimate" = numeric(), "se" = numeric(), "statistic" = numeric(),
                       "p.value" = numeric(), "Sig0" = numeric(), "ciLow" = numeric(),
                       "ciUp" = numeric(), "stcoef" = numeric(), "stciLow" = numeric(),
                       "stciUp" = numeric(), "Rle" = numeric(), "Rls" = numeric(), "Rlp" = numeric(), n.total = numeric())

  if(!is.null(orig)){


    if(length(orig)==6){


      #orig <-   list(m.1 = 12.83, m.2 = 10.78, sd.1 = 1.86, sd.2 = 3.15, n.1 = 18, n.2 = 18)

      # Get effect sizes of original study
      orig.est <- summary(metafor::escalc(measure="MD", m1i=orig$m.1, sd1i=orig$sd.1, n1i=orig$n.1, m2i=orig$m.2, sd2i=orig$sd.2, n2i=orig$n.2))
      # Get standardised effect sizes of original study
      orig.st.est <- summary(metafor::escalc(measure="SMD", m1i=orig$m.1, sd1i=orig$sd.1, n1i=orig$n.1, m2i=orig$m.2, sd2i=orig$sd.2, n2i=orig$n.2))

      estimate <- orig.est$yi[1]
      orig.vec <- round(with(orig.est, relevance::inference(estimate = yi[1], se = sei, df = orig$n.1+orig$n.2-1,
                                                 stcoef = orig.st.est$yi[1])), digitsForRounding)
      orig.vec["n.total"] <- sum(orig$n.1, orig$n.2)

      output["Original", names(orig.vec)] <- orig.vec

    }else{
        output["Original", names(orig)] <- orig
      }



  }

  unfortunate.standardisation <- function(df){
    df <- df %>%
      dplyr::mutate(
        stcoef = estimate,
        stciLow = ciLow,
        stciUp = ciUp,
        Rle = estimate/ relevance.threshold,
        Rls = ciLow/ relevance.threshold,
        Rlp = ciUp/ relevance.threshold
      )
    df
  }

  # Check if predictor is not numeric.
  if(unf.std <- length(unique(dat[[treatment_var]]))<3 & family == "binomial"){
    output <- unfortunate.standardisation(output)
  }


  # Create formula based on input
  formula.temp <- as.formula(paste0(response_var, " ~ ", paste(treatment_var, collapse = " + "),  "| Location"))

  # Fit all models with this function from lmer
  model.list <- lme4::lmList(formula.temp, as.data.frame(dat), family = family)
  # Get relevance inference results
  model.inference <- lapply(model.list, FUN = function(x){
    return.value <- round(inference(x)[var.of.interest,], digitsForRounding)
    return.value["n.total"] <- nobs(x)
    return.value

    })

  if(!is.null(orig)){
    # Calculate power for each lab given the estimated effect size of the original study and the x vector of the replication study
    power.temp <- t(sapply(labs, FUN =
                             function(l) power.f(x = as.vector(subset(dat, subset = Location == l)[[treatment_var]]),
                                                 std.effect = output$stcoef[1], threshold = .1, family = family)))
  }

  # Bind power-results and results to the output dfs

  output[labs, ] <- do.call("rbind", model.inference)
  if(!is.null(orig)){
    output[labs, c("p.power", "r.power")] <- power.temp
  }
  # Return output vector
  if(unf.std){
    output <- unfortunate.standardisation(output)
  }
  output
}






# create difference table ------------------------------------------------------
##' WIP
##' Calculates two difference tables with two distinct methods to calculate CI.
##' wald: Uses pooled SD to create normal CI
##' newcombe: Uses newcombe method to create CI
##' @param object Study object.
##' @param coverage_prob Double. Coverage probability of the CIs. Defaults to .95
##' @returns CI. numeric vector containing the difference in effect sizes and the
##' lower and upper boundaries of a confidence interval for this difference.
##' @author Stefan Thoma, adapted from Lorenz Herger
##' @export
#object <- scales
effect_differences.study <- function(object, coverage_prob = .95, standardise = NULL) {
  # check if original study effect or table is supplied.
  check_slot.study(object, c("original", "table"))


  ## standardise?
  # IF family is gaussian I standardise, if not I don't.
  if(is.null(standardise)){
    if(object@variables$measure=="SMD"){
      standardise <- TRUE
    } else {standardise <- FALSE}
  }

  # extract information
  table <- object@table
  labs <- rownames(table)[-1]
  original <- table["Original",]
  quantile <- 1 - (1 - coverage_prob) / 2
  replications <- table[-1,]

  # function to calculate pooled sd
  pool.sd <- function(n, s){
    sqrt(
      {sum({n-1}*{s^2})} / {sum(n)-length(n)}

    )
  }

  # define function to calculate differences with CI
  diff <- function(lab){
    difference <- - as.numeric(lab[["estimate"]]) + as.numeric(original[["estimate"]])
    joint_se <- sqrt(sum(c(as.numeric(original[["se"]]), as.numeric(lab[["se"]]))^2))
    lower <- difference - qnorm(quantile) * joint_se
    upper <- difference + qnorm(quantile) * joint_se
    CI <-data.frame("estimate" = difference, "ciLow" = lower, "ciUp" = upper)
    return(CI)
  }


  # alternative method by newcombe
  diff.newcombe <- function(lab){

    difference <- - {est2 <- as.numeric(lab[["estimate"]])} + {est1 <- as.numeric(original[["estimate"]])}

    l1 <- as.numeric(original[["ciLow"]])
    l2 <- as.numeric(lab[["ciLow"]] )
    u1 <- as.numeric(original[["ciUp"]])
    u2 <- as.numeric(lab[["ciUp"]] )

    lower <- difference - sqrt((est1-l1)^2 + (u2-est2)^2)
    upper <- difference + sqrt((est2-l2)^2 + (u1-est1)^2)
    CI <-data.frame("estimate" = difference, "ciLow" = lower, "ciUp" = upper)



    return(CI)
  }


  ## standardise
  standardise.ed <- function(lab){
#    print(lab)
#lab <- replications[1,]
    if(standardise & object@variables$measure== "SMD"){
      sd.original <- original$estimate/original$stcoef
      sd.lab <- lab["estimate"]/lab["stcoef"]
      pooled.sd <- pool.sd(n = c(original$n.total, lab["n.total"]),
                           s = c(sd.original, sd.lab))
    } else if(!standardise){
      pooled.sd  <-  1
    }

  return(pooled.sd)
  }

  standardise.factors <- apply(replications[,!names(replications) %in% "type"], 1, FUN = standardise.ed, simplify = TRUE)


  # call function diff
  diff.wald.list <- apply((replications), 1, FUN = diff, simplify = TRUE)
  diff.newcombe.list <- apply((replications), 1, FUN = diff.newcombe, simplify = TRUE)
  # bind results
  diff.wald.table <- do.call("rbind", diff.wald.list)
  diff.newcombe.table <- do.call("rbind", diff.newcombe.list)
  # add factor
  diff.wald.table["factor"] <- standardise.factors
  diff.newcombe.table["factor"] <- standardise.factors
  # standardise
  diff.wald.table[c("stcoef", "stciLow", "stciUp")] <- diff.wald.table[c("estimate", "ciLow", "ciUp")]/standardise.factors
  diff.newcombe.table[c("stcoef", "stciLow", "stciUp")] <- diff.newcombe.table[c("estimate", "ciLow", "ciUp")]/standardise.factors


  # add relevance measures:
  diff.wald.table[c("Rle", "Rls", "Rlp")]<- diff.wald.table[c("stcoef", "stciLow", "stciUp")] / object@relevance.threshold
  diff.newcombe.table[c("Rle", "Rls", "Rlp")]<- diff.newcombe.table[c("stcoef", "stciLow", "stciUp")] / object@relevance.threshold

  return(list("wald" = diff.wald.table, "newcombe" = diff.newcombe.table))
}





# get CI for R^2 ------------------------------------------------------
##' We want to get the conditional R^2 with confidence interval separately for each
##' lab but also for the MEMo.
##'
##' Conditional R^2 refers to:
##' Conditional R_GLMM² is interpreted as a variance explained by the entire model,
##' including both fixed and random effects,
##' and is calculated according to the equation:
##' R_GLMM(c)² = (σ_f² + σ_α²) / (σ_f² + σ_α² + σ_ε²)
##'
##' This is opposed to the marginal r^2, see MuMIn package.
##'
##'
##' @param object Study object.
##' @param coverage_prob double.
##' @returns CI. numeric vector containing estimate and CI
##' @author Stefan Thoma
#object <-  lobue.study
#dat <- lobue.study@data.revised
r_squared.study <- function(object,  coverage_prob = .95, bootstrap.type = "norm"){

  # prepare formulas as string
  vars <- object@variables
  formula1 <- paste(vars$dv, "~", vars$iv)
  formula2 <- paste(vars$dv, "~", object@var.of.interest)

  # setup return df
  output <- data.frame(estimate = numeric(),
                       se = numeric(),
                       ciLow = numeric(),
                       ciUp = numeric(),
                       n.total = numeric(),
                       type = character()
                      #Rle = numeric(),
                      #Rls = numeric(),
                      #Rlp = numeric()
                         )




  # create bootstrap functions
  ## first for glm
  stat.f.lm <- function(dat, i){
    dat2 <- dat[i, ]
    m1 <- lm(as.formula(formula1), data  = dat2)
    m2 <- lm(as.formula(formula2), data = dat2)
    #MuMIn::r.squaredGLMM(m1, m2)[1,"R2c"]
    #attr(MuMIn::r.squaredLR(m2, m1), "adj.r.squared")
    rsq::rsq(m1, adj = TRUE) - rsq::rsq(m2, adj = TRUE)
  }

  ## now for MEMo
  stat.f.MEMo <- function(dat, i){

    dat2 <- dat[i, ]
    m1 <- lme4::lmer(as.formula(paste(formula1, "(1 | Location)", sep = " + ")), data  = dat2)
    m2 <- lme4::lmer(as.formula(paste(formula2, "(1 | Location)", sep = " + ")), data = dat2)
    #MuMIn::r.squaredGLMM(m1, m2)[1,"R2c"]
    rsq::rsq.lmm(m1, adj = TRUE)$fixed - rsq::rsq.lmm(m2, adj = TRUE)$fixed

    #attr(MuMIn::r.squaredLR(m2, m1), "adj.r.squared")
  }

  # create function to calculate statistic for all labs
  # This function calls stat.f.lm
  for.all.labs <- function(dat){
    #dat <- object@data.revised
    out.list <- by(data = dat, INDICES = dat$Location, FUN = function(x){

        boot_o <- boot::boot(data = x, statistic = stat.f.lm, R = 1000)
        boot_ci <- boot::boot.ci(boot_o, conf = coverage_prob, type = bootstrap.type)
        out <- c(boot_ci$t0, sd(boot_o$t), tail(c(boot_ci[[4]]), 2), nrow(x))
        names(out) <- names(output)[-6]
        out
      })
    tibble::as_tibble(do.call(rbind, out.list))

  }




  # Bootstrap for each type of manylab project
  if(object@manyLabs=="ml5"){
    # MEmo
    # Estimate & Bootstrap CI
    MEMo_boot <- boot::boot(data = object@data.revised, statistic = stat.f.MEMo, R = 1000)
    MEMo_boot_ci <- boot::boot.ci(MEMo_boot, conf = coverage_prob, type = bootstrap.type)
    out.all <- data.frame(rbind(c(round(c(MEMo_boot_ci$t0, sd(MEMo_boot$t), tail(c(MEMo_boot_ci[[4]]), 2)), digits = digitsForRounding), nrow(object@data.revised))))
    out.all["type"] <- "fixef.rand.coef"
    names(out.all) <- names(output)

    rwnms.rev <- unique(object@data.revised$Location)
    out.rev <- as.data.frame(round(for.all.labs(object@data.revised), digits = digitsForRounding))
    out.rev["type"] <- "revision"
    rownames(out.rev) <- rwnms.rev

    rwnms.rep <- unique(object@data.replication$Location)
    out.rep <- as.data.frame(round(for.all.labs(object@data.replication), digits = digitsForRounding))
    out.rep["type"] <- "replication"
    rownames(out.rep) <- rwnms.rep


    output <- rbind(output, out.rev, out.rep, "All" = out.all)
  } else if(object@manyLabs == "ml1"){
    # MEmo
    # Estimate & Bootstrap CI
    MEMo_boot <- boot::boot(data = object@data.replication, statistic = stat.f.MEMo, R = 1000)
    MEMo_boot_ci <- boot::boot.ci(MEMo_boot, conf = coverage_prob, type = bootstrap.type)
    out.all <- data.frame(rbind(c(round(c(MEMo_boot_ci$t0, sd(MEMo_boot$t), tail(c(MEMo_boot_ci[[4]]), 2)), digits = digitsForRounding), nrow(object@data.replication))))
    out.all["type"] <- "fixef.rand.coef"
    names(out.all) <- names(output)


    # create output for each type of study
    #out.all <- c(round(c(MEMo_boot_ci$t0, sd(MEMo_boot$t), tail(c(MEMo_boot_ci[[4]]), 2)), digits = digitsForRounding),
    #             nrow(object@data.revised),
    #             type = "fixef.rand.coef")
    rwnms.rep <- unique(object@data.replication$Location)
    out.rep <- as.data.frame(round(for.all.labs(object@data.replication), digits = digitsForRounding))
    out.rep["type"] <- "replication"
    rownames(out.rep) <- rwnms.rep

    # one ring to bind them
    output <- rbind(output, out.rev, out.rep, "All" = out.all)
  }

  # add original effect size. (I should make an external function to do this...)
  output <- output %>% tibble::add_row(data.frame(object@original), .before =  TRUE)
  # add rowname "Original" and type "original"
  rownames(output) <- c("Original", rownames(output)[-1])
  output <- output %>% dplyr::mutate(
    type = ifelse(rownames(output)=="Original", "original", type)
  )


  output[c("Rle",
          "Rls",
          "Rlp",
          "stcoef",
          "stciLow",
          "stciUp")] <- c(output[["estimate"]]/object@relevance.threshold,
                     (output[["ciLow"]])/object@relevance.threshold,
                     (output[["ciUp"]])/object@relevance.threshold,
                     output[["estimate"]],
                     output[["ciLow"]],
                     output[["ciUp"]])


  return(output)
}


# Standardise Effects of glmer-model -------------------------------------------
##' Standardise Effects of glmer-model
##' This function is used by the create_overview_table function from replication-fn.R
##' @param mod any glmer model (probabyl RandomCoefficients Model)
##' @return standardised effects + standardised CI based on
##' @author Stefan Thoma



stnd.beta.glmer <- function(mod) {
  family <- family(mod)

  #mod <- alb5@mixedModels$Models$RandomCoefficients
  b <- lme4::fixef(mod)[-1]
  se <- lme4:::summary.merMod(mod)$coefficients[-1, "Std. Error"]

  statistic <- lme4:::summary.merMod(mod)$coefficients[-1,3]


  # Prediction interval according to @borenstein2009
  df <- nrow(lme4::ranef(mod)[[1]])-2
  t <- qt(p = 0.975, df = df)
  vcov <- as.data.frame(print(lme4::VarCorr(mod), comp = "Variance"))
  var.b <- as.numeric(na.omit(vcov$vcov[vcov$var1==names(b)]))
  lo.pred <- b - t* sqrt(se^2 + var.b)
  up.pred <- b + t* sqrt(se^2 + var.b)


  if(family$family=="binomial"){
    ci <- confint(mod, parm = names(b), method = "Wald") # Profiling did not work for this example, at least not for all models tested.
  } else{
    ci <- confint(mod, parm = names(b))
  }
  sd.x <- apply(x <- as.matrix(lme4::getME(mod,"X")[,-1]),2,sd)
  sd.y <- sd(lme4::getME(mod,"y"))

  if(family$family=="binomial"){
    if(length(unique(x))>2){
      factor <- sd.x#/1.6683
    } else {factor <- 1}
  } else if(family$family == "gaussian"){
    if(length(unique(x))>2){
      factor <- sd.x/sd.y
    } else {factor <- 1/sd(resid(mod))}
  } else error("function for family", family$family, "not defined")


  data.frame("estimate" = b,"se" = se, "statistic" = statistic, "ciLow" = ci[, 1], "ciUp" = ci[, 2], "stcoef" = b*factor, "stciLow" = ci[, 1]*factor, "stciUp" = ci[, 2]*factor,
             "predUp" = up.pred, "predLow" = lo.pred, "stpredUp" = up.pred*factor, "stpredLow" = lo.pred*factor)
}

# Assigns success label to table ----------------------------------------------
##' This function takes as input the output of relevance_overview_tables
##' It creates a new column that indicates whether a replication has been successful or not.
##'
##' @param table list containing the full data and other infos from the manyLabs experiment
##' to .95
##' @param threshold relevance threshold. Needed if Rlp & Rls are not in table.
##' @returns vector indicating whether repl. has been successful or not.
##' @author Stefan Thoma
##'
##'

success <- function(table, threshold = NULL){

  if(is.null(threshold)){

    vec <- with(table, ifelse(Rls>1, "Rlv",
                              ifelse(stciLow>0 & Rlp>1, "Amb.Sig",
                                     ifelse(stciLow<0 & Rlp>1, "Amb",
                                            ifelse(stciLow>0 & Rlp<1, "Ngl.Sig",
                                                   ifelse(stciLow<0 & dplyr::between(Rlp, 0, 1), "Ngl",
                                                          ifelse(stciUp<0, "Ctr", NA)
                                                          )
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
                                ifelse(lo>0 & hi<threshold, "Ngl.Sig",
                                       ifelse(lo<0 & dplyr::between(hi, 0, threshold), "Ngl",
                                              ifelse(hi<0, "Ctr", NA)
                                       )
                                )
                         )
                  )

    )
  }

}




# Power function --------------------------------------
##' Calculates Power given a relvance table object.
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




##  plot both ------------------------------------------------------------
##' This function plots the relevance_overview_table and the difference table
##' @param object a study object
##' @param standardize should plot be based on standardized values.
##' @param coverage_probability Double. Coverage probability of the CIs. Defaults
##' @param cutoff.est is forwarded as cutoff for estimate_plot function
##' @param cutoff.diff is forwarded as cutoff for difference_plot function
##' to .95
##' @import patchwork
##' @author Stefan Thoma
##' @export
plot_both.study <- function(object, standardise = TRUE, coverage_probability = .95, cutoff.est = NULL, cutoff.diff = NULL){

  measure = object@variables$measure

 # Decide which difference table to use
  if(measure=="SMD"){
    diff.type <- "wald"
  } else if(measure == "OR"){
    diff.type <- "newcombe"
  } else{stop("function not defined for such measure not defined")}


  # adjust difference table of object:
  object@difference.table[[diff.type]] <- rbind(NA, object@difference.table[[diff.type]])


  est.ylab <- c(paste(rownames(object@table), ": Rle = ", signif(object@table$Rle, digits = digitsForRounding), "\n[", signif(object@table$Rls, digits = digitsForRounding), "; ", signif(object@table$Rlp, digits = digitsForRounding), "]", sep = ""))
  diff.ylab <- c("", paste("Diff.: ", signif(object@difference.table[[diff.type]][["stcoef"]], digits = digitsForRounding), "\n[", c(signif(object@difference.table[[diff.type]][["stciLow"]], digits = digitsForRounding)),"; ", c(signif(object@difference.table[[diff.type]][["stciUp"]], digits = digitsForRounding)), "]", sep = ""))[-2]

  rownames(object@difference.table[[diff.type]]) <- diff.ylab
  rownames(object@table) <- est.ylab


  est.plot <- plot_estimates.study(object, cutoff = cutoff.est, standardise = standardise) +
    ggplot2::theme(legend.position = "right")

  diff.plot <- plot_difference.study(object, ci.type = diff.type, cutoff = cutoff.diff, standardise = standardise) +
  # ggplot2::theme(legend.position = "right", axis.text.y = ggplot2::element_blank())
    ggplot2::theme(legend.position = "none")


# patchwork::wrap_plots(est.plot, diff.plot, guides = "collect")

  (est.plot + diff.plot ) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = 'A')


}


#' creates latex table for publication.
#'
#' The resulting table may still need to be adjusted, but serves as a sceleton.
#'
#' @param object a study object
#' @param type which type of difference table should be returned.
#'             defaults to "wald"
#' @param path where table should be saved.
#'
#' @return last called table in latex format. It does save BOTH tables (if available) to path.
#' @export
#'

create_pub_tables.study <- function(object, type = "wald", path = NULL){
  # check if slot table is empty.
  check_slot.study(object, slot = "table")

  #extract information
  est.table <- object@table

  # create file name
  file.est <- paste(path, object@name, "_est.txt", sep = "")

  if(both <- !rlang::is_empty(object@difference.table)){
    diff.table  <- object@difference.table[[type]]
    # file name for diff table:
    file.diff <- paste(path, object@name, "_diff.txt", sep = "")
  } else{print("diff table not supplied")}


  # delete me!
  table <- est.table



  # function creating each table
  prepare <- function(table, file){
    table$category <- success(table)
    if(!(identical(table$estimate, table$stcoef ))){
      vars <- c("Location", "Est.", "Std.Est.", "Rle")
    } else{
      vars <- c("Location", "Est.", "Rle")
    }
    if("p.power" %in% names(table)){ #check if table is diff table or est table
      vars <- c("type", vars, "n.total", "p.power", "r.power")
    }


    table <- table %>%
      tibble::rownames_to_column(var = "Location") %>%
      #dplyr::select(vars) %>%
      dplyr::filter(Location != "Original1") %>%
      dplyr::mutate_if(is.numeric, round,  digits = digitsForRounding) %>%
      dplyr::mutate("Est." = paste(estimate, " [", ciLow, "; ", ciUp, "]", sep = ""),
                    "Std.Est." = paste(stcoef, " [", stciLow, ";", stciUp, "]", sep = ""),
                    "Rle" = paste(Rle, " [", Rls, "; ", Rlp, "]", sep = "")) %>%
      dplyr::select(vars)


    # where should extra lines be drawn?

    extra.lines <-ifelse("type" %in% names(table),
      c(sapply(unique(table$type), function(x){min(which(table$type==x))}) -1, nrow(table)),
      c(0, nrow(table)))

      xt <- xtable::xtable(table, digits = digitsForRounding)
      xtable::print.xtable(xt, include.rownames = FALSE, booktabs = TRUE, hline.after = extra.lines, type = "latex", file = file)
      }

# apply function to estimate table
prepare(est.table, file = file.est)


if(both){
  #apply function to diff table
  prepare(diff.table, file = file.diff)
  }

}


