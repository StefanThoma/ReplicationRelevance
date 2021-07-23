###################################
###    Mixed Effect Models    ###
###################################
##' Classical Meta-Analysis procedure
##' @param input_table single element of the list effect_tables (contain effect's estimate and sd associated with it)
##' @return Overall estimated effect and CI for it. Moreover Q-value and associated p-values are reported
##' to assess heterogeneity in the estimated effects
##' @author Federico Rogai
meta_analysis <- function(input_table , level=.95, calculate_OR=FALSE) {
  meta_results <- rma(yi, vi, data=input_table, level=level, method="REML")
  effectsize <- as.vector(meta_results$beta)
  CI <- c(lower_bound = meta_results$ci.lb, upper_bound = meta_results$ci.ub)
  q_test <- c(q_statistic = round(meta_results$QE, 2), 
              p_value = round(meta_results$QEp,5))
  lis <- list(effect=effectsize, confint=CI, Q_test = q_test, ratio_vars=meta_results$se.tau2/(meta_results$se))
}


##########################################
# Function which fits two mixed models and a fixed effect model and returns objects containing them.
##' @param variables List of two elements. The first one is the dependent variable, the second one is the treatment variable
##' @param dat dataset
##' @param methodFamily which family does the error distribution belong to? ("gaussian", "binary" or "poisson")
##' @param transformation name of the transformation that should be applied to the response (e.g. log). 
##' @return list with three objects of fitted models. The first is the "empty model", i.e. a model without random effects
##' but still with the fixed treatment effect. The second is the model with the random interaction. The third is the model
##' with an interaction between treatment and (random) university effect.
##' @author Federico Rogai

MixedModels<-function(variables, dat, methodFamily, transformation=NULL){
  yVar<-variables[[1]]
  xVar<-variables[[2]]
  if(!is.null(transformation)) {
    FUN <- match.fun(transformation) 
    dat[, yVar]<-FUN(dat[,yVar])
  }
  if(tolower(methodFamily)!="gaussian"){
  fit.empty<-glm(formula=as.formula(paste0(yVar," ~", xVar)), dat, family=methodFamily, na.action = na.omit)
  fit0 <- glmer(formula=as.formula(paste0(yVar," ~", xVar, "+ (1|Location)")), dat, family=methodFamily)
  fit1 <- glmer(formula=as.formula(paste0(yVar," ~", xVar, "+ (1+", xVar ,"|Location)")), dat, family=methodFamily)
  }
  else{
    fit.empty<-aov(formula=as.formula(paste0(yVar," ~", xVar)), dat, na.action=na.omit)
    fit0 <- lmer(formula=as.formula(paste0(yVar," ~", xVar, "+ (1|Location)")), dat)
    if(xVar=="sex"){ # We had errors at convergence: this is an ad hoc solution
      fit1 <- lmer(formula=as.formula(paste0(yVar," ~", xVar, "+ (1+", xVar ,"|Location)")), dat, 
                   control = lmerControl(optimizer ="Nelder_Mead"))
    }
    else{    fit1 <- lmer(formula=as.formula(paste0(yVar," ~", xVar, "+ (1+", xVar ,"|Location)")), dat)
    
    if(isSingular(fit1)){
      fit1 <- lmer(formula=as.formula(paste0(yVar," ~", xVar, " + (1 | Location) + (" , xVar ," - 1 | Location)")), dat)
      sing <- isSingular(fit1)
    }
    }
  }
  list(fit.empty,fit0, fit1, "complexfit.singular" = sing)
}

Example<-FALSE # An example to evaluate if the functions below work
if(Example){
memo.ex<-MixedModels(list("sunkDV", "sunkgroup"),df,"gaussian")
memo.ex1<-MixedModels(list("scales", "scalesgroup"),df,"binomial")
}
##########################################
##'  Function computes a likelihood ratio test (LRT) between the model with and without random intercept. If the p-value for this test is
##'  smaller than 0.05, we compute the LRT between the model with the random slope and the one without this term (but with the 
##'  random intercept). If the larger model delivers a superior fit to the data, we print the summary of that model. Otherwise we
##'  report the summary of the model we had chosen before.
##'  Note: to compute the LRT we need to have the ML estimate. anova() automatically refits the models,
##'  whose parameters were estimated using REML
##' @param modelsToCompare Output of function MixedModels, i.e. list with three model fits which are to be compared
##' @return It returns a list with five elements.
##' The first element is the model we have chosen (model with 1-merely fixed effects, 2- random intercept only,
##' 3- with random intercept and random slope)
##' The second element reports the p-value of the hypothesis we have tested in order to choose the model
##' The third element the orginal fit for the "winning model" as from MixedModels
##' The fourth returned argument is the summary of the chosen model
##' The fifth element are the confidence interval of the winning model
##' @author Federico Rogai


LRTfunction<-function(modelsToCompare){
#modelsToCompare <- alb5$heterogen$mm  
  
  pval0<-anova(modelsToCompare[[3]], modelsToCompare[[1]], test = "LRT")[["Pr(>Chisq)"]][[2]]
  pval1<-anova(modelsToCompare[[2]], modelsToCompare[[3]], test = "LRT")[["Pr(>Chisq)"]][[2]]
  if(pval0<0.05 & pval1<0.05){
      testResult<- cat("Random interaction needed! \n Larger model fits data better than
                       purely fixed effect \n and random intercept only models. \n
                       p-values of LRT",round(pval0,digitsForRounding), " \n
                       respectively",round(pval1,digitsForRounding), ".")
      finalModel<-3
    }
  else{
    pval2<-anova(modelsToCompare[[2]], modelsToCompare[[1]], test = "LRT")[["Pr(>Chisq)"]][[2]]
    if(pval1>0.05 & pval2<0.05){
    testResult<- cat("Random intercept only is the preferred model. \n
                       p-value of LRT",round(pval1,digitsForRounding), " comparing it to the random slope model \n
                       respectively",round(pval2,digitsForRounding)," to the fixed effects only model")
    finalModel<-2
    }
    if(pval2>0.05 & pval0>0.05){
      testResult<- cat("Pure fixed effect model is the preferred model. \n
                       p-value of LRT",round(pval0,digitsForRounding), " comparing it to the random slope model \n
                       respectively",round(pval2,digitsForRounding)," to the random intercept only model")
      finalModel<-1
    }
    else{ cat("WARNING: contraddictory evidence from the LRTs")}
  }
  winningModel<-summary(modelsToCompare[[finalModel]])
  
  
  ## get intraclass correlation (icc)
  icc <- as.data.frame(icc_specs(modelsToCompare[[3]]))
  
  if(nrow(icc) == 4){
    icc <- replace_values(icc, "grp", c("Location", "Location.1", "Location.2", "Residual"), c("Intercept", "Slope", "Intercept.Slope.cov", "Residual"))
  } else if(nrow(icc) == 3){
    icc <- replace_values(icc, "grp", c("Location", "Location.1", "Residual"), c("Intercept", "Slope", "Residual"))
  }
  
  list("finalModel" = finalModel, "testResult" = testResult, modelsToCompare[finalModel], list("Summary of best model", winningModel,confint(modelsToCompare[[finalModel]], method="Wald")), "icc" = icc)
}

if(Example){
wM<-LRTfunction(memo.ex)
wM1<-LRTfunction(memo.ex1)
}
##########################################
##'  This function tests if a model in which the correlation between the two random effects (random intercept and random slope)
##'  is removed fits the data as good as one in which this term is present. If it does, it is compared to the random intercept 
##'  model only. 
##' @param ModelCorr it is the output from the function MixedModels. 
##' It is passed immediately to LRT since we need the fit of the "winning model". 
##' @return Returns an output from anova(model.1, model.2), which tells us
##' a. Is the model without correlation as good as the one in which we allow for it?
##' b. if the answer for a) is "Yes", we use a LRT as implemented in anova() to compare this model to the one
##' with the random Intercept only
##' @author Federico Rogai
TestCorrRE<-function(ModelCorr){
if(LRTfunction(ModelCorr)[[1]]!=3)  {"Model does not contain random slope"}
else{
  fitCorr<-LRTfunction(ModelCorr)[[3]][[1]] # Extracting the fitted object 
  treatment<-colnames(fitCorr@frame[2])     # Treatment name
  university<-colnames(fitCorr@frame[3])    # University (referrer)
  
  # Updating the model to make the RE uncorrelated
  fitUncorr<-update(fitCorr, formula= as.formula(paste0(".~", treatment, "+ (1 + ", treatment, "||", university,")")))
  LRTcorr<-anova(fitUncorr, fitCorr)
  if(LRTcorr[["Pr(>Chisq)"]][[2]]>0.05){
    fitREinterc<-update(fitCorr, formula= as.formula(paste0(".~", treatment, "+ (1 |", university,")")))
    return(list("The model without correlation fits the data as well as the one where we allow for it",
                "Comparing that model to the Random Intercept only model we obtain the following result",
                anova(fitREinterc, fitUncorr)))
  }
  list("The model without correlation does NOT fit the data as well as the one where we allow for it",
       LRTcorr)
}
}
if(Example){
TestCorrRE(memo.ex)
TestCorrRE(memo.ex1)
}

##########################################
##' Moderation analysis for lmer objects. We compute an F-test in anove (SS of type III) and use the p-values
##' computed by the package lmerTest to assess the significance of the interaction terms. Unfortunately such an 
##' approach cannot be used for glmer objects (lmerTest does not return p-values for this method). 
##' For the sake of consistency we use the function CI_FE_moderation to do that.
##' FUNCTION NOT USED IN THE END
##' @param ModelCorr it is the output from the function LRTfunction (we need the fit of the "winning model")
##' @return Which model is chosen? p-values supporting the decision
##' @author Federico Rogai
LmerMA<-function(fitModer){
  library(lmerTest)
  F.test<-anova(fitModer)
  indexF<-length(F.test[["Pr(>F)"]])
  if(F.test[["Pr(>F)"]][(indexF-1)]<0.05 & F.test[["Pr(>F)"]][(indexF)]<0.05 ){
    return(list("Both interaction terms are significant at the 0.05 level", F.test[["Pr(>F)"]][(indexF-1):indexF]))
  }
  if(F.test[["Pr(>F)"]][(indexF-1)]<0.05){
    return(list("The interaction term with the US- not US binary variable is significant at the 0.05 level",
                F.test[["Pr(>F)"]][(indexF-1)]))
  }
  if(F.test[["Pr(>F)"]][(indexF)]<0.05 ){
    return(list("The interaction term with the Lab- Online binary variable is significant at the 0.05 level",
                F.test[["Pr(>F)"]][(indexF)]))
  }
  else{
    list("None of the interaction terms is significant at the 0.05 level", F.test[["Pr(>F)"]])
  }
}

##########################################
##' Moderation analysis. We construct the Wlad CI for the fixed effects. If any interaction term's coefficient does 
##' cover zero, we keep the term in the model 
##' @param ModelCorr it is the output from the function LRTfunction (we need the fit of the "winning model")
##' @return Which model is chosen? p-values supporting the decision
##' @author Federico Rogai
CI_FE_moderation<-function(fitModer){
  CI<-confint(fitModer, method="Wald")
  indexCI<-nrow(CI)
  if((sign(CI[indexCI-1,1])==sign(CI[indexCI-1,2])) & (sign(CI[indexCI,1])==sign(CI[indexCI,2])) ){
    return(list("Both interaction terms are significant at the 0.05 level",CI[(indexCI-1):indexCI,1:2]))
  }
  if((sign(CI[indexCI-1,1])==sign(CI[indexCI-1,2])) ){
    return(list("The interaction term with the US- not US binary variable is significant at the 0.05 level",
                CI))
  }
  if( (sign(CI[indexCI,1])==sign(CI[indexCI,2])) ){
    return(list("The interaction term with the Lab- Online binary variable is significant at the 0.05 level",
                CI))
  }
  else{
    list("None of the interaction terms is significant at the 0.05 level",CI)
  }
}

##########################################
##'  This function inserts the fixed effects predictors "us_or_international" and "lab_or_online" in the "winning 
##'  model" from the Likelihood ratio test. These two predictors are separately made interact with the treatment
##'  effect to evaluate if they are to be considered moderators.
##' @param ModelCorr it is the output from the function LRTfunction (we need the fit of the "winning model")
##' @return See CI_FE_moderation
##' @author Federico Rogai
moderationAnalysis<-function(ModelCorr){
  bestFit<-LRTfunction(ModelCorr)[[3]][[1]] # Extracting the fitted object 
  treatment<-colnames(bestFit@frame[2])     # Treatment name
  
  fitModer<-update(bestFit, formula= as.formula(paste0(".~. +us_or_international+ lab_or_online+",
                                                       treatment, ":us_or_international +",
                                                       treatment, ":lab_or_online")))
    CI_FE_moderation(fitModer)
}

if(Example){
moderationAnalysis(memo.ex)
moderationAnalysis(memo.ex1)
}

##########################################
##'  Interaction plot for one Research Question
##' @param variables one element of the list use_var_list, containing response and explanatory variable
##' @param dat dataset used
##' @return Interaction plot
##' @author Federico Rogai
interactPlot<-function(variables, dat){ 
  y_var<-variables[[1]]
  x_var<-variables[[2]]
  set.seed(15)
  uni_sampled<-sample(unique(dat$sample),length(unique(dat$sample))/2)
  df1<-subset(dat, !is.na(get(x_var,dat)) & !is.na(get(y_var,dat))  & (dat$sample %in% uni_sampled))
  ggplot(df1, aes(x = get(x_var,df1), y = get(y_var,df1))) +
    stat_summary(fun.y = mean, geom = "line", aes(x = get(x_var,df1), y = get(y_var,df1),group = referrer)) +
    theme_bw()+
    xlab(x_var) +
    ylab(y_var)
}

if(FALSE){
interactPlot(use_vars_list[[1]], df)
}

##########################################
##'  Interaction plot for one Research Question (within colored lines for moderator)
##' @param variables one element of the list use_var_list, containing response and explanatory variable
##' @param dat dataset used
##' @return Interaction plot
##' @author Federico Rogai
moderationInteractionPlot<-function(variables, dat, moderator){ 
  y_var<-variables[[1]]
  x_var<-variables[[2]]
  mod<-get(moderator,dat)
  moderator_0<-sample(unique(dat$sample[mod==0]),length(unique(dat$sample[mod==0]))/2)
  moderator_1<-sample(unique(dat$sample[mod==1]),length(unique(dat$sample[mod==1]))/2)
  joint_sample<-c(moderator_0, moderator_1)
  df1<-subset(dat, !is.na(get(x_var,dat)) & !is.na(get(y_var,dat))  & (dat$sample %in% joint_sample))
  ggplot(df1, aes(x = get(x_var,df1), y = get(y_var,df1), color=factor(get(moderator,df1))))+
    stat_summary(fun.y = mean, geom = 'line', aes(x =  get(x_var,df1), y = get(y_var,df1), group = Location)) +
    theme_bw()+
    theme(legend.position="bottom") +
    scale_x_discrete( expand = c(0, 0))     +      
    labs(colour = paste(as.character(moderator)))+
    xlab(x_var) +
    ylab(y_var)
}


if(FALSE){
  moderationInteractionPlot(use_vars_list[[9]], df, "us_or_international")
}

##########################################
##'  Diagnostic Plot for "winning model" as defined by the likelihood ratio test. 
##' @param winningModel best model, according to LRT
##' @param percent used to avoid overplotting in Tukey-Anscombe plot. 
##' @return Tukey-Anscombe plot, Q-Q plot for RE and residuals' Q-Q plot.
##' @author Federico Rogai
diagnosticPlot<-function(winningModel, percent = 1){
  if(!{percent> 0 & percent<=1}) stop("percentage has to be between 0 and 1")
  
  model<-winningModel[[3]][[1]]
  y_name<-names(model@frame)[[1]]
  x_name<-names(model@frame)[[2]]
  show_less_points<-sample(length(resid(model)),floor(percent * length(resid(model)))) # Possible to show less points if 
                                                                      # one wants to "see more"
  # Fitted values vs residuals
  plot(resid(model)[show_less_points]~fitted(model)[show_less_points], xlab=x_name, ylab=y_name, main= "Tukey Anscombe Plot" )
  lines(loess.smooth(fitted(model)[show_less_points],resid(model)[show_less_points] , family="gaussian"), col="blue")
  # Q-Q plot of random effects
  print(qqmath(ranef(model)))
  ## Q-Q plot of residuals
  qqnorm(resid(model)[show_less_points])  
}
if(FALSE){
  diagnosticPlot(wM)
  diagnosticPlot(wM1)
}

##' Standardise Effects of lmer-model
##' This function is used by the create_overview_table function from replication-fn.R
##' @param lmer best model, according to LRT
##' @return standardised effects + standardised CI based on 
##' @author Stefan Thoma

#mod <- object@mixedModels$Models$RandomCoefficients

stnd.beta.glmer <- function(mod) {
  family <- family(mod)
  
  #mod <- alb5$heterogen$lrt[[3]][[1]]
  b <- fixef(mod)[-1]
  ci <- confint(mod, parm = names(b))
  sd.x <- apply(x <- as.matrix(getME(mod,"X")[,-1]),2,sd)
  sd.y <- sd(getME(mod,"y"))
  
  if(family$family=="binomial"){
    if(length(unique(x))>2){
    factor <- sd.x#/1.6683
    } else {factor <- 1}
  } else if(family$family == "gaussian"){
    if(length(unique(x))>2){
      factor <- sd.x/sd.y
    } else {factor <- 1/sd(resid(mod))}
  } else error("function for family", family$family, "not defined")
  
  
  data.frame("estimate" = b, "ciLow" = ci[, 1], "ciUp" = ci[, 2], "stcoef" = b*factor, "stciLow" = ci[, 1]*factor, "stciUp" = ci[, 2]*factor)
}

