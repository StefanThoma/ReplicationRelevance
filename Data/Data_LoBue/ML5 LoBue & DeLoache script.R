### Many Labs 5
### LoBue & DeLoache (2008) replication


# installing and loading required packages

install.packages("lme4") # for fitting mixed effect models
install.packages("MuMIn") # for calculationg pseudo R2
install.packages("foreign") # for easier data importing
install.packages("ggplot2") # for data visualization
install.packages("psych") # for descriptive statistics

library(lme4)
library(MuMIn)
library(foreign)
library(ggplot2)
library(psych)


# selecting the data
dataset_LoBue_raw <- read.spss(file.choose(), to.data.frame=TRUE)

# database cleaning
dataset_LoBue_cleaned <- subset(dataset_LoBue_raw, RT.correct <= 8.3 & number_errors < 8 &
                                is.na(dataset_LoBue_raw$snake_experience) == FALSE,
                                select=Site:number_errors)

# creating a new variable - site country

levels(dataset_LoBue_cleaned$Site)

dataset_LoBue_cleaned$country [dataset_LoBue_cleaned$Site == "BG"] <- "Serbia"
dataset_LoBue_cleaned$country [dataset_LoBue_cleaned$Site == "NI"] <- "Serbia"
dataset_LoBue_cleaned$country [dataset_LoBue_cleaned$Site == "NS"] <- "Serbia"
dataset_LoBue_cleaned$country [dataset_LoBue_cleaned$Site == "LY"] <- "US"

dataset_LoBue_cleaned$country <- as.factor(dataset_LoBue_cleaned$country)

levels(dataset_LoBue_cleaned$country)


# overriding the default contrast option in R (contr.treatment) with sum to zero contrast option
options(contrasts = c("contr.sum", "contr.poly")) 



# Descriptives

ggplot(dataset_LoBue_cleaned, aes(x=child_parent, y=RT.correct, fill=target_stimulus)) + 
   stat_summary(fun.y = mean, geom = "bar", position = "dodge", color="black") + 
   scale_fill_grey() + 
   xlab("Participant age") + 
   ylab("RT(s)") +
   stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2) 



describeBy(dataset_LoBue_cleaned$RT.correct, group = dataset_LoBue_cleaned$target_stimulus)

describeBy(dataset_LoBue_cleaned$RT.correct, group = dataset_LoBue_cleaned$child_parent)

describeBy(dataset_LoBue_cleaned$RT.correct, group = list(dataset_LoBue_cleaned$target_stimulus, 
                                                       dataset_LoBue_cleaned$child_parent))
sd(dataset_LoBue_cleaned$RT.correct)


# descriptive statistics by protocol

describeBy(dataset_LoBue_cleaned$RT.correct, group = dataset_LoBue_cleaned$protocol)

dataset_LoBue_cleaned_RPP <- subset(dataset_LoBue_cleaned, protocol == "NP") 
dataset_LoBue_cleaned_revised <- subset(dataset_LoBue_cleaned, protocol == "RP") 

describeBy(dataset_LoBue_cleaned_RPP$age_months, group = dataset_LoBue_cleaned_RPP$Site)
describeBy(dataset_LoBue_cleaned_RPP$age_months, group = list(dataset_LoBue_cleaned_RPP$country, 
                                                       dataset_LoBue_cleaned_RPP$child_parent))

table <- table(dataset_LoBue_cleaned_RPP$gender, dataset_LoBue_cleaned_RPP$country)
prop.table(table, 2)

levels(dataset_LoBue_cleaned$child_parent)

data_children <- subset(dataset_LoBue_cleaned_RPP, child_parent == "child")
data_parents <- subset(dataset_LoBue_cleaned_RPP, child_parent == "parent")

table <- table(data_children$gender, data_children$country)
prop.table(table, 2)

table <- table(dataset_LoBue_cleaned_RPP$gender, dataset_LoBue_cleaned_RPP$country)
prop.table(table, 2)



describeBy(dataset_LoBue_cleaned_revised$age_months, group = dataset_LoBue_cleaned_revised$Site)
describeBy(dataset_LoBue_cleaned_revised$age_months, group = list(dataset_LoBue_cleaned_revised$country, 
                                                              dataset_LoBue_cleaned_revised$child_parent))

table <- table(dataset_LoBue_cleaned_revised$gender, dataset_LoBue_cleaned_revised$country)
prop.table(table, 2)

data_children <- subset(dataset_LoBue_cleaned_revised, child_parent == "child")
data_parents <- subset(dataset_LoBue_cleaned_revised, child_parent == "parent")

table <- table(data_children$gender, data_children$country)
prop.table(table, 2)

table <- table(dataset_LoBue_cleaned_revised$gender, dataset_LoBue_cleaned_revised$country)
prop.table(table, 2)


# Confirmatory analysis - direct replication of LoBue & DeLoache (2008) and RPP (2015)
# Three main factors and their interactions as fixed effects, 
# the three-way interaction nested within site as random effects 

model_confirmatory_initial <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + 
                                  (1 + target_stimulus:child_parent:snake_experience | Site), 
                                data=dataset_LoBue_cleaned)
# model failed to converge

# removing the random effect of snake experience and its interactions

model_confirmatory2 <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + 
                              (1 + target_stimulus:child_parent | Site), 
                            data=dataset_LoBue_cleaned)
# model failed to converge


# removing the target_stimulus*child_parent interaction from the random effects
# keeping both target_stimulus and child_parent random effects

model_confirmatory3 <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + 
                              (1 + target_stimulus | Site) + 
                            (0 + child_parent | Site),
                            data=dataset_LoBue_cleaned)
# model failed to converge


# removed the random effect of child_parent

model_confirmatory3a <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + 
                                (1 + target_stimulus | Site), 
                             data=dataset_LoBue_cleaned)
# model failed to converge


# removed target_stimulus random effect
# keeping only the random intercept

model_confirmatory_final <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + 
                               (1 | Site), data=dataset_LoBue_cleaned)
# model did converge - selected as final model


summary(model_confirmatory_final)
r.squaredLR(model_confirmatory_final)

summary(model_confirmatory_final)$coefficients
summary(model_confirmatory_final)$r.squared

# https://www.r-bloggers.com/three-ways-to-get-parameter-specific-p-values-from-lmer/
# extract coefficients
coefs <- data.frame(coef(summary(model_confirmatory_final)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs


# final model without main effect of target_stimulus and all its interactions

model_comparison_A <- lmer(RT.correct ~ child_parent*snake_experience + 
                                (1 | Site), data=dataset_LoBue_cleaned)

summary(model_comparison_A)
r.squaredLR(model_comparison_A)

anova(model_comparison_A, model_confirmatory_final)
r.squaredLR(model_confirmatory_final, model_comparison_A)



# a new baseline model which only includes main target_stimulus effect
# and child_parent:snake_experience interaction

model_comparison_B <- lmer(RT.correct ~ target_stimulus + child_parent*snake_experience + 
                              (1 | Site), data=dataset_LoBue_cleaned)

summary(model_comparison_B)
r.squaredLR(model_comparison_B)

coefs <- data.frame(coef(summary(model_comparison_B)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

anova(model_comparison_B, model_comparison_A)
r.squaredLR(model_comparison_B, model_comparison_A)


# comparison of model including fixed effects of snake_experience, child_parent and target_stimulus
# as well as snake_experience*child_parent interaction with a model that also includes
# target_stimulus*snake_experiece, target_stimulus*child_parent and the three-way interaction

anova(model_comparison_B, model_confirmatory_final)
r.squaredLR(model_confirmatory_final, model_comparison_B)





# Testing for moderation by protocol version
# Creating a model with protocol added as a fixed factor 
# (along with all interactions with the remaining three factors)

model_protocol <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience*protocol + 
                          (1 | Site), data=dataset_LoBue_cleaned)
summary(model_protocol)
r.squaredLR(model_protocol)

# Comparing the model including protocol with the baseline confirmatory model

anova(model_confirmatory_final, model_protocol)
r.squaredLR(model_protocol, model_confirmatory_final)



# comparing the model including protocol without interactions with baseline

model_comparison_C <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + protocol + 
                         (1 | Site), data=dataset_LoBue_cleaned)
summary(model_comparison_C)

anova(model_confirmatory_final, model_comparison_C)
r.squaredLR(model_comparison_C, model_confirmatory_final)


# model including protocol and target_stimulus:protocol interaction

model_comparison_D <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + protocol + 
                           target_stimulus:protocol +
                           (1 | Site), data=dataset_LoBue_cleaned)
summary(model_comparison_D)

anova(model_comparison_C, model_comparison_D)
r.squaredLR(model_comparison_D, model_comparison_C)


# model including protocol and child_parent:protocol interaction

model_comparison_E <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + protocol + 
                           child_parent:protocol +
                           (1 | Site), data=dataset_LoBue_cleaned)

summary(model_comparison_E)

anova(model_comparison_C, model_comparison_E)
r.squaredLR(model_comparison_E, model_comparison_C)



# Analyses by protocol


describeBy(dataset_LoBue_cleaned$RT.correct, group = dataset_LoBue_cleaned$protocol)

dataset_LoBue_cleaned_RPP <- subset(dataset_LoBue_cleaned, protocol == "NP") 
dataset_LoBue_cleaned_revised <- subset(dataset_LoBue_cleaned, protocol == "RP") 


describeBy(dataset_LoBue_cleaned_RPP$RT.correct, group = dataset_LoBue_cleaned_RPP$target_stimulus)
sd(dataset_LoBue_cleaned_RPP$RT.correct)

describeBy(dataset_LoBue_cleaned_revised$RT.correct, group = dataset_LoBue_cleaned_revised$target_stimulus)
sd(dataset_LoBue_cleaned_revised$RT.correct)


model_confirmatory_final_RPP <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + 
                              (1 | Site), data=dataset_LoBue_cleaned_RPP)
summary(model_confirmatory_final_RPP)
coefs <- data.frame(coef(summary(model_confirmatory_final_RPP)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

model_comparison_A_RPP <- lmer(RT.correct ~ child_parent*snake_experience + 
                                (1 | Site), data=dataset_LoBue_cleaned_RPP)

model_comparison_B_RPP <- lmer(RT.correct ~ target_stimulus + child_parent*snake_experience + 
                                (1 | Site), data=dataset_LoBue_cleaned_RPP)
summary(model_comparison_B_RPP)
coefs <- data.frame(coef(summary(model_comparison_B_RPP)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

# model with target stimulus and all interactions compared to child_parent*snake experience
anova(model_comparison_A_RPP, model_confirmatory_final_RPP)
r.squaredLR(model_confirmatory_final_RPP, model_comparison_A_RPP)

# model with target stimulus but no interactions compared to child_parent*snake experience
anova(model_comparison_A_RPP, model_comparison_B_RPP)
r.squaredLR(model_comparison_B_RPP, model_comparison_A_RPP)




model_confirmatory_final_revised <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + 
                                  (1 | Site), data=dataset_LoBue_cleaned_revised)
summary(model_confirmatory_final_revised)
coefs <- data.frame(coef(summary(model_confirmatory_final_revised)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

model_comparison_A_revised <- lmer(RT.correct ~ child_parent*snake_experience + 
                                    (1 | Site), data=dataset_LoBue_cleaned_revised)

model_comparison_B_revised <- lmer(RT.correct ~ target_stimulus + child_parent*snake_experience + 
                                    (1 | Site), data=dataset_LoBue_cleaned_revised)
summary(model_comparison_B_revised)
coefs <- data.frame(coef(summary(model_comparison_B_revised)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs


# model with target stimulus and all interactions compared to child_parent*snake experience
anova(model_comparison_A_revised, model_confirmatory_final_revised)
r.squaredLR(model_confirmatory_final_revised, model_comparison_A_revised)

# model with target stimulus but no interactions compared to child_parent*snake experience
anova(model_comparison_A_revised, model_comparison_B_revised)
r.squaredLR(model_comparison_B_revised, model_comparison_A_revised)



### exploratory analyses by country

# creating a new variable - site country

levels(dataset_LoBue_cleaned$Site)

dataset_LoBue_cleaned$country [dataset_LoBue_cleaned$Site == "BG"] <- "Serbia"
dataset_LoBue_cleaned$country [dataset_LoBue_cleaned$Site == "NI"] <- "Serbia"
dataset_LoBue_cleaned$country [dataset_LoBue_cleaned$Site == "NS"] <- "Serbia"
dataset_LoBue_cleaned$country [dataset_LoBue_cleaned$Site == "LY"] <- "US"

dataset_LoBue_cleaned$country <- as.factor(dataset_LoBue_cleaned$country)

is.factor(dataset_LoBue_cleaned$country)
levels(dataset_LoBue_cleaned$country)


model_country_main_effect <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + 
                                    country + 
                                   (1 | Site), data=dataset_LoBue_cleaned)

model_country_interaction <- lmer(RT.correct ~ target_stimulus*child_parent*snake_experience + 
                                    country + target_stimulus:country +
                                    (1 | Site), data=dataset_LoBue_cleaned)


# comparing model with country:target_stimulus interaction to model with country only
anova(model_country_main_effect, model_country_interaction)
r.squaredLR(model_country_interaction, model_country_main_effect)

describeBy(dataset_LoBue_cleaned$RT.correct, group = list(dataset_LoBue_cleaned$target_stimulus, 
                                                          dataset_LoBue_cleaned$country))



### analyses without Novi Sad
### all analyses use the same code as for the full sample, 
### only participants from Novi Sad are excluded

levels(dataset_LoBue_cleaned$Site)
dataset_LoBue <- subset(dataset_LoBue_cleaned, Site != "NS",
                                select=Site:number_errors)

levels(dataset_LoBue$Site)


