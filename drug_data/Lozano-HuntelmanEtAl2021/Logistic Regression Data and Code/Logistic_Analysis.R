library(readr)
library(aod)
library(GJRM)
library(tidyr)

setwd("~/Desktop/Lab_Stuff/Ch3_Antagonisms_3-4-5/materials/DATA/Log_It_Analysis")

X5Drugs <- read_csv("SINGLES_LOGREG5.csv")
X4Drugs <- read_csv("SINGLES_LOGREG4.csv")
X3Drugs <- read_csv("SINGLES_LOGREG3.csv")

XX5Drugs <- read_csv("DOUBLES_LOGREG5.csv")
XX4Drugs <- read_csv("DOUBLES_LOGREG4.csv")
XX3Drugs <- read_csv("DOUBLES_LOGREG3.csv")

X5MOA <- read_csv("SINGLES_MOA5.csv")
X4MOA <- read_csv("SINGLES_MOA4.csv")
X3MOA <- read_csv("SINGLES_MOA3.csv")

XX5MOA <- read_csv("DOUBLES_MOA5.csv")
XX4MOA <- read_csv("DOUBLES_MOA4.csv")
XX3MOA <- read_csv("DOUBLES_MOA3.csv")

#drug logit

X5drug_logit<- glm(suppression~AMP+CPR+DOX+ERY+FOX+FUS+STR+TMP+0, data = X5Drugs, family = "binomial")
summary(X3drug_logit)


XX5drug_logit<- glm(suppression~(AMP_CPR)+(AMP_DOX)+(AMP_ERY)+(AMP_FOX)+(AMP_FUS)+(AMP_STR)+(AMP_TMP)+(CPR_DOX)+(CPR_ERY)+(CPR_FOX)+(CPR_FUS)+(CPR_STR)+(CPR_TMP)+
                      (DOX_ERY)+(DOX_FOX)+(DOX_FUS)+(DOX_STR)+(DOX_TMP)+(ERY_FOX)+(ERY_FUS)+(ERY_STR)+(ERY_TMP)+(FOX_FUS)+(FOX_STR)+
                      (FOX_TMP)+(FUS_STR)+(FUS_TMP)+(STR_TMP)+0, data = XX5Drugs, family = "binomial")
summary(XX3drug_logit)


# P-values
coef(summary(X3drug_logit))[,4]
coef(summary(XX3drug_logit))[,4]


# CIs using profiled log-likelihood
confint(X5drug_logit,level = (1-0.00625))

confint(XX5drug_logit, level =(1-0.00179))


# CIs using SE
confint.default(X3drug_logit)

exp(cbind(OR = coef(X3drug_logit), confint(X3drug_logit)))


# use for 
wald.test(Sigma = vcov(X3drug_logit), b = coef(X3drug_logit), Terms = 1:8)
vcov(X3drug_logit)

wald.test(Sigma = vcov(XX3drug_logit), b = coef(XX3drug_logit), Terms = 1:28)

#Odds-Ratio
exp(coef(X3drug_logit))
exp(coef(XX3drug_logit))

#main mech log_it
names(X5MOA)[names(X5MOA) == "30s"] <- "S30s"
names(X5MOA)[names(X5MOA) == "50s"] <- "S50s"

X5moa_logit<- glm(SUPR~C_Wall+F_Acid+DNA_gyr+S30s+S50s+0, data = X5MOA, family = "binomial")
summary(X5moa_logit)


XX5moa_logit<- glm(suppression~C_Wall_F_Acid+C_Wall_DNA_gyr+C_Wall_S30s+C_Wall_S50s+F_Acid_DNA_gyr+F_Acid_S30s+F_Acid_S50s+
              DNA_gyr_S30s+DNA_gyr_S50s+S30s_S50s+C_Wall_C_Wall+S30s_S30s+S50s_S50s+0, data = XX5MOA, family = "binomial")
summary(XX4moa_logit)


# P-values
coef(summary(X4moa_logit))[,4]
coef(summary(XX4moa_logit))[,4]


# CIs using profiled log-likelihood
confint(X5moa_logit,level = (1-0.01))

confint(XX5moa_logit,level = (1-0.0039))


# CIs using SE
confint.default(X4moa_logit)

exp(cbind(OR = coef(X4moa_logit), confint(X4moa_logit)))


# use for 
wald.test(Sigma = vcov(X4moa_logit), b = coef(X4moa_logit), Terms = 1:8)
vcov(X4moa_logit)

wald.test(Sigma = vcov(XX4moa_logit), b = coef(XX4moa_logit), Terms = 1:28)

#Odds-Ratio
exp(coef(X4moa_logit))
exp(coef(XX4moa_logit))

