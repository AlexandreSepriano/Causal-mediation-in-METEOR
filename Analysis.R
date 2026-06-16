#<<<<<<<<<<<########################################################>>>>>>>>>>>#
#<<<<<<<<<<<########## CAUSAL MEDIATION ANALYSIS ######################>>>>>>>>#
#<<<<<<<<<<<########################################################>>>>>>>>>>>#


#########################################################################
####################### Libraries  ######################################
#########################################################################


library(mediation)
library(haven)
library(ltmle)
library(dplyr)


#########################################################################
####################### Load the data  ##################################
#########################################################################


setwd("your path")
causalmediation <- read_dta("causalmediation.dta")


#########################################################################
####################### Time-fixed mediation gformula ###################
#######################     Supplementary Box S1      ###################
#########################################################################


model.m1 <- lm(asdastotalt2 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew
               , data=causalmediation)

model.y1 <- lm(basfitotalt2 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew + asdastotalt2
               , data=causalmediation)

set.seed(1234)
out.1 <- mediate(model.m1, model.y1, sims = 1000,  treat = "bionew", mediator = "asdastotalt2", 
                 control.value = 0, treat.value = 1, INT = FALSE)
summary(out.1)
plot(out.1)


## Standard errors
se_d0  <- sd(out.1$d0.sims)
se_z0  <- sd(out.1$z0.sims)
se_tau <- sd(out.1$tau.sims)
se_n0  <- sd(out.1$n0.sims)

effects <- data.frame(
  Effect   = c("ACME", "ADE", "Total Effect", "Proportion Mediated"),
  Estimate = c(out.1$d0, out.1$z0, out.1$tau.coef, out.1$n0),
  SE       = c(se_d0, se_z0, se_tau, se_n0),
  CI.lower = c(out.1$d0.ci[1], out.1$z0.ci[1], out.1$tau.ci[1], out.1$n0.ci[1]),
  CI.upper = c(out.1$d0.ci[2], out.1$z0.ci[2], out.1$tau.ci[2], out.1$n0.ci[2]))

effects



#########################################################################
####################### Time-fixed mediation gformula ###################
########## Ignoring  post-treatment mediator-outcome confounders ########
#######################     Supplementary Box S2      ###################
#########################################################################


model.m2 <- lm(asdastotalt2 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew
               , data=causalmediation)

model.y2 <- lm(basfitotalt2 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew*asdastotalt2
               , data=causalmediation)

set.seed(1234)
out.2 <- mediate(model.m2, model.y2, sims = 1000,  treat = "bionew", mediator = "asdastotalt2", 
                 control.value = 0, treat.value = 1, INT = TRUE)
summary(out.2)
plot(out.2)


test.TMint(out.2, conf.level = .95)

## Standard errors

# ACME
se_d0      <- sd(out.2$d0.sims)
se_d1      <- sd(out.2$d1.sims)
se_d_avg   <- sd(out.2$d.avg.sims)
# ADE
se_z0      <- sd(out.2$z0.sims)
se_z1      <- sd(out.2$z1.sims)
se_z_avg   <- sd(out.2$z.avg.sims)
# Total effect
se_tau     <- sd(out.2$tau.sims)
# Proportion mediated
se_n0      <- sd(out.2$n0.sims)
se_n1      <- sd(out.2$n1.sims)
se_n_avg   <- sd(out.2$n.avg.sims)
effects <- data.frame(
  Effect = c(
    "ACME (control)", "ACME (treated)", "ACME (average)",
    "ADE (control)",  "ADE (treated)",  "ADE (average)",
    "Total Effect",
    "Prop. Mediated (control)", "Prop. Mediated (treated)", "Prop. Mediated (average)" ),
  Estimate = c(
    out.2$d0, out.2$d1, out.2$d.avg,
    out.2$z0, out.2$z1, out.2$z.avg,
    out.2$tau.coef,
    out.2$n0, out.2$n1, out.2$n.avg ),
  
  SE = c(
    se_d0, se_d1, se_d_avg,
    se_z0, se_z1, se_z_avg,
    se_tau,
    se_n0, se_n1, se_n_avg),
  
  CI.lower = c(
    out.2$d0.ci[1], out.2$d1.ci[1], out.2$d.avg.ci[1],
    out.2$z0.ci[1], out.2$z1.ci[1], out.2$z.avg.ci[1],
    out.2$tau.ci[1],
    out.2$n0.ci[1], out.2$n1.ci[1], out.2$n.avg.ci[1]),
  
  CI.upper = c(
    out.2$d0.ci[2], out.2$d1.ci[2], out.2$d.avg.ci[2],
    out.2$z0.ci[2], out.2$z1.ci[2], out.2$z.avg.ci[2],
    out.2$tau.ci[2],
    out.2$n0.ci[2], out.2$n1.ci[2], out.2$n.avg.ci[2])
)

effects



#########################################################################
####################### TMLE for the Total Treatment Effect #############
#######################     Supplementary Box S5      ###################
#########################################################################

#################################### with post-treatment confounders
#>>>>>>>>>>># This version does not correspond to the algorithm in box 5 

tmledata <- causalmediation
tmledata <- causalmediation %>% select(-id)
tmledata <- tmledata[, c("mny", "asasmri", "hla", "age", "sex", "comorbbin", "pertvt1", 
                         "ibdbl", "emmtvt1", "comedtvt1", "asdastotalt1", "basfitotalt1", 
                         "bionew", "pertvt2", "emmtvt2", "comedtvt2", "asdastotalt2", 
                         "basfitotalt2")]


notreat <- function(row) c(0)
alltreat <- function(row) c(1) 

Lnodes <- c("pertvt2","emmtvt2","comedtvt2","asdastotalt2")
Anodes <- c("bionew")
Ynodes <- c("basfitotalt2")
Qform <- c("pertvt2"      = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew",
           "emmtvt2"      = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew",
           "comedtvt2"    = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew",
           "asdastotalt2" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew + pertvt2 + emmtvt2 + comedtvt2",
           "basfitotalt2" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew + asdastotalt2")
gform <- c("bionew" = "bionew ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1")
EY.11 <- ltmle(tmledata, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes, 
               Qform = Qform, gform = gform, 
               rule=alltreat, estimate.time = FALSE)
EY.00 <- ltmle(tmledata, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes, 
               Qform = Qform, gform = gform, 
               rule=notreat, estimate.time = FALSE)
ATE <- ltmle(tmledata, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes,
             Qform = Qform, gform = gform,
             rule = list(alltreat, notreat), estimate.time = FALSE)
print(summary(ATE))


#################################### without post-treatment confounders but with ASDAS at 6M (mediator)
#>>>>>>>>>>># This version does not correspond to the algorithm in box 5 

tmledata_nomoc<- tmledata %>% select(-pertvt2,-emmtvt2, -comedtvt2)

notreat <- function(row) c(0)
alltreat <- function(row) c(1) 

Lnodes <- c("asdastotalt2")
Anodes <- c("bionew")
Ynodes <- c("basfitotalt2")
Qform <- c("asdastotalt2" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew",
           "basfitotalt2" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew + asdastotalt2")
gform <- c("bionew" = "bionew ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1")
EY.11 <- ltmle(tmledata_nomoc, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes, 
               Qform = Qform, gform = gform, 
               rule=alltreat, estimate.time = FALSE)
EY.00 <- ltmle(tmledata_nomoc, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes, 
               Qform = Qform, gform = gform, 
               rule=notreat, estimate.time = FALSE)
ATE <- ltmle(tmledata_nomoc, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes,
             Qform = Qform, gform = gform,
             rule = list(alltreat, notreat), estimate.time = FALSE)
print(summary(ATE))



#################################### without post-treatment confounders and without ASDAS at 6M (mediator)
#>>>>>>>>>>># This version corresponds to the algorithm in box 5 

tmledata_nomoc_nomed<- tmledata_nomoc %>% select(-asdastotalt2)

notreat <- function(row) c(0)
alltreat <- function(row) c(1) 

Lnodes <- NULL
Anodes <- c("bionew")
Ynodes <- c("basfitotalt2")
Qform <- c("basfitotalt2" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1 + bionew")
gform <- c("bionew" = "bionew ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + asdastotalt1 + basfitotalt1")
EY.11 <- ltmle(tmledata_nomoc_nomed, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes, 
               Qform = Qform, gform = gform, 
               rule=alltreat, estimate.time = FALSE)
EY.00 <- ltmle(tmledata_nomoc_nomed, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes, 
               Qform = Qform, gform = gform, 
               rule=notreat, estimate.time = FALSE)
ATE <- ltmle(tmledata_nomoc_nomed, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes,
             Qform = Qform, gform = gform,
             rule = list(alltreat, notreat), estimate.time = FALSE)
print(summary(ATE))


##### Compute E-values
beta <- summary(ATE)$effect.measures$ATE$estimate
sd_y <- sd(tmledata_nomoc_nomed$basfitotalt2)

rr <- exp(0.91 * abs(beta / sd_y))
evalue_point <- rr + sqrt(rr * (rr - 1))
evalue_point


##### Compute diagnostics

g_trunc_lb <- 0.025
g_trunc_ub <- 0.975

ps <- ATE$cum.g[, dim(ATE$cum.g)[2], 1] 
pct_gtrunc <- 100 * mean(ps < g_trunc_lb | ps > g_trunc_ub)

IC <- ATE$IC
if (is.matrix(IC) && ncol(IC) == 2) IC <- IC[,1] - IC[,2]
if (is.list(IC)) IC <- IC$IC1 - IC$IC2
sd_IC <- sd(IC)

cat("% g-truncated:", round(pct_gtrunc, 1), "%\n")
cat("SD of influence curve:", round(sd_IC, 3), "\n")






