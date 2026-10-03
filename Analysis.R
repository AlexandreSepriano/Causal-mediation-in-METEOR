#<<<<<<<<<<<########################################################>>>>>>>>>>>#
#<<<<<<<<<<<########## CAUSAL MEDIATION ANALYSIS IN METEOR #########>>>>>>>>>>>#
#<<<<<<<<<<<##########   bDMARDs -> ASDAS -> BASFI at 6 months  ####>>>>>>>>>>>#
#<<<<<<<<<<<########################################################>>>>>>>>>>>#

#===============================================================================
# This script is the R version of Stata's Analysis.do. Every estimator in the
# supplementary tables is computed twice: once by hand, following the algorithm
# written in the corresponding Supplementary Box, and once with the R
# package that implements it.
#
# PROGRAMS AND THE SUPPLEMENTARY TABLES THEY CORRESPOND TO 
#
#   Table  Content                                 Box    Stata                     R
#   -----  --------------------------------------  -----  ------------------------  -----------------------
#   S2     Total effect                            S1     manual, medeff, gformula  manual, mediation
#   S3     Mediation ignoring MOC                  S2     manual, medeff, gformula  manual, mediation
#   S4     Two mediators, main paths               S2.1   manual x2                 manual x2
#   S5     Two mediators, all paths                S2.1   manual x2                 manual x2, paths
#   S6     Total effect, MSM with IPTW             S3     manual, teffects ipw      manual, ipw + survey
#   S7     Total effect, TMLE                      S4     manual, teffects aipw     manual, ltmle
#   S8     Balance, original vs pseudo-population  S3     msmbalanceV1 (run 5)      manual (section 5)
#   S9     Mediation within strata of baseline     S2     manual (run 7)            manual (section 7)
#   S9.1   Two mediators within strata of baseline S2.1   manual (run 7)            manual (section 7)
#   S10    Mediation accounting for MOC            S6     manual, gformula          manual
#
#   S4 and S5 each hold two specifications, with and without the two significant
#   interactions in the outcome model, so one table shows both side by side:
#     S4  main paths (Direct, PSE_M1, PSE_M2, Total), without and with interaction
#     S5  all paths (the four above plus the two splits of PSE_M1 and the three
#         mediator-scale effects), without and with interaction
#
#   Notes on R
#     - There is no R package for the stand-alone time-fixed g-formula, so the
#       total effect from the mediation package is used in S2.
#     - Stata's gformula does mediation with post-treatment confounders; the R
#       gfoRmula package does not, so S10 has no R package row.
#     - The paths package is run without interaction only
#
# Results are appended in estimates_r.csv and diagnostics_r.csv in the Tables
# folder and combined into estimates_r.xlsx, one sheet per supplementary table,
# in the same way as estimates_stata.xlsx.
#===============================================================================


#########################################################################
####################### Libraries  ######################################
#########################################################################

library(writexl)   # export to the Excel file
library(mediation) # medeff's counterpart (Boxes S1 and S2)
library(paths)     # edge g-formula for ordered mediators (Box S2.1)
library(ltmle)     # TMLE (Box S4)
library(ipw)       # inverse probability of treatment weights (Box S3)
library(survey)    # marginal structural model with a robust SE (Box S3)
library(ggplot2)   # Figures

#########################################################################
####################### Settings  #######################################
#########################################################################

SEED  <- 1234   # one seed for every bootstrap, as in Analysis.do
NBOOT <- 1000   # bootstrap replications
DEC   <- 2      # decimals for the estimates in the excel file
DECD  <- 3      # decimals for the diagnostics and p-values
PCUT  <- 0.20   # significance level for keeping an interaction

PATH   <- "C:/Users/alexa/OneDrive/Work/Projects/Causal_axSpA/METEOR/Data/Mediation manuscript/"
DATA   <- paste0(PATH, "Datasets/")
TABLES <- paste0(PATH, "Tables/")
FIGURES <- paste0(PATH, "Figures/")


#########################################################################
####################### Prepare the dataset  ############################
#########################################################################

#===============================================================================
# Dataset:   originalfulllong800.csv, the same file Analysis.do reads.
#
# Cohort:    population2 == 1 & allcompleters == 1, then the 6-month row only.
#            n = 419 over two visits: complete ASDAS and BASFI at both visits
#            and complete baseline covariates. High disease activity at baseline
#            is already applied at the database.
#
# Treatment: bionew. On the 6-month row it records the bDMARD use between
#            baseline to 6 months.
#===============================================================================

d <- read.csv(paste0(DATA, "originalfulllong800.csv"), header = TRUE, sep = ",")

d$population3 <- d$population2          # same role as population3, over two visits
d <- d[d$population3 == 1 & d$allcompleters == 1 & d$t == 6, ]
stopifnot(nrow(d) == 419)
cat("patients in the analysis:", nrow(d), "\n")

# ASDAS without CRP, used by the two-mediator analyses
d$asdastotalt1_pro <- d$asdastotalt1 - 0.579 * log(d$crpt1 + 1)
d$asdastotalt2_pro <- d$asdastotalt2 - 0.579 * log(d$crpt2 + 1)

#####>>>>>> Variables (the globals of Stata's Analysis.do)

W    <- c("age", "sex", "comorbbin", "mny", "asasmri", "hla",
          "pertvt1", "ibdbl", "emmtvt1", "comedtvt1", "asdastotalt1", "basfitotalt1")
A    <- "bionew"
M    <- "asdastotalt2"
MOC  <- c("pertvt2", "emmtvt2", "comedtvt2")
Y    <- "basfitotalt2"

# Two-mediator setting: CRP -> ASDAS-PRO -> BASFI. asdastotalt1 is replaced by
# the baseline values of the two mediators, so it is dropped from W here.
W2   <- c("age", "sex", "comorbbin", "mny", "asasmri", "hla",
          "pertvt1", "ibdbl", "emmtvt1", "comedtvt1", "basfitotalt1")
M1   <- "crpt2"
M2   <- "asdastotalt2_pro"
bM1  <- "crpt1"
bM2  <- "asdastotalt1_pro"

ANA_N <- nrow(d)
SD_Y  <- sd(d[[Y]])

stopifnot(!anyNA(d[, c(W, A, M, MOC, Y, M1, M2, bM1, bM2)]))


#########################################################################
####################### Helpers  ########################################
#########################################################################

#####>>>>>> Formula built from character vectors: the model specifications

f <- function(lhs, rhs) as.formula(paste(lhs, "~", paste(rhs, collapse = " + ")))

#####>>>>>> Non-parametric percentile bootstrap over patients

boot_pct <- function(est, data, R = NBOOT, seed = SEED) {
  set.seed(seed)
  point <- est(data)
  reps  <- matrix(NA_real_, nrow = R, ncol = length(point),
                  dimnames = list(NULL, names(point)))
  n <- nrow(data)
  for (b in seq_len(R)) {
    rs <- data[sample.int(n, n, replace = TRUE), ]
    v  <- try(est(rs), silent = TRUE)
    if (!inherits(v, "try-error")) reps[b, ] <- v
  }
  data.frame(
    effect   = names(point),
    estimate = as.numeric(point),
    se       = apply(reps, 2, sd, na.rm = TRUE),
    lci      = apply(reps, 2, quantile, 0.025, na.rm = TRUE),
    uci      = apply(reps, 2, quantile, 0.975, na.rm = TRUE),
    reps     = apply(reps, 2, function(x) sum(!is.na(x))),
    row.names = NULL, stringsAsFactors = FALSE)
}

#####>>>>>> E-value for a continuous outcome (Supplementary Box S5)

evalue_md <- function(beta, sd_y = SD_Y) {
  rr <- exp(0.91 * abs(beta / sd_y))
  rr + sqrt(rr * (rr - 1))
}

#####>>>>>> HC1 robust standard error, matching Stata's "robust"

robust_se <- function(m, w = NULL) {
  X <- model.matrix(m)
  u <- residuals(m)
  if (is.null(w)) w <- rep(1, nrow(X))
  bread <- solve(crossprod(X * w, X))
  meat  <- crossprod(X * w * u)
  n <- nrow(X); k <- ncol(X)
  sqrt(diag(bread %*% meat %*% bread * n / (n - k)))
}


#########################################################################
####################### Collect the estimates  ##########################
#########################################################################

#===============================================================================
# The counterpart of Programs/saveestV1.do. Each estimate is appended to
# estimates_r.csv and each diagnostic to diagnostics_r.csv; save_sheets() then
# rebuilds estimates_r.xlsx with one sheet per supplementary table.
#===============================================================================

EST_FILE  <- paste0(TABLES, "estimates_r.csv")
DIAG_FILE <- paste0(TABLES, "diagnostics_r.csv")
XLS_FILE  <- paste0(TABLES, "estimates_r.xlsx")

savest <- function(tbl, est, eff, b, se = NA, lo = NA, hi = NA,
                   ev = NA, reps = NBOOT, n = ANA_N) {
  # Some estimators report an interval but no standard error. It is backed out
  # of the interval so the SE column is never empty: se = (uci - lci)/(2*1.96)
  if (is.na(se) && !is.na(lo) && !is.na(hi)) se <- (hi - lo) / (2 * 1.96)
  row <- data.frame(table = tbl, estimator = est, effect = eff,
                    estimate = b, se = se, lci = lo, uci = hi,
                    evalue = ev, reps = reps, n = n, stringsAsFactors = FALSE)
  write.table(row, EST_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(EST_FILE), append = file.exists(EST_FILE))
  invisible(row)
}

#####>>>>>> Write a whole boot_pct() result in one call

savebs <- function(tbl, est, bt, ev = NA) {
  for (i in seq_len(nrow(bt)))
    savest(tbl, est, bt$effect[i], bt$estimate[i], bt$se[i], bt$lci[i], bt$uci[i],
           if (i == 1) ev else NA, bt$reps[i])
  invisible(bt)
}

savediag <- function(tbl, est, diag, val) {
  row <- data.frame(table = tbl, estimator = est, diagnostic = diag,
                    value = val, stringsAsFactors = FALSE)
  write.table(row, DIAG_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(DIAG_FILE), append = file.exists(DIAG_FILE))
  invisible(row)
}

#####>>>>>> Rebuild the excel file, one sheet per supplementary table.
#####>>>>>> The last rows written for an estimator, so re-running an
#####>>>>>> estimator replaces its rows instead of adding to them.
#####>>>>>>
#####>>>>>> write_xlsx cannot replace a single sheet: it rewrites the whole
#####>>>>>> file. So this always rebuilds every table that has rows in the two
#####>>>>>> csv files, which is what makes it safe to call after any section.
#####>>>>>> The csv files are the record; the excel file is only a view of them.

save_sheets <- function() {
  tables <- c("S2","S3","S4","S5","S6","S7","S10","T2")
  fmt  <- function(x, k) ifelse(is.na(x), "", sprintf(paste0("%.", k, "f"), x))
  last <- function(df, keys) df[!duplicated(df[keys], fromLast = TRUE), , drop = FALSE]

  e <- if (file.exists(EST_FILE))  read.csv(EST_FILE,  stringsAsFactors = FALSE) else NULL
  g <- if (file.exists(DIAG_FILE)) read.csv(DIAG_FILE, stringsAsFactors = FALSE) else NULL
  if (!is.null(e)) e <- last(e, c("table","estimator","effect"))
  if (!is.null(g)) g <- last(g, c("table","estimator","diagnostic"))

  sheets <- list()
  for (tb in tables) {
    blocks <- list()
    if (!is.null(e) && any(e$table == tb)) {
      x <- e[e$table == tb, ]
      ci <- fmt(x$estimate, DEC)
      has <- !is.na(x$lci) & !is.na(x$uci)
      ci[has] <- sprintf("%s (%s; %s)", ci[has], fmt(x$lci[has], DEC), fmt(x$uci[has], DEC))
      blocks$estimates <- data.frame(
        estimator = x$estimator, effect = x$effect, est_95CI = ci,
        se = fmt(x$se, DEC), estimate = fmt(x$estimate, DEC),
        lci = fmt(x$lci, DEC), uci = fmt(x$uci, DEC),
        evalue = fmt(x$evalue, DEC), reps = x$reps, n = x$n,
        stringsAsFactors = FALSE)
    }
    if (!is.null(g) && any(g$table == tb)) {
      x <- g[g$table == tb, ]
      blocks$diagnostics <- data.frame(
        estimator = x$estimator, diagnostic = x$diagnostic,
        value = fmt(x$value, DECD), stringsAsFactors = FALSE)
    }
    if (!length(blocks)) next
    # the two blocks share one sheet, separated by a blank row
    out <- blocks$estimates
    if (!is.null(blocks$diagnostics)) {
      pad <- function(df, k) { while (ncol(df) < k) df[[paste0("v", ncol(df) + 1)]] <- ""; df }
      k   <- max(ncol(blocks$estimates), ncol(blocks$diagnostics), 3)
      nm  <- if (is.null(out)) names(pad(blocks$diagnostics, k)) else names(pad(out, k))
      dia <- pad(blocks$diagnostics, k); names(dia) <- nm
      hdr <- as.data.frame(as.list(c("estimator", "diagnostic", "value",
                                     rep("", k - 3))), stringsAsFactors = FALSE)
      names(hdr) <- nm
      blk <- as.data.frame(as.list(rep("", k)), stringsAsFactors = FALSE); names(blk) <- nm
      out <- if (is.null(out)) rbind(hdr, dia) else
             rbind(pad(out, k), blk, hdr, dia)
    }
    sheets[[tb]] <- out
    cat(sprintf("sheet %s: %d estimates, %d diagnostics\n", tb,
                if (is.null(blocks$estimates)) 0 else nrow(blocks$estimates),
                if (is.null(blocks$diagnostics)) 0 else nrow(blocks$diagnostics)))
  }
  if (length(sheets)) write_xlsx(sheets, XLS_FILE)
  cat("workbook written:", XLS_FILE, "\n")
  invisible(sheets)
}


#########################################################################
############## 1 - SUPPLEMENTARY BOX S1: total effect  ##################
##############        time-fixed g-formula (manual)    ##################
##############                          -> Table S2    ##################
#########################################################################

#===============================================================================
# Step 1  outcome model and mediator model on the observed data
# Step 2  simulate the mediator under both regimens, then the outcome
# Step 3  the total effect is the difference of the two means
# Step 4  percentile bootstrap
#===============================================================================

setA <- function(z, a) { z[[A]] <- a; z }

gf_total <- function(z) {

  ## Step 1
  my <- lm(f(Y, c(M, A, W)), z)   # outcome model
  mm <- lm(f(M, c(A, W)),    z)   # mediator model
  at <- glm(f(A, W), binomial, z) # treatment model, for the diagnostics only

  ## Step 2.1  counterfactual ASDAS at 6 months
  z1 <- setA(z, 1); z1[[M]] <- predict(mm, z1)
  z0 <- setA(z, 0); z0[[M]] <- predict(mm, z0)

  ## Step 2.2  counterfactual BASFI at 6 months
  m11 <- mean(predict(my, z1))
  m00 <- mean(predict(my, z0))

  ## Natural course, for the estimated mean of the observed regimen
  zp <- z
  zp[[A]] <- rbinom(nrow(z), 1, predict(at, z, type = "response"))
  zp[[M]] <- predict(mm, zp)

  ## Step 3
  c(ATE       = m11 - m00,
    Mean_Y1M1 = m11,
    Mean_Y0M0 = m00,
    Mean_A    = mean(z[[A]]),
    Mean_Ap   = mean(zp[[A]]),
    Mean_M    = mean(z[[M]]),
    Mean_Mp   = mean(zp[[M]]),
    Mean_Y    = mean(z[[Y]]),
    Mean_Yp   = mean(predict(my, zp)))
}

bt1 <- boot_pct(gf_total, d)
print(bt1, digits = 3)

ev1 <- evalue_md(bt1$estimate[bt1$effect == "ATE"])
cat("E-value (point estimate):", round(ev1, 3), "\n")

savebs("S2", "G-formula (manual)", bt1, ev1)


#########################################################################
############## 1.1 - SUPPLEMENTARY BOX S1: total effect  ################
##############               (R mediation package)      ################
##############                           -> Table S2    ################
#########################################################################

model.m1 <- lm(f(M, c(A, W)),    data = d)
model.y1 <- lm(f(Y, c(M, A, W)), data = d)

set.seed(SEED)
out.1 <- mediate(model.m1, model.y1, sims = NBOOT, treat = A, mediator = M,
                 control.value = 0, treat.value = 1, INT = FALSE)
summary(out.1)

savest("S2", "mediation", "ACME",  out.1$d0,       sd(out.1$d0.sims),
       out.1$d0.ci[1],  out.1$d0.ci[2])
savest("S2", "mediation", "ADE",   out.1$z0,       sd(out.1$z0.sims),
       out.1$z0.ci[1],  out.1$z0.ci[2])
savest("S2", "mediation", "Total", out.1$tau.coef, sd(out.1$tau.sims),
       out.1$tau.ci[1], out.1$tau.ci[2],
       evalue_md(out.1$tau.coef))
savest("S2", "mediation", "PM",    out.1$n0,       sd(out.1$n0.sims),
       out.1$n0.ci[1],  out.1$n0.ci[2])

save_sheets()


#########################################################################
############## 2 - SUPPLEMENTARY BOX S2: mediation      #################
##############      ignoring post-treatment MOC (manual) ################
##############                          -> Table S3     #################
#########################################################################

#####>>>>>> Step 0  test for a treatment-mediator interaction.
#####>>>>>> The product term is fitted once on the observed data; its p-value
#####>>>>>> decides whether the interaction is carried into the outcome models.

int_fit <- lm(f(Y, c(M, A, W, paste0(M, ":", A))), d)
p_inter <- summary(int_fit)$coefficients[paste0(M, ":", A), "Pr(>|t|)"]
INTER   <- p_inter < PCUT

cat(sprintf("\nTreatment-mediator interaction (%s x %s): p = %.4f\n", A, M, p_inter))
cat(if (INTER) "  -> interaction kept in the outcome model\n"
    else       sprintf("  -> no interaction (p >= %.2f): no A x M term\n", PCUT))

savediag("S3", "G-formula ignoring MOC (manual)", "p interaction AxM", p_inter)

gf_med <- function(z) {

  ## Step 1
  rhs <- if (INTER) c(paste0(M, "*", A), W) else c(M, A, W)
  my  <- lm(f(Y, rhs),      z)
  mm  <- lm(f(M, c(A, W)),  z)

  ## Step 2.1  counterfactual mediator under each regimen
  Mt <- predict(mm, setA(z, 1))
  Mc <- predict(mm, setA(z, 0))

  ## Step 2.2  the four counterfactual outcomes
  yam <- function(a, mv) { q <- setA(z, a); q[[M]] <- mv; mean(predict(my, q)) }
  y11 <- yam(1, Mt); y10 <- yam(1, Mc)
  y01 <- yam(0, Mt); y00 <- yam(0, Mc)

  ## Step 3  the mediation formula
  NIE <- y01 - y00; TIE <- y11 - y10
  NDE <- y10 - y00; TDE <- y11 - y01

  c(ATE = y11 - y00,
    NIE = NIE, TIE = TIE, AIE = (NIE + TIE) / 2,
    NDE = NDE, TDE = TDE, ADE = (NDE + TDE) / 2,
    Mean_Y1M1 = y11, Mean_Y1M0 = y10, Mean_Y0M1 = y01, Mean_Y0M0 = y00)
}

bt2 <- boot_pct(gf_med, d)
print(bt2, digits = 3)

savebs("S3", "G-formula ignoring MOC (manual)", bt2)


#########################################################################
############## 2.1 - SUPPLEMENTARY BOX S2: mediation    #################
##############               (R mediation package)      #################
##############                          -> Table S3     #################
#########################################################################

model.m2 <- lm(f(M, c(A, W)), data = d)
model.y2 <- lm(f(Y, if (INTER) c(paste0(M, "*", A), W) else c(M, A, W)), data = d)

set.seed(SEED)
out.2 <- mediate(model.m2, model.y2, sims = NBOOT, treat = A, mediator = M,
                 control.value = 0, treat.value = 1, INT = INTER)
summary(out.2)

savest("S3", "mediation", "ACME_control", out.2$d0, sd(out.2$d0.sims), out.2$d0.ci[1], out.2$d0.ci[2])
savest("S3", "mediation", "ACME_treated", out.2$d1, sd(out.2$d1.sims), out.2$d1.ci[1], out.2$d1.ci[2])
savest("S3", "mediation", "ADE_control",  out.2$z0, sd(out.2$z0.sims), out.2$z0.ci[1], out.2$z0.ci[2])
savest("S3", "mediation", "ADE_treated",  out.2$z1, sd(out.2$z1.sims), out.2$z1.ci[1], out.2$z1.ci[2])
savest("S3", "mediation", "AIE",   out.2$d.avg, sd(out.2$d.avg.sims), out.2$d.avg.ci[1], out.2$d.avg.ci[2])
savest("S3", "mediation", "ADE",   out.2$z.avg, sd(out.2$z.avg.sims), out.2$z.avg.ci[1], out.2$z.avg.ci[2])
savest("S3", "mediation", "Total", out.2$tau.coef, sd(out.2$tau.sims), out.2$tau.ci[1], out.2$tau.ci[2])

save_sheets()


#########################################################################
############## 3 - SUPPLEMENTARY BOX S2.1: two ordered  #################
##############      mediators, CRP -> ASDAS-PRO -> BASFI ################
##############                    -> Tables S4 and S5    ################
#########################################################################

#===============================================================================
# Edge g-formula. The nested counterfactuals are built in causal order:
#   M1 (CRP) depends on treatment; M2 (ASDAS-PRO) depends on treatment and M1;
#   BASFI depends on treatment and both mediators.
#
#   Y_000  everyone untreated, both mediators at their untreated values
#   Y_100  outcome sees treatment, mediators still untreated -> direct effect
#   Y_mid  M1 treated, M2 still at its fully untreated value
#   Y_110  M1 treated and M2 responds to M1, but not to treatment directly
#   Y_111  everything treated
#===============================================================================

edge_gf <- function(z, extra = NULL) {

  ## Step 1  models in causal order. extra carries the interaction terms, and
  ## goes into the outcome model only - the two mediator models never take one,
  ## exactly as in run 3.2 of Analysis.do.
  my  <- lm(f(Y,  c(A, M1, M2, W2, bM1, bM2, extra)), z)   # outcome
  mm2 <- lm(f(M2, c(A, M1,     W2, bM1, bM2)), z)   # mediator 2
  mm1 <- lm(f(M1, c(A,         W2, bM1, bM2)), z)   # mediator 1

  ## Step 2.1  counterfactual M1
  M1_1 <- predict(mm1, setA(z, 1))
  M1_0 <- predict(mm1, setA(z, 0))

  ## Step 2.2  counterfactual M2 at the three combinations needed
  m2at <- function(a, m1v) { q <- setA(z, a); q[[M1]] <- m1v; predict(mm2, q) }
  M2_A0_M1_0 <- m2at(0, M1_0)
  M2_A0_M1_1 <- m2at(0, M1_1)
  M2_A1_M1_1 <- m2at(1, M1_1)

  ## Step 2.3  counterfactual BASFI at the nested combinations
  yat <- function(a, m1v, m2v) {
    q <- setA(z, a); q[[M1]] <- m1v; q[[M2]] <- m2v; mean(predict(my, q))
  }
  m000 <- yat(0, M1_0, M2_A0_M1_0)
  m100 <- yat(1, M1_0, M2_A0_M1_0)
  mmid <- yat(1, M1_1, M2_A0_M1_0)
  m110 <- yat(1, M1_1, M2_A0_M1_1)
  m111 <- yat(1, M1_1, M2_A1_M1_1)

  ## Step 3  path-specific effects
  c(Direct        = m100 - m000,   # A -> Y
    PSE_M1        = m110 - m100,   # A -> M1 ~> Y (combined)
    PSE_M2        = m111 - m110,   # A -> M2 -> Y
    Total         = m111 - m000,
    PSE_M1_direct = mmid - m100,   # A -> M1 -> Y, skipping M2
    PSE_M1_via_M2 = m110 - mmid,   # A -> M1 -> M2 -> Y
    ## Step 4  effects on the mediators themselves, not on Y
    ATE_M1   = mean(M1_1) - mean(M1_0),
    ATE_M2   = mean(M2_A1_M1_1) - mean(M2_A0_M1_0),
    M1_on_M2 = mean(M2_A0_M1_1) - mean(M2_A0_M1_0))
}

#####>>>>>> Step 0  screen the four candidate interactions one at a time

it <- function(fml, term) {
  m <- lm(fml, d); summary(m)$coefficients[term, "Pr(>|t|)"]
}
p1 <- it(f(Y,  c(A, M1, M2, W2, bM1, bM2, paste0(A, ":", M1))), paste0(A, ":", M1))
p2 <- it(f(Y,  c(A, M1, M2, W2, bM1, bM2, paste0(A, ":", M2))), paste0(A, ":", M2))
p3 <- it(f(Y,  c(A, M1, M2, W2, bM1, bM2, paste0(M1, ":", M2))), paste0(M1, ":", M2))
p4 <- it(f(M2, c(A, M1,     W2, bM1, bM2, paste0(A, ":", M1))), paste0(A, ":", M1))

cat(sprintf("\nInteraction tests (significant at p < %.2f)\n", PCUT))
pline <- function(i, lab, p) cat(sprintf("%d. %-46s p = %.4f\n", i, lab, p))
pline(1, paste(A,  "x", M1, "in the outcome model"),  p1)
pline(2, paste(A,  "x", M2, "in the outcome model"),  p2)
pline(3, paste(M1, "x", M2, "in the outcome model"),  p3)
pline(4, paste(A,  "x", M1, "in the mediator model"), p4)

#####>>>>>> The two that pass the screen, refit together. Both product terms
#####>>>>>> contain M1, so they compete for the same variance and each p-value
#####>>>>>> rises; the joint test is what says whether there is a signal at all.

INT2 <- c(paste0(A, ":", M1), paste0(M1, ":", M2))

both <- lm(f(Y, c(A, M1, M2, W2, bM1, bM2, INT2)), d)
p5 <- summary(both)$coefficients[INT2[1], "Pr(>|t|)"]
p6 <- summary(both)$coefficients[INT2[2], "Pr(>|t|)"]
p7 <- anova(lm(f(Y, c(A, M1, M2, W2, bM1, bM2)), d), both)[2, "Pr(>F)"]

cat("\nBoth retained interactions in one model\n")
pline(5, paste(A,  "x", M1), p5)
pline(6, paste(M1, "x", M2), p6)
pline(7, "joint test of the two", p7)

EGF <- "Edge g-formula: all paths with interaction (manual)"
savediag("S5", EGF, "p interaction AxM1 outcome model",  p1)
savediag("S5", EGF, "p interaction AxM2 outcome model",  p2)
savediag("S5", EGF, "p interaction M1xM2 outcome model", p3)
savediag("S5", EGF, "p interaction AxM1 mediator model", p4)
savediag("S5", EGF, "p interaction AxM1 both in model",  p5)
savediag("S5", EGF, "p interaction M1xM2 both in model", p6)
savediag("S5", EGF, "p interaction joint test of the two", p7)


#####>>>>>> Step 1 to 4  the decomposition, both specifications, one seed, so
#####>>>>>> the two sets of rows are computed on the same bootstrap resamples

bt3  <- boot_pct(function(z) edge_gf(z),       d)   # without the interactions
bt3i <- boot_pct(function(z) edge_gf(z, INT2), d)   # with both interactions
print(bt3,  digits = 3)
print(bt3i, digits = 3)

#####>>>>>> Table S4 takes the four main paths, Table S5 the whole set, each
#####>>>>>> with and without the interactions so the two sit in one table

main <- c("Direct", "PSE_M1", "PSE_M2", "Total")
savebs("S4", "Edge g-formula: main paths without interaction (manual)", bt3 [bt3$effect  %in% main, ])
savebs("S4", "Edge g-formula: main paths with interaction (manual)",    bt3i[bt3i$effect %in% main, ])
savebs("S5", "Edge g-formula: all paths without interaction (manual)",  bt3)
savebs("S5", "Edge g-formula: all paths with interaction (manual)",     bt3i)


#########################################################################
############## 3.1 - SUPPLEMENTARY BOX S2.1                #############
##############               (R paths package)             #############
##############                     -> Tables S4 and S5     #############
#########################################################################

#===============================================================================
# The paths package fits the same nested outcome models and reports the pure
# (Type I and Type II) decompositions. Without interactions, both target the
# same effects as the manual code above.
#===============================================================================

fit_m0 <- lm(f(Y, c(A,         W2, bM1, bM2)), data = d)
fit_m1 <- lm(f(Y, c(A, M1,     W2, bM1, bM2)), data = d)
fit_m2 <- lm(f(Y, c(A, M1, M2, W2, bM1, bM2)), data = d)

set.seed(SEED)
paths_out <- paths(a = A, y = Y, m = list(M1, M2),
                   models = list(fit_m0, fit_m1, fit_m2),
                   data = d, nboot = NBOOT, conf_level = 0.95)
summary(paths_out)

## summary.paths returns the four decompositions in $estimates. The pure
## Type I table is the one the manual code above computes; the rownames carry
## the effect labels and the columns are Estimate, Std. Err., lower CI, upper CI,
## P-value. Both pure decompositions (Type I and Type II) are saved.
ests <- summary(paths_out)$estimates
for (dec in c("pure_t1", "pure_t2")) {
  pe <- ests[[dec]]
  lab <- paste("paths package without interaction:", sub("_", " ", dec))
  for (i in seq_len(nrow(pe)))
    savest("S5", lab, rownames(pe)[i], pe[i, 1], pe[i, 2], pe[i, 3], pe[i, 4])
}

save_sheets()


#########################################################################
############## 4 - SUPPLEMENTARY BOX S6: mediation       ################
##############        accounting for MOC (manual)        ################
##############                         -> Table S10       ################
#########################################################################

#===============================================================================
# The three post-treatment mediator-outcome confounders are simulated forward
# from treatment, then the mediator is simulated given treatment and the
# simulated MOC, then the outcome. Linear predictions of the probability of each
# MOC (conditional expectations) are used instead of stochastic draws, to reduce
# Monte-Carlo variance, as stated under Supplementary Box S6.
#===============================================================================

p_interMOC <- summary(lm(f(Y, c(M, A, MOC, W, paste0(M, ":", A))), d)
                      )$coefficients[paste0(M, ":", A), "Pr(>|t|)"]
cat(sprintf("\nTreatment-mediator interaction (%s x %s): p = %.4f\n", A, M, p_interMOC))
savediag("S10", "G-formula with MOC (manual)", "p interaction AxM", p_interMOC)

gf_moc <- function(z) {

  ## Step 1
  my  <- lm(f(Y, c(M, A, MOC, W)), z)                       # outcome
  mm  <- lm(f(M, c(A, MOC, W)),    z)                       # mediator
  ## emmtvt2 is rare (about 7% of patients), so some bootstrap resamples are
  ## perfectly separated and glm warns that fitted probabilities hit 0 or 1.
  ## The warnings are expected and harmless here: the conditional expectation
  ## is used as it stands, never as the probability of a draw. Note that Stata's
  ## logit instead DROPS perfectly predicted observations.
  
  mc  <- lapply(MOC, function(v) glm(f(v, c(A, W)), binomial, z))  # the three MOC
  names(mc) <- MOC

  ## Step 2.1  counterfactual MOC, as conditional expectations
  moc_at <- function(a) {
    q <- setA(z, a)
    for (v in MOC) q[[v]] <- predict(mc[[v]], q, type = "response")
    q
  }
  q1 <- moc_at(1); q0 <- moc_at(0)

  ## Step 2.2  counterfactual mediator, each with its own MOC
  Mt <- predict(mm, q1)
  Mc <- predict(mm, q0)

  ## Step 2.3  the four counterfactual outcomes; the MOC take the value of the
  ## mediator's regime (am), as in Box S6 and Analysis.do
  yam <- function(a, am) {
    q <- if (am == 1) q1 else q0
    q <- setA(q, a); q[[M]] <- if (am == 1) Mt else Mc
    mean(predict(my, q))
  }
  y11 <- yam(1, 1); y10 <- yam(1, 0)
  y01 <- yam(0, 1); y00 <- yam(0, 0)

  NIE <- y01 - y00; TIE <- y11 - y10
  NDE <- y10 - y00; TDE <- y11 - y01

  c(ATE = y11 - y00,
    NIE = NIE, TIE = TIE, AIE = (NIE + TIE) / 2,
    NDE = NDE, TDE = TDE, ADE = (NDE + TDE) / 2,
    Mean_Y1M1 = y11, Mean_Y1M0 = y10, Mean_Y0M1 = y01, Mean_Y0M0 = y00)
}

bt4 <- boot_pct(gf_moc, d)
print(bt4, digits = 3)

savebs("S10", "G-formula with MOC (manual)", bt4)
save_sheets()



#########################################################################
############## 5 - SUPPLEMENTARY BOX S3: marginal        ################
##############        structural model with IPTW         ################
##############                         -> Table S6       ################
#########################################################################

#===============================================================================
# Step 1  denominator and numerator treatment models
# Step 2  stabilised IPTW = numerator / denominator
# Step 3  weighted regression of BASFI on bDMARD, the marginal structural model
# Step 4  HC1 robust standard error, matching Stata's "regress ..., robust"
#===============================================================================

##### Step 1
den <- glm(f(A, W), binomial, d)
num <- glm(f(A, "1"), binomial, d)

##### Step 2
ps  <- predict(den, type = "response")
nps <- predict(num, type = "response")

d$denom  <- ps  * d[[A]] + (1 - ps)  * (1 - d[[A]])
d$numer  <- nps * d[[A]] + (1 - nps) * (1 - d[[A]])
d$uiptw  <- 1 / d$denom
d$siptw  <- d$numer / d$denom

diagw <- function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))
DIAG  <- rbind("Propensity score" = diagw(ps),
               "uIPTW"            = diagw(d$uiptw),
               "sIPTW"            = diagw(d$siptw))
cat("\nIPTW diagnostics\n"); print(round(DIAG, 3))

##### Step 3 and 4  the MSM in the whole population
msm  <- lm(f(Y, A), d, weights = d$siptw)
b    <- coef(msm)[A]
se   <- robust_se(msm, d$siptw)[A]
ev5  <- evalue_md(b)

cat(sprintf("\nMSM: %.3f (%.3f; %.3f), robust SE %.3f, e-value %.3f\n",
            b, b - 1.96 * se, b + 1.96 * se, se, ev5))

MSMm <- "MSM with sIPTW (manual)"
savest("S6", MSMm, "ATE", b, se, b - 1.96 * se, b + 1.96 * se, ev5, NA)
savediag("S6", MSMm, "sIPTW mean", mean(d$siptw))
savediag("S6", MSMm, "sIPTW min",  min(d$siptw))
savediag("S6", MSMm, "sIPTW max",  max(d$siptw))
savediag("S6", MSMm, "sIPTW SD",   sd(d$siptw))

##### Positivity: area of common support (trimming)
acsmax <- min(max(ps[d[[A]] == 1]), max(ps[d[[A]] == 0]))
acsmin <- max(min(ps[d[[A]] == 1]), min(ps[d[[A]] == 0]))
acs    <- ps >= acsmin & ps <= acsmax
n_out  <- sum(!acs)
pct_out <- 100 * mean(!acs)
cat(sprintf("ACS bounds: %.4f to %.4f; outside: %d (%.1f%%)\n",
            acsmin, acsmax, n_out, pct_out))

msm_acs <- lm(f(Y, A), d[acs, ], weights = d$siptw[acs])
b_a  <- coef(msm_acs)[A]
se_a <- robust_se(msm_acs, d$siptw[acs])[A]

MSMa <- "MSM with sIPTW within ACS (manual)"
savest("S6", MSMa, "ATE", b_a, se_a, b_a - 1.96 * se_a, b_a + 1.96 * se_a, NA, NA,
       n = sum(acs))
savediag("S6", MSMa, "N outside ACS",   n_out)
savediag("S6", MSMa, "pct outside ACS", pct_out)

##### Positivity: truncation at the 1st and 99th percentiles
qs <- quantile(d$siptw, c(0.01, 0.99), type = 2)
d$siptw_t <- pmin(pmax(d$siptw, qs[1]), qs[2])
n_trunc   <- sum(d$siptw != d$siptw_t)

msm_tr <- lm(f(Y, A), d, weights = d$siptw_t)
b_t  <- coef(msm_tr)[A]
se_t <- robust_se(msm_tr, d$siptw_t)[A]

MSMt <- "MSM with sIPTW truncated P1-P99 (manual)"
savest("S6", MSMt, "ATE", b_t, se_t, b_t - 1.96 * se_t, b_t + 1.96 * se_t, NA, NA)
savediag("S6", MSMt, "N truncated sIPTW", n_trunc)

cat(sprintf("within ACS: %.3f (SE %.3f); truncated (%d): %.3f (SE %.3f)\n",
            b_a, se_a, n_trunc, b_t, se_t))

##### Balance before and after weighting (Supplementary Table S8)
##### Means of the baseline confounders in the treated and the untreated, in the
##### original population and in the pseudo-population created by the sIPTW.
##### SMD as in Programs/msmbalanceV1.do (Jackson, Epidemiology 2016): the
##### unweighted and the weighted difference in means are divided by the same
##### unweighted pooled SD, sqrt((var treated + var untreated)/2), with
##### var = p(1-p) for binary variables.

BAL <- c("basfitotalt1", "asdastotalt1", "age", "sex", "comorbbin", "mny",
         "asasmri", "hla", "pertvt1", "ibdbl", "emmtvt1", "comedtvt1")

tr <- d[[A]] == 1
un <- d[[A]] == 0

balance_row <- function(v) {
  x  <- d[[v]]
  m1 <- mean(x[tr])
  m0 <- mean(x[un])
  if (all(x %in% c(0, 1))) {
    s1 <- m1 * (1 - m1)
    s0 <- m0 * (1 - m0)
  } else {
    s1 <- var(x[tr])
    s0 <- var(x[un])
  }
  den <- sqrt((s1 + s0) / 2)
  wm1 <- weighted.mean(x[tr], d$siptw[tr])
  wm0 <- weighted.mean(x[un], d$siptw[un])
  data.frame(variable        = v,
             mean_treated    = m1,
             mean_untreated  = m0,
             smd_unweighted  = (m1 - m0) / den,
             wmean_treated   = wm1,
             wmean_untreated = wm0,
             smd_weighted    = (wm1 - wm0) / den,
             stringsAsFactors = FALSE)
}

S8 <- do.call(rbind, lapply(BAL, balance_row))

## Group sizes for the column headers: patients in the original population and
## sum of the sIPTW in the pseudo-population
S8n <- data.frame(group      = c("treated", "untreated", "total"),
                  n_original = c(sum(tr), sum(un), nrow(d)),
                  n_pseudo   = c(sum(d$siptw[tr]), sum(d$siptw[un]), sum(d$siptw)),
                  mean_siptw = c(NA, NA, mean(d$siptw)),
                  stringsAsFactors = FALSE)

cat("\nSupplementary Table S8: original population vs pseudo-population (sIPTW)\n")
print(S8, digits = 3)
print(S8n, digits = 4)

write_xlsx(list(balance = S8, groups = S8n),
           paste0(TABLES, "Supplementary_Table_S8_r.xlsx"))


#########################################################################
############## 5.1 - SUPPLEMENTARY BOX S3: MSM with IPTW ################
##############           (R ipw and survey packages)     ################
##############                         -> Table S6       ################
#########################################################################

#===============================================================================
# The package counterpart of the manual code above, and of Stata's teffects ipw.
#
# ipwpoint() fits the same two treatment models and returns numerator/denominator
# for a point exposure, so its weights are the stabilised IPTW of Step 2 - they
# come out identical to d$siptw, which the check below asserts.
#
# svyglm() then fits the weighted regression and reports a sandwich standard
# error, treating the weights as known. It divides by n-1 where Stata's robust
# option divides by n-k.
#
# ipwpoint deparses its denominator argument instead of evaluating it, so the
# formula has to be written out here rather than built from W. The stopifnot
# below is what guards against the two drifting apart: the weights are compared
# with d$siptw, which IS built from W, so any divergence stops the script.
#===============================================================================

w_ipw <- ipwpoint(exposure    = bionew,
                  family      = "binomial",
                  link        = "logit",
                  numerator   = ~ 1,
                  denominator = ~ age + sex + comorbbin + mny + asasmri + hla +
                                  pertvt1 + ibdbl + emmtvt1 + comedtvt1 +
                                  asdastotalt1 + basfitotalt1,
                  data        = d)

d$siptw_ipw <- w_ipw$ipw.weights
stopifnot(isTRUE(all.equal(d$siptw_ipw, unname(d$siptw))))

des <- svydesign(ids = ~ 1, weights = ~ siptw_ipw, data = d)
msm_pkg <- svyglm(f(Y, A), design = des)

b_p  <- coef(msm_pkg)[A]
se_p <- sqrt(diag(vcov(msm_pkg)))[A]
ci_p <- confint(msm_pkg, A)

cat(sprintf("\nMSM (ipw + survey): %.3f (%.3f; %.3f), SE %.3f\n",
            b_p, ci_p[1], ci_p[2], se_p))

MSMp <- "MSM with sIPTW (ipw and survey)"
savest("S6", MSMp, "ATE", b_p, se_p, ci_p[1], ci_p[2], evalue_md(b_p), NA)
savediag("S6", MSMp, "sIPTW mean", mean(d$siptw_ipw))
savediag("S6", MSMp, "sIPTW min",  min(d$siptw_ipw))
savediag("S6", MSMp, "sIPTW max",  max(d$siptw_ipw))

save_sheets()


#########################################################################
############## 6 - SUPPLEMENTARY BOX S4: targeted        ################
##############        maximum likelihood estimation      ################
##############                         -> Table S7       ################
#########################################################################

#===============================================================================
# Step 1  outcome model and treatment model
# Step 2  clever covariate, fluctuation model, targeted predictions
# Step 3  the ATE is the difference of the two targeted means
# Step 4  variance from the efficient influence function
#===============================================================================

##### Step 1
qmod <- lm(f(Y, c(A, W)), d)
gmod <- glm(f(A, W), binomial, d)

##### Step 2.1  predicted outcomes under each regimen
Y1   <- predict(qmod, setA(d, 1))
Y0   <- predict(qmod, setA(d, 0))
Yobs <- predict(qmod, d)

##### Step 2.2  propensity score
PS <- predict(gmod, type = "response")

##### Step 2.3  clever covariate
H  <- d[[A]] / PS - (1 - d[[A]]) / (1 - PS)

##### Step 2.4  targeting step
eps <- coef(lm(d[[Y]] ~ offset(Yobs) + H - 1))["H"]
Y1_t   <- Y1   + eps / PS
Y0_t   <- Y0   - eps / (1 - PS)
Yobs_t <- Yobs + eps * H

##### Step 3
ATE_tmle <- mean(Y1_t) - mean(Y0_t)

##### Step 4  efficient influence function
EIF    <- H * (d[[Y]] - Yobs_t) + Y1_t - Y0_t - ATE_tmle
se_ATE <- sqrt(sum(EIF^2) / nrow(d)^2)
lb <- ATE_tmle - 1.96 * se_ATE
ub <- ATE_tmle + 1.96 * se_ATE
ev6 <- evalue_md(ATE_tmle)

cat(sprintf("\nTMLE: %.3f (%.3f; %.3f), SE %.3f, e-value %.3f\n",
            ATE_tmle, lb, ub, se_ATE, ev6))

TM <- "TMLE (manual)"
savest("S7", TM, "ATE", ATE_tmle, se_ATE, lb, ub, ev6, NA)
savediag("S7", TM, "pct g-truncated",    100 * mean(PS < 0.025 | PS > 0.975))
savediag("S7", TM, "SD influence curve", sd(EIF))
savediag("S7", TM, "Mean Y1", mean(Y1_t))
savediag("S7", TM, "Mean Y0", mean(Y0_t))


#########################################################################
############## 6.1 - SUPPLEMENTARY BOX S4: TMLE          ################
##############               (R ltmle package)           ################
##############                         -> Table S7       ################
#########################################################################

#===============================================================================
# Column order is the causal order: ltmle takes every column to the left of a
# node as its parents. No post-treatment node is included, so this is the
# time-fixed total effect of Supplementary Box S4.
#
# ltmle rescales a continuous outcome to [0, 1] and fluctuates on the logit
# scale, while the manual code in Box S4 fluctuates linearly. The two therefore
# do not have to always agree. The manual code is the one that reproduces the
# Stata manual TMLE exactly.
#===============================================================================

tmledata <- d[, c(W, A, Y)]

Qform <- setNames(paste("Q.kplus1 ~", paste(c(W, A), collapse = " + ")), Y)
gform <- setNames(paste(A, "~", paste(W, collapse = " + ")), A)

set.seed(SEED)
ATEl <- ltmle(tmledata, Anodes = A, Lnodes = NULL, Ynodes = Y,
              Qform = Qform, gform = gform,
              abar = list(1, 0), estimate.time = FALSE)
sl <- summary(ATEl)
print(sl)

e <- sl$effect.measures$ATE
savest("S7", "ltmle", "ATE", e$estimate, e$std.dev, e$CI[1], e$CI[2],
       evalue_md(e$estimate), NA)

##### Diagnostics, computed manually after ltmle
psl <- ATEl$cum.g[, dim(ATEl$cum.g)[2], 1]
ICl <- ATEl$IC
if (is.list(ICl)) ICl <- ICl[[1]]
if (is.matrix(ICl) && ncol(ICl) == 2) ICl <- ICl[, 1] - ICl[, 2]

savediag("S7", "ltmle", "pct g-truncated",    100 * mean(psl < 0.025 | psl > 0.975))
savediag("S7", "ltmle", "SD influence curve", sd(ICl))
savediag("S7", "ltmle", "Mean Y1", sl$effect.measures$treatment$estimate)
savediag("S7", "ltmle", "Mean Y0", sl$effect.measures$control$estimate)

save_sheets()


#########################################################################
############## 7 - SUPPLEMENTARY TABLES S9 AND S9.1:     ################
##############      mediation within levels of baseline  ################
##############      characteristics (manual)             ################
#########################################################################

#===============================================================================
# The R counterpart of run 7 of Analysis.do.
#   Table S9    one mediator (ASDAS): gf_med() of section 2, within each stratum
#   Table S9.1  two mediators (CRP, then ASDAS-PRO): edge_gf() of section 3,
#               main paths without interaction, within each stratum
# The interaction p-value tests effect modification of the total effect: the
# bDMARD x characteristic term in the outcome model without the mediator.
# Within a stratum the characteristic is constant, so it is removed from W and
# W2 while the stratum is analysed, and W and W2 are set back afterwards.
#===============================================================================

if (!exists("gf_med") || !exists("INTER") || !exists("edge_gf"))
  stop("Run sections 2 and 3 first: gf_med(), INTER and edge_gf() are defined there")

W_main  <- W
W2_main <- W2

S9  <- NULL
S91 <- NULL

for (s in c("crpelevatedt1", "asasmri", "mny")) {

  ## Interaction with treatment on the total effect (no mediator in the model)
  rhs    <- unique(c(A, W_main, s, paste0(A, ":", s)))
  p_int  <- summary(lm(f(Y, rhs), d))$coefficients[paste0(A, ":", s), "Pr(>|t|)"]
  cat(sprintf("\nInteraction %s x %s (total effect): p = %.4f\n", s, A, p_int))

  for (lv in c(0, 1)) {

    z  <- d[d[[s]] == lv, ]
    W  <- setdiff(W_main,  s)
    W2 <- setdiff(W2_main, s)

    b1 <- boot_pct(gf_med, z)                   # one mediator (Table S9)
    b2 <- boot_pct(function(q) edge_gf(q), z)   # two mediators (Table S9.1)

    g1 <- function(e) b1[b1$effect == e, ]
    g2 <- function(e) b2[b2$effect == e, ]

    S9 <- rbind(S9, data.frame(
      subgroup = s, level = lv, N = nrow(z), REPS = min(b1$reps),
      ATE = g1("ATE")$estimate, ATE_LCL = g1("ATE")$lci, ATE_UCL = g1("ATE")$uci,
      ADE = g1("ADE")$estimate, ADE_LCL = g1("ADE")$lci, ADE_UCL = g1("ADE")$uci,
      AIE = g1("AIE")$estimate, AIE_LCL = g1("AIE")$lci, AIE_UCL = g1("AIE")$uci,
      P_INTER = p_int))

    S91 <- rbind(S91, data.frame(
      subgroup = s, level = lv, N = nrow(z), REPS = min(b2$reps),
      TOTAL   = g2("Total")$estimate,  TOTAL_LCL   = g2("Total")$lci,  TOTAL_UCL   = g2("Total")$uci,
      DIRECT  = g2("Direct")$estimate, DIRECT_LCL  = g2("Direct")$lci, DIRECT_UCL  = g2("Direct")$uci,
      PSE_CRP = g2("PSE_M1")$estimate, PSE_CRP_LCL = g2("PSE_M1")$lci, PSE_CRP_UCL = g2("PSE_M1")$uci,
      PSE_PRO = g2("PSE_M2")$estimate, PSE_PRO_LCL = g2("PSE_M2")$lci, PSE_PRO_UCL = g2("PSE_M2")$uci,
      ATE_CRP = g2("ATE_M1")$estimate, ATE_CRP_LCL = g2("ATE_M1")$lci, ATE_CRP_UCL = g2("ATE_M1")$uci,
      ATE_PRO = g2("ATE_M2")$estimate, ATE_PRO_LCL = g2("ATE_M2")$lci, ATE_PRO_UCL = g2("ATE_M2")$uci))
  }
}

W  <- W_main    # back to the main confounder sets
W2 <- W2_main

print(S9,  digits = 3)
print(S91, digits = 3)

write_xlsx(S9,  paste0(TABLES, "Supplementary_Table_S9_r.xlsx"))
write_xlsx(S91, paste0(TABLES, "Supplementary_Table_S9_1_r.xlsx"))



#<<<<<<<<<<<########################################################>>>>>>>>>>>#
#<<<<<<<<<<<##########   SENSITIVITY ANALYSES (Table 2,     ########>>>>>>>>>>>#
#<<<<<<<<<<<##########   Figure 3 and Box S7)               ########>>>>>>>>>>>#
#<<<<<<<<<<<########################################################>>>>>>>>>>>#

#===============================================================================
# The R counterpart of runs 8 to 8.4 of Analysis.do. Every estimate goes to the
# sheet T2 of estimates_r.xlsx, with the same estimator labels as the sheet T2
# of estimates_stata.xlsx.
#
#   8    Complete vs missing 6-month data; positivity of being observed (n = 481)
#   8.1  Treatment = at least 1 month of bDMARD: single and two mediators (n = 419)
#   8.2  Treatment = at least 3 months of bDMARD: single and two mediators (n = 419)
#   8.3  Total effect, MSM with sIPTW x sIPCW, manual        (n = 481, 419 observed)
#   8.3.1 Total effect, MSM with censoring weights, ipw + survey packages
#   8.4  Total effect, TMLE with censoring, manual           (n = 481, 419 observed)
#   8.4.1 Total effect, TMLE with censoring, ltmle package
#
# Runs 8.1 and 8.2 reuse gf_med() (Box S2) and edge_gf() (Box S2.1) defined
# above: only the treatment variable changes, so A is switched to the new
# definition, the two functions are run, and A is set back to bionew.
#===============================================================================


#########################################################################
############## 8 - SENSITIVITY: complete vs missing      ################
##############      6-month data, n = 481                ################
##############                         -> sheet T2       ################
#########################################################################

#===============================================================================
# Eligible patients with complete baseline covariates (W): n = 481.
# Censored = ASDAS or BASFI missing at 6 months (62), so the 419 observed are
# exactly the patients of the main analysis: runs 8.3 and 8.4 reweight them to
# represent all 481.
#===============================================================================

e <- read.csv(paste0(DATA, "originalfulllong800.csv"), header = TRUE, sep = ",")
e <- e[e$t == 6, ]
e <- e[complete.cases(e[, W]), ]
stopifnot(nrow(e) == 481)
N_ELIG <- nrow(e)

## complete ASDAS and BASFI at 6 months (the analysis cohort)
e$complete6 <- 0
e$complete6[!is.na(e$asdastotalt2) & !is.na(e$basfitotalt2)] <- 1
stopifnot(sum(e$complete6) == 419)

## ASDAS and BASFI observed at 6 months (the uncensored: the main analysis)
e$uncensored <- 0
e$uncensored[e$complete6 == 1] <- 1
stopifnot(sum(e$uncensored == 0) == 62)

## ibdbl is left out of the censoring model: every patient with baseline IBD has
## ASDAS and BASFI at 6 months, so it predicts being observed perfectly (Stata's logit
## drops it). This holds while no IBD patient is censored.
stopifnot(sum(e$uncensored == 0 & e$ibdbl == 1) == 0)
Wc <- setdiff(W, "ibdbl")

## Positivity of complete 6-month data, P(complete | treatment, baseline)
pc6 <- predict(glm(f("complete6", c(A, W)), binomial, e), type = "response")

cat(sprintf("\nEligible %d; complete 6-month data %d; missing %d (treated %d, untreated %d)\n",
            N_ELIG, sum(e$complete6), sum(e$complete6 == 0),
            sum(e$complete6 == 0 & e[[A]] == 1), sum(e$complete6 == 0 & e[[A]] == 0)))
cat(sprintf("Min P(complete): treated %.3f, untreated %.3f\n",
            min(pc6[e[[A]] == 1]), min(pc6[e[[A]] == 0])))

MD <- "Missing 6-month data (descriptive)"
savediag("T2", MD, "N eligible with complete baseline", N_ELIG)
savediag("T2", MD, "N complete 6-month data",           sum(e$complete6))
savediag("T2", MD, "N missing 6-month data",            sum(e$complete6 == 0))
savediag("T2", MD, "N missing among treated",           sum(e$complete6 == 0 & e[[A]] == 1))
savediag("T2", MD, "N missing among untreated",         sum(e$complete6 == 0 & e[[A]] == 0))
savediag("T2", MD, "Min P(complete) treated",           min(pc6[e[[A]] == 1]))
savediag("T2", MD, "Min P(complete) untreated",         min(pc6[e[[A]] == 0]))


#########################################################################
############## 8.1 - SENSITIVITY: at least 1 month of    ################
##############        bDMARD, single and two mediators   ################
##############                         -> sheet T2       ################
#########################################################################

#===============================================================================
# Same models and specification as the main analysis (no interaction), only
# the treatment variable changes: bionew1m = at least 30 days of exposure.
#===============================================================================

## gf_med() and INTER come from section 2, edge_gf() from section 3
if (!exists("gf_med") || !exists("INTER") || !exists("edge_gf"))
  stop("Run sections 2 and 3 first: gf_med(), INTER and edge_gf() are defined there")

A <- "bionew1m"
stopifnot(sum(d[[A]]) == 112)

bt81s <- boot_pct(gf_med, d)                    # single mediator (Box S2)
bt81p <- boot_pct(function(z) edge_gf(z), d)    # two mediators (Box S2.1)
print(bt81s, digits = 3)
print(bt81p, digits = 3)

S81 <- "Single mediator - bDMARD at least 1 month (manual)"
P81 <- "Two mediators - bDMARD at least 1 month (manual)"
savebs("T2", S81, bt81s[bt81s$effect %in% c("ATE", "AIE", "ADE", "Mean_Y1M1", "Mean_Y0M0"), ])
savebs("T2", P81, bt81p[bt81p$effect %in% c("Direct", "PSE_M1", "PSE_M2", "Total"), ])
savediag("T2", S81, "N treated",   sum(d[[A]] == 1))
savediag("T2", S81, "N untreated", sum(d[[A]] == 0))
savediag("T2", P81, "N treated",   sum(d[[A]] == 1))
savediag("T2", P81, "N untreated", sum(d[[A]] == 0))

A <- "bionew" # back to the main definition


#########################################################################
############## 8.2 - SENSITIVITY: at least 3 months of   ################
##############        bDMARD, single and two mediators   ################
##############                         -> sheet T2       ################
#########################################################################

A <- "bionew3m"
stopifnot(sum(d[[A]]) == 81)

bt82s <- boot_pct(gf_med, d)                    # single mediator (Box S2)
bt82p <- boot_pct(function(z) edge_gf(z), d)    # two mediators (Box S2.1)
print(bt82s, digits = 3)
print(bt82p, digits = 3)

S82 <- "Single mediator - bDMARD at least 3 months (manual)"
P82 <- "Two mediators - bDMARD at least 3 months (manual)"
savebs("T2", S82, bt82s[bt82s$effect %in% c("ATE", "AIE", "ADE", "Mean_Y1M1", "Mean_Y0M0"), ])
savebs("T2", P82, bt82p[bt82p$effect %in% c("Direct", "PSE_M1", "PSE_M2", "Total"), ])
savediag("T2", S82, "N treated",   sum(d[[A]] == 1))
savediag("T2", S82, "N untreated", sum(d[[A]] == 0))
savediag("T2", P82, "N treated",   sum(d[[A]] == 1))
savediag("T2", P82, "N untreated", sum(d[[A]] == 0))

A <- "bionew" # back to the main definition

save_sheets()


#########################################################################
############## 8.3 - SENSITIVITY: total effect, MSM with ################
##############        sIPTW x censoring weights (manual) ################
##############          -> sheet T2, Figure 3, Box S7    ################
#########################################################################

#===============================================================================
# Step 1  treatment models (denominator, numerator) and censoring models
#         (denominator: treatment + baseline; numerator: treatment), all 481
# Step 2  sIPTW = P(A) / P(A | W);  sIPCW = P(observed | A) / P(observed | A, W)
# Step 3  MSM: regression of BASFI on bDMARD in the 419 uncensored, weighted by
#         sIPTW x sIPCW; HC1 robust standard error, as Stata's "robust"
#===============================================================================

##### Step 1
den_a <- glm(f(A, W),                    binomial, e)
num_a <- glm(f(A, "1"),                  binomial, e)
den_c <- glm(f("uncensored", c(A, Wc)),  binomial, e)
num_c <- glm(f("uncensored", A),         binomial, e)

##### Step 2.1  stabilised IPTW
ps_e  <- predict(den_a, type = "response")
nps_e <- predict(num_a, type = "response")
e$siptw <- (nps_e * e[[A]] + (1 - nps_e) * (1 - e[[A]])) /
           (ps_e  * e[[A]] + (1 - ps_e)  * (1 - e[[A]]))

##### Step 2.2  stabilised IPCW (probability of being observed under the
##### treatment actually received)
pc_e    <- predict(den_c, type = "response")
npc_e   <- predict(num_c, type = "response")
e$sipcw <- npc_e / pc_e
e$w_tc  <- e$siptw * e$sipcw

obs <- e$uncensored == 1

cat(sprintf("\nMin P(observed): treated %.3f, untreated %.3f\n",
            min(pc_e[e[[A]] == 1]), min(pc_e[e[[A]] == 0])))
cat(sprintf("sIPCW (uncensored): mean %.3f, min %.3f, max %.3f\n",
            mean(e$sipcw[obs]), min(e$sipcw[obs]), max(e$sipcw[obs])))
cat(sprintf("sIPTW x sIPCW (uncensored): mean %.3f, min %.3f, max %.3f\n",
            mean(e$w_tc[obs]), min(e$w_tc[obs]), max(e$w_tc[obs])))

##### Step 3  MSM in the uncensored
msm_c <- lm(f(Y, A), e[obs, ], weights = e$w_tc[obs])
b_c   <- coef(msm_c)[A]
se_c  <- robust_se(msm_c, e$w_tc[obs])[A]

cat(sprintf("MSM with censoring weights: %.3f (%.3f; %.3f), robust SE %.3f, n = %d (%d observed)\n",
            b_c, b_c - 1.96 * se_c, b_c + 1.96 * se_c, se_c, N_ELIG, sum(obs)))

MC <- "MSM with sIPTW and sIPCW (manual)"
savest("T2", MC, "ATE", b_c, se_c, b_c - 1.96 * se_c, b_c + 1.96 * se_c, NA, NA, n = N_ELIG)
savediag("T2", MC, "N censored",                sum(!obs))
savediag("T2", MC, "N censored treated",        sum(!obs & e[[A]] == 1))
savediag("T2", MC, "N censored untreated",      sum(!obs & e[[A]] == 0))
savediag("T2", MC, "Min P(observed) treated",   min(pc_e[e[[A]] == 1]))
savediag("T2", MC, "Min P(observed) untreated", min(pc_e[e[[A]] == 0]))
savediag("T2", MC, "sIPCW mean",                mean(e$sipcw[obs]))
savediag("T2", MC, "sIPCW max",                 max(e$sipcw[obs]))
savediag("T2", MC, "sIPTW x sIPCW max",         max(e$w_tc[obs]))


#########################################################################
############## 8.3.1 - SENSITIVITY: MSM with censoring   ################
##############          weights (ipw and survey)         ################
##############                         -> sheet T2       ################
#########################################################################

#===============================================================================
# ipwpoint() is called twice: once for treatment (as in section 5.1) and once
# for being observed, with treatment in the numerator model. Its weights are the
# sIPTW and sIPCW of the manual code, which the two checks below assert. As in
# section 5.1, ipwpoint deparses its formula arguments, so they are written out.
# svyglm() then fits the MSM in the uncensored with the product of the weights.
#===============================================================================

w_a <- ipwpoint(exposure    = bionew,
                family      = "binomial",
                link        = "logit",
                numerator   = ~ 1,
                denominator = ~ age + sex + comorbbin + mny + asasmri + hla +
                                pertvt1 + ibdbl + emmtvt1 + comedtvt1 +
                                asdastotalt1 + basfitotalt1,
                data        = e)

w_c <- ipwpoint(exposure    = uncensored,
                family      = "binomial",
                link        = "logit",
                numerator   = ~ bionew,
                denominator = ~ bionew + age + sex + comorbbin + mny + asasmri + hla +
                                pertvt1 + emmtvt1 + comedtvt1 +
                                asdastotalt1 + basfitotalt1,
                data        = e)

stopifnot(isTRUE(all.equal(w_a$ipw.weights,      unname(e$siptw))))
stopifnot(isTRUE(all.equal(w_c$ipw.weights[obs], unname(e$sipcw[obs]))))

e$w_tc_pkg <- w_a$ipw.weights * w_c$ipw.weights

des_c <- svydesign(ids = ~ 1, weights = ~ w_tc_pkg, data = e[obs, ])
msm_c_pkg <- svyglm(f(Y, A), design = des_c)

b_cp  <- coef(msm_c_pkg)[A]
se_cp <- sqrt(diag(vcov(msm_c_pkg)))[A]
ci_cp <- confint(msm_c_pkg, A)

cat(sprintf("\nMSM with censoring weights (ipw + survey): %.3f (%.3f; %.3f), SE %.3f\n",
            b_cp, ci_cp[1], ci_cp[2], se_cp))

savest("T2", "MSM with sIPTW and sIPCW (ipw and survey)", "ATE",
       b_cp, se_cp, ci_cp[1], ci_cp[2], NA, NA, n = N_ELIG)

save_sheets()


#########################################################################
############## 8.4 - SENSITIVITY: total effect, TMLE     ################
##############        with censoring (manual)            ################
##############          -> sheet T2, Figure 3, Box S7    ################
#########################################################################

#===============================================================================
# As Box S4, with three changes:
#   - the outcome model is fitted in the uncensored
#   - the clever covariate carries the probability of being observed:
#       H = C x [A / (PS x pc1) - (1 - A) / ((1 - PS) x pc0)],  C = 1 if observed,
#     so it is zero for the censored
#   - the fluctuation is fitted in the uncensored; the counterfactual
#     predictions and the effect are averaged over all 481
# The censoring model is den_c of run 8.3 (treatment + baseline without ibdbl).
#===============================================================================

##### Step 1  outcome model in the uncensored; PS on all 481
qmod_c <- lm(f(Y, c(A, W)), e[obs, ])

##### Step 2.1  predicted outcomes under each regimen, all 481
Y1c   <- predict(qmod_c, setA(e, 1))
Y0c   <- predict(qmod_c, setA(e, 0))
Yobsc <- predict(qmod_c, e)

##### Step 2.2  PS, and probability of being observed if treated / untreated
PSc <- ps_e
pc1 <- predict(den_c, setA(e, 1), type = "response")
pc0 <- predict(den_c, setA(e, 0), type = "response")

##### Step 2.3  clever covariate (zero for the censored)
Hc <- rep(0, nrow(e))
Hc[obs] <- (e[[A]] / (PSc * pc1) - (1 - e[[A]]) / ((1 - PSc) * pc0))[obs]

##### Step 2.4  targeting step, in the uncensored
eps_c <- coef(lm(e[[Y]][obs] ~ offset(Yobsc[obs]) + Hc[obs] - 1))[1]
Y1c_t   <- Y1c   + eps_c / (PSc * pc1)
Y0c_t   <- Y0c   - eps_c / ((1 - PSc) * pc0)
Yobsc_t <- Yobsc + eps_c * Hc

##### Step 3  the effect, averaged over all 481
ATE_tmle_c <- mean(Y1c_t) - mean(Y0c_t)

##### Step 4  efficient influence function; its first term is zero for the censored
EIFc <- Y1c_t - Y0c_t - ATE_tmle_c
EIFc[obs] <- EIFc[obs] + (Hc * (e[[Y]] - Yobsc_t))[obs]
se_tmle_c <- sqrt(sum(EIFc^2) / nrow(e)^2)
lb_c <- ATE_tmle_c - 1.96 * se_tmle_c
ub_c <- ATE_tmle_c + 1.96 * se_tmle_c

gc <- PSc * pc1 * e[[A]] + (1 - PSc) * pc0 * (1 - e[[A]])

cat(sprintf("\nTMLE with censoring: %.3f (%.3f; %.3f), SE %.3f, n = %d (%d observed)\n",
            ATE_tmle_c, lb_c, ub_c, se_tmle_c, N_ELIG, sum(obs)))

TC <- "TMLE with censoring (manual)"
savest("T2", TC, "ATE", ATE_tmle_c, se_tmle_c, lb_c, ub_c, NA, NA, n = N_ELIG)
savediag("T2", TC, "N censored",          sum(!obs))
savediag("T2", TC, "Mean Y1",             mean(Y1c_t))
savediag("T2", TC, "Mean Y0",             mean(Y0c_t))
savediag("T2", TC, "Min g x P(observed)", min(gc))
savediag("T2", TC, "SD influence curve",  sd(EIFc))


#########################################################################
############## 8.4.1 - SENSITIVITY: TMLE with censoring  ################
##############          (R ltmle package)                ################
##############                         -> sheet T2       ################
#########################################################################

#===============================================================================
# The censoring node C comes between treatment and outcome. It must be a factor
# ("censored" / "uncensored"), built with BinaryToCensoring(); ltmle then fits
# the outcome model in the uncensored and adds the probability of being
# observed to the treatment mechanism. gform has one formula per A or C node,
# in causal order. As in section 6.1, ltmle rescales BASFI to [0, 1] and
# fluctuates on the logit scale, so it can differ from the manual code.
# The manual code is the one that reproduces Stata exactly.
#===============================================================================

tmle_c_data <- e[, c(W, A, Y)]
tmle_c_data$C <- BinaryToCensoring(is.uncensored = e$uncensored == 1)
tmle_c_data[[Y]][e$uncensored == 0] <- NA   # BASFI of the 15 without ASDAS is not used
tmle_c_data   <- tmle_c_data[, c(W, A, "C", Y)]   # causal order: W, A, C, Y

Qform_c <- setNames(paste("Q.kplus1 ~", paste(c(W, A), collapse = " + ")), Y)
gform_c <- c(paste(A,   "~", paste(W, collapse = " + ")),
             paste("C", "~", paste(c(A, Wc), collapse = " + ")))

set.seed(SEED)
ATElc <- ltmle(tmle_c_data, Anodes = A, Cnodes = "C", Lnodes = NULL, Ynodes = Y,
               Qform = Qform_c, gform = gform_c,
               abar = list(1, 0), estimate.time = FALSE)
slc <- summary(ATElc)
print(slc)

elc <- slc$effect.measures$ATE
savest("T2", "TMLE with censoring (ltmle)", "ATE", elc$estimate, elc$std.dev,
       elc$CI[1], elc$CI[2], NA, NA, n = N_ELIG)

save_sheets()


#########################################################################
####################### Rebuild the excel file  #########################
#########################################################################

# Every section above already rebuilds the whole file, so this is only needed
# after deleting it, after changing DEC or DECD, or after editing the csv files
# by hand. It never re-estimates anything: it just re-reads the two csv files.

save_sheets()


#<<<<<<<<<<<########################################################>>>>>>>>>>>#
#<<<<<<<<<<<##########                                      ########>>>>>>>>>>>#
#<<<<<<<<<<<##########   FIGURES FOR THE MAIN MANUSCRIPT    ########>>>>>>>>>>>#
#<<<<<<<<<<<##########   (code for figures only)            ########>>>>>>>>>>>#
#<<<<<<<<<<<##########                                      ########>>>>>>>>>>>#
#<<<<<<<<<<<########################################################>>>>>>>>>>>#

#===============================================================================
# Nothing below estimates anything. This section only draws the figures of the
# main manuscript from the estimates that Stata wrote to Tables/estimates_stata.csv
# (Analysis.do). The Stata estimates are the main results of the paper; the R
# estimates above are for double checking only and are NOT used here.
#
# It uses only library(ggplot2) and TABLES and FIGURES from the Settings at the
# top, and reads no object estimated in the sections above.
#
# Figure 2  Panel A: total effect decomposed into the average indirect effect
#                    (through ASDAS) and the average direct effect.
#                    Box S2, manual g-formula ignoring MOC -> Table S3
#           Panel B: total effect decomposed into the direct effect and the
#                    path-specific effects through CRP and through ASDAS-PRO.
#                    Box S2.1, manual edge g-formula, main paths without
#                    interaction -> Table S4
#
# As in save_sheets(), the csv file is appended at every run, so the last row
# written for each table/estimator/effect is the current estimate.
#===============================================================================

#####>>>>>> Read the Stata estimates and keep the last row of each estimate

st <- read.csv(paste0(TABLES, "estimates_stata.csv"), stringsAsFactors = FALSE)
st <- st[!duplicated(st[c("table", "estimator", "effect")], fromLast = TRUE), ]

#####>>>>>> Panel A: one mediator (ASDAS), Box S2, Table S3

EST_A <- "G-formula ignoring MOC (manual)"

a_total    <- st[st$table == "S3" & st$estimator == EST_A & st$effect == "ATE", ]
a_direct   <- st[st$table == "S3" & st$estimator == EST_A & st$effect == "ADE", ]
a_indirect <- st[st$table == "S3" & st$estimator == EST_A & st$effect == "AIE", ]

stopifnot(nrow(a_total) == 1, nrow(a_direct) == 1, nrow(a_indirect) == 1)

a_total$label    <- "Total effect"
a_direct$label   <- "Direct effect"
a_indirect$label <- "Indirect effect via ASDAS"

a_total$type    <- "Total"
a_direct$type   <- "Direct"
a_indirect$type <- "Indirect"

a_total$order    <- 1
a_direct$order   <- 2
a_indirect$order <- 3

panA <- rbind(a_total, a_direct, a_indirect)
panA$panel <- "A. One mediator: ASDAS"

#####>>>>>> Panel B: two ordered mediators (CRP -> ASDAS-PRO), Box S2.1, Table S4
#####>>>>>> Main paths, without interaction. For the with-interaction
#####>>>>>> specification, replace "without" by "with" in EST_B.

EST_B <- "Edge g-formula: main paths without interaction (manual)"

b_total  <- st[st$table == "S4" & st$estimator == EST_B & st$effect == "Total", ]
b_direct <- st[st$table == "S4" & st$estimator == EST_B & st$effect == "Direct", ]
b_crp    <- st[st$table == "S4" & st$estimator == EST_B & st$effect == "PSE_M1", ]
b_pro    <- st[st$table == "S4" & st$estimator == EST_B & st$effect == "PSE_M2", ]

stopifnot(nrow(b_total) == 1, nrow(b_direct) == 1, nrow(b_crp) == 1, nrow(b_pro) == 1)

b_total$label  <- "Total effect"
b_direct$label <- "Direct effect"
b_crp$label    <- "Path-specific effect via CRP"
b_pro$label    <- "Path-specific effect via ASDAS-PRO"

b_total$type  <- "Total"
b_direct$type <- "Direct"
b_crp$type    <- "Indirect"
b_pro$type    <- "Indirect"

b_total$order  <- 1
b_direct$order <- 2
b_crp$order    <- 3
b_pro$order    <- 4

panB <- rbind(b_total, b_direct, b_crp, b_pro)
panB$panel <- "B. Two ordered mediators: CRP, then ASDAS-PRO"

#####>>>>>> One data frame for both panels

fig2 <- rbind(panA, panB)

# Text shown on the right: estimate (95% CI), two decimals as in the tables.
# -0.00 is shown as -0.00, as in Supplementary Table S5.
fig2$text <- sprintf("%.2f (%.2f; %.2f)", fig2$estimate, fig2$lci, fig2$uci)

# Rows from top to bottom in the order above. Labels are made unique per panel
# (the panel name is attached and then removed from the axis by labeller).
fig2$row <- paste(fig2$panel, fig2$label, sep = "__")
fig2$row <- factor(fig2$row, levels = rev(fig2$row[order(fig2$panel, fig2$order)]))

fig2$panel <- factor(fig2$panel, levels = c("A. One mediator: ASDAS",
                                            "B. Two ordered mediators: CRP, then ASDAS-PRO"))
fig2$type  <- factor(fig2$type, levels = c("Total", "Direct", "Indirect"))

#####>>>>>> Draw
#####>>>>>> Colours as in Figure 1: green = direct effect, blue = indirect
#####>>>>>> (mediation) pathway; black = total effect.

COL <- c("Total" = "black", "Direct" = "#1B5E3F", "Indirect" = "#2E5283")

XMIN  <- -1.2      # x-axis range (BASFI units)
XMAX  <-  0.6
XTEXT <-  0.75     # where the estimate (95% CI) column starts
XTITL <- -2.40     # where the panel titles start (left of the row labels)

p2 <- ggplot(fig2, aes(x = estimate, y = row, colour = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lci, xmax = uci), width = 0.18, linewidth = 0.6,
                orientation = "y") +
  geom_point(aes(shape = type), size = 3) +
  geom_text(aes(x = XTEXT, label = text), hjust = 0, colour = "black", size = 3.6) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_x_continuous(breaks = seq(-1.2, 0.6, by = 0.2),
                     labels = sprintf("%.1f", seq(-1.2, 0.6, by = 0.2))) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_shape_manual(values = c("Total" = 18, "Direct" = 16, "Indirect" = 16), guide = "none") +
  coord_cartesian(xlim = c(XMIN, XMAX), clip = "off") +
  labs(x = NULL, y = NULL) +     # no axis title: the figure title goes in the manuscript
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor    = element_blank(),
        strip.text.y        = element_blank(),
        strip.placement     = "outside",
        axis.text.y         = element_text(colour = "black", size = 10.5),
        panel.spacing.y     = unit(2.4, "lines"),
        plot.margin         = margin(t = 30, r = 150, b = 8, l = 8))

# Panel titles and the column header, written above the top row of each
# panel (the total effect), instead of in the facet strip

titles <- data.frame(panel = factor(levels(fig2$panel), levels = levels(fig2$panel)),
                     row   = c(paste("A. One mediator: ASDAS", "Total effect", sep = "__"),
                               paste("B. Two ordered mediators: CRP, then ASDAS-PRO",
                                     "Total effect", sep = "__")))

p2 <- p2 +
  geom_text(data = titles, aes(x = XTITL, y = row, label = panel), inherit.aes = FALSE,
            hjust = 0, vjust = -2.8, fontface = "bold", size = 3.9) +
  geom_text(data = titles[1, ], aes(x = XTEXT, y = row), label = "Estimate (95% CI)",
            inherit.aes = FALSE, hjust = 0, vjust = -2.8, fontface = "bold", size = 3.6)

print(p2)

ggsave(paste0(FIGURES, "Figure2_mediation.png"), p2, width = 9, height = 5.5, dpi = 600, bg = "white")
ggsave(paste0(FIGURES, "Figure2_mediation.pdf"), p2, width = 9, height = 5.5, bg = "white")



#===============================================================================
# Figure 3  Total effect of bDMARDs on BASFI at 6 months across the g-methods.
#           Panel A: main analysis, patients with complete data (n = 419):
#                    g-formula (Box S1/S2 -> Table S3), MSM with IPTW (Box S3 ->
#                    Table S6) and TMLE (Box S4 -> Table S7).
#           Panel B: sensitivity analysis for missing 6-month data, the 481
#                    eligible patients with complete baseline data: MSM with IPTW
#                    and censoring weights (run 8.3) and TMLE with censoring
#                    (run 8.4) -> sheet T2.
#
# Reads its own copy of the Stata estimates, so it runs without the Figure 2 code.
#===============================================================================

#####>>>>>> Read the Stata estimates and keep the last row of each estimate

st <- read.csv(paste0(TABLES, "estimates_stata.csv"), stringsAsFactors = FALSE)
st <- st[!duplicated(st[c("table", "estimator", "effect")], fromLast = TRUE), ]

#####>>>>>> Panel A: main analysis (n = 419)

f3_gf   <- st[st$table == "S3" & st$estimator == "G-formula ignoring MOC (manual)" & st$effect == "ATE", ]
f3_msm  <- st[st$table == "S6" & st$estimator == "MSM with sIPTW (manual)"         & st$effect == "ATE", ]
f3_tmle <- st[st$table == "S7" & st$estimator == "TMLE (manual)"                   & st$effect == "ATE", ]

stopifnot(nrow(f3_gf) == 1, nrow(f3_msm) == 1, nrow(f3_tmle) == 1)

f3_gf$label   <- "G-formula"
f3_msm$label  <- "MSM with IPTW"
f3_tmle$label <- "TMLE"

f3_gf$order   <- 1
f3_msm$order  <- 2
f3_tmle$order <- 3

f3A <- rbind(f3_gf, f3_msm, f3_tmle)
f3A$panel <- "A. Main analysis: complete data (n = 419)"

#####>>>>>> Panel B: censoring weights (n = 481)

f3_msmc  <- st[st$table == "T2" & st$estimator == "MSM with sIPTW and sIPCW (manual)" & st$effect == "ATE", ]
f3_tmlec <- st[st$table == "T2" & st$estimator == "TMLE with censoring (manual)"      & st$effect == "ATE", ]

stopifnot(nrow(f3_msmc) == 1, nrow(f3_tmlec) == 1)

f3_msmc$label  <- "MSM with IPTW and censoring weights"
f3_tmlec$label <- "TMLE with censoring"

f3_msmc$order  <- 1
f3_tmlec$order <- 2

f3B <- rbind(f3_msmc, f3_tmlec)
f3B$panel <- "B. Censoring weights: all eligible patients (n = 481)"

#####>>>>>> One data frame for both panels, with a fixed vertical position per row
#####>>>>>> (top to bottom: panel A rows, then panel B rows)

f3A$y <- c(7, 6, 5)
f3B$y <- c(2.5, 1.5)

fig3 <- rbind(f3A, f3B)

fig3$text <- sprintf("%.2f (%.2f; %.2f)", fig3$estimate, fig3$lci, fig3$uci)

#####>>>>>> Layout, in BASFI units on the x axis
#####>>>>>> Everything (panel titles, row labels, estimates) is drawn inside the
#####>>>>>> plotting area, so nothing depends on the width of the system font
#####>>>>>> and nothing can be cut at the edge of the image.

XLEFT3  <- -3.45    # left edge of the figure: panel titles start here
XLAB3   <- -1.50    # row labels end here (right-aligned)
XMIN3   <- -1.40    # effect axis
XMAX3   <-  0.20
XTEXT3  <-  0.35    # estimate (95% CI) column starts here
XRIGHT3 <-  1.45    # right edge of the figure

YTOP_A <- 7.4; YBOT_A <- 4.6   # vertical extent of the grid lines, panel A
YTOP_B <- 2.9; YBOT_B <- 1.1   # vertical extent of the grid lines, panel B
YTIT_A <- 7.9                  # panel title and column header, panel A
YTIT_B <- 3.4                  # panel title, panel B

BREAKS3 <- seq(XMIN3, XMAX3, by = 0.2)

grid3 <- rbind(data.frame(x = BREAKS3, y0 = YBOT_A, y1 = YTOP_A),
               data.frame(x = BREAKS3, y0 = YBOT_B, y1 = YTOP_B))

zero3 <- data.frame(x = 0, y0 = c(YBOT_A, YBOT_B), y1 = c(YTOP_A, YTOP_B))

titles3 <- data.frame(x = XLEFT3, y = c(YTIT_A, YTIT_B),
                      label = c("A. Main analysis: complete data (n = 419)",
                                "B. Censoring weights: all eligible patients (n = 481)"))

#####>>>>>> Draw (black: total effect, as in Figure 2)

p3 <- ggplot(fig3, aes(x = estimate, y = y)) +
  geom_segment(data = grid3, aes(x = x, xend = x, y = y0, yend = y1),
               inherit.aes = FALSE, colour = "grey92", linewidth = 0.5) +
  geom_segment(data = zero3, aes(x = x, xend = x, y = y0, yend = y1),
               inherit.aes = FALSE, colour = "grey55", linetype = "dashed") +
  geom_errorbar(aes(xmin = lci, xmax = uci), width = 0.25, linewidth = 0.6,
                orientation = "y", colour = "black") +
  geom_point(shape = 18, size = 3.5, colour = "black") +
  geom_text(aes(x = XLAB3, label = label), hjust = 1, colour = "black", size = 3.7) +
  geom_text(aes(x = XTEXT3, label = text), hjust = 0, colour = "black", size = 3.6) +
  geom_text(data = titles3, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, fontface = "bold", size = 3.9) +
  annotate("text", x = XTEXT3, y = YTIT_A, label = "Estimate (95% CI)",
           hjust = 0, fontface = "bold", size = 3.6) +
  scale_x_continuous(breaks = BREAKS3, labels = sprintf("%.1f", BREAKS3)) +
  coord_cartesian(xlim = c(XLEFT3, XRIGHT3), ylim = c(0.9, 8.2), expand = FALSE) +
  labs(x = NULL, y = NULL) +     # no axis title: the figure title goes in the manuscript
  theme_minimal(base_size = 11) +
  theme(panel.grid   = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_text(colour = "grey30", size = 10),
        plot.margin  = margin(t = 8, r = 8, b = 8, l = 8))

print(p3)

ggsave(paste0(FIGURES, "Figure3_gmethods.png"), p3, width = 9, height = 4, dpi = 600, bg = "white")
ggsave(paste0(FIGURES, "Figure3_gmethods.pdf"), p3, width = 9, height = 4, bg = "white")

#<<<<<<<<<<<########################################################>>>>>>>>>>>#
#<<<<<<<<<<<##########   END OF THE FIGURES SECTION         ########>>>>>>>>>>>#
#<<<<<<<<<<<########################################################>>>>>>>>>>>#
