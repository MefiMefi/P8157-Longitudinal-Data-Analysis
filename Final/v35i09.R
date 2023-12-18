##################################################################
# Title: R script file for JSS paper "JM: An R Package for the   #
#        Joint Modelling of Longitudinal and Time-to-Event Data" #
# Author: Dimitris Rizopoulos                                    #
##################################################################


# load packages 'JM' and 'lattice'
library("JM")
library("lattice")


#####################
# Descriptive Plots #
#####################

# longitudinal outcome
xyplot(sqrt(CD4) ~ obstime | drug, group = patient, data = aids, 
    xlab = "Months", ylab = expression(sqrt("CD4")), col = 1, type = "l")

# survival outcome
plot(survfit(Surv(Time, death) ~ drug, data = aids.id), conf.int = FALSE, 
    mark.time = TRUE, col = c("black", "red"), lty = 1:2, 
    ylab = "Survival", xlab = "Months") 
legend("topright", c("ddC", "ddI"), lty = 1:2, col = c("black", "red"), 
    bty = "n")


###################
# Naive Cox Model #
###################

td.Cox <- coxph(Surv(start, stop, event) ~ drug + sqrt(CD4), data = aids)
summary(td.Cox)


###############
# Joint Model #
###############

# first we fit a linear mixed effects model with lme(),
fitLME <- lme(sqrt(CD4) ~ obstime + obstime:drug, random = ~ obstime | patient, data = aids)
# and a Cox model with coxph() -- you need 'x = TRUE' in the call to coxph()
fitSURV <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# fit the joint model using function jointModel()
fit.JM <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "piecewise-PH-GH")
summary(fit.JM)

# likelihood ratio test for no treatment effect in the survival process:
# we fit the joint model under the null hypothesis
fitSURV2 <- coxph(Surv(Time, death) ~ 1, data = aids.id, x = TRUE)
fit.JM2 <- jointModel(fitLME, fitSURV2, timeVar = "obstime", method = "piecewise-PH-GH")
# the likelihood ratio test is performed with the anova() method
anova(fit.JM2, fit.JM)


###################
# Residuals Plots #
###################

# the standard call to the plot() method for objects of class 'jointModel'
par(mfrow = c(2, 2))
plot(fit.JM)

# a useful function used in the residual plots below
plotResid <- function (x, y, ...) {
    plot(x, y, ...)
    lines(lowess(x, y), col = "red", lwd = 2)
    abline(h = 0, lty = 3, col = "grey", lwd = 2)
}

par(mfrow = c(2, 2))

# Subject-Specific Residuals vs Fitted Values
resSubY <- residuals(fit.JM, process = "Longitudinal", type = "stand-Subject")
fitSubY <- fitted(fit.JM, process = "Longitudinal", type = "Subject")
plotResid(fitSubY, resSubY, xlab = "Fitted Values", ylab = "Residuals",
    main = "Subject-Specific Residuals vs Fitted Values")

# Marginal Residuals vs Fitted Values
resMargY <- residuals(fit.JM, process = "Longitudinal", type = "stand-Marginal")
fitMargY <- fitted(fit.JM, process = "Longitudinal", type = "Marginal")
plotResid(fitMargY, resMargY, xlab = "Fitted Values", ylab = "Residuals",
    main = "Marginal Residuals vs Fitted Values")

# Martingale Residuals vs Fitted Values
resMartT <- residuals(fit.JM, process = "Event", type = "Martingale")
fitSubY <- fitted(fit.JM, process = "Longitudinal", type = "EventTime")
plotResid(fitSubY, resMartT, xlab = "Fitted Values", ylab = "Residuals",
    main = "Martingale Residuals vs Fitted Values")

# Cox-Snell Residuals
resCST <- residuals(fit.JM, process = "Event", type = "CoxSnell")
sfit <- survfit(Surv(resCST, death) ~ 1, data = aids.id)
plot(sfit, mark.time = FALSE, conf.int = TRUE, lty = 1:2, 
    xlab = "Cox-Snell Residuals", ylab = "Survival Probability", 
    main = "Survival Function of Cox-Snell Residuals")
curve(exp(-x), from = 0, to = max(aids.id$Time), add = TRUE, col = "red", lwd = 2)


#################################
# Multiple-Imputation Residuals #
#################################

set.seed(123) # we set the seed for reproducibility
# to calculate the MI residuals, set logical argument 'MI' to TRUE
res.MI <- residuals(fit.JM, process = "Longitudinal", type = "stand-Marginal", MI = TRUE)

# Extract components from the returned list:
# fitted values corresponding to missing Y's
fitMargY.miss <- res.MI$fitted.valsM
# multiply imputed residuals corresponding to missing Y's
resMargY.miss <- res.MI$resid.valsM

# Put together the multiply imputed residuals and the observed residuals, and
# the fitted values that correspond to the observed and imputed y_i^m, respectively
M <- ncol(resMargY.miss) # number of imputations
resMargY.MI <- c(resMargY, resMargY.miss)
fitMargY.MI <- c(fitMargY, rep(fitMargY.miss, M))

# data frame with the above defined objects and corresponding case weights
dat.resid <- data.frame(
    resid = resMargY.MI,
    fitted = fitMargY.MI,
    weight = c(rep(1, length(resMargY)), rep(1/M, length(resMargY.miss)))
)
fit.loess <- loess(resid ~ fitted, data = dat.resid, weights = weight)
nd <- data.frame(fitted = seq(min(fitMargY.MI), max(fitMargY.MI), 
        length.out = 100))
prd.loess <- predict(fit.loess, nd)


# Plot of the multiply-imputed standardized marginal residuals
plot(range(fitMargY.MI), range(resMargY.MI), type = "n",
    xlab = "Fitted Values", ylab = "MI Standardized Marginal Residuals")
abline(h = 0, lty = 2)
points(rep(fitMargY.miss, M), resMargY.miss, cex = 0.5, col = "grey")
points(fitMargY, resMargY)
lines(lowess(fitMargY, resMargY), lwd = 2)
lines(nd$fit, prd.loess, col = "grey", lwd = 2)


########################################
# Prediction of Survival Probabilities #
########################################

set.seed(123) # we set the seed for reproducibility
# a data frame with 4 patients who have not died by the end of the study
ND <- aids[aids$patient %in% c("7", "15", "117", "303"), ]
# calculate predicted survival probabilities
predSurv <- survfitJM(fit.JM, newdata = ND, idVar = "patient", last.time = "Time")
predSurv

# plot of the predicted survival probabilities
# for the 4 patients
par(mfrow = c(2, 2))
plot(predSurv, conf.int = TRUE)

# plot of predicted survival probabilities for Patient '7'
# in different scales
par(mfrow = c(2, 2))
plot(predSurv, which = "7", conf.int = TRUE)
plot(predSurv, which = "7", conf.int = TRUE, fun = log, ylab = "log Survival")
plot(predSurv, which = "7", conf.int = TRUE, fun = function (x) -log(x), 
        ylab = "Cumulative Risk")


# plot of the predicted survival probabilities for Patient '7', including
# his observed square root CD4 cell count values
plot(predSurv, estimator = "median", which = "7", conf.int = TRUE, include.y = TRUE)
