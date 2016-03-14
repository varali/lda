glm.nb.fit <- function(x, y, w = NULL, start = NULL, etastart=NULL, mustart=NULL, offset=NULL, method="glm.fit",control = glm.control(...), contrasts = NULL, ..., init.theta, link = log, intercept=TRUE)
{
#Provides functionality of glm.fit() for negative binomial regression models
#  
# based on glm.nb() of MASS 7-3.11
# @author: Thomas Rusch
# @version: 1.4.2011
#  
# @input:
# x,y: ‘x’ is a design matrix of dimension ‘n * p’, and ‘y’ is a vector of observations of length ‘n’. Make sure that x has the intercept (if any) in the first column.
#
# w: an optional vector of ‘prior weights’ to be used in the fitting process.  Should be ‘NULL’ or a numeric vector.
#  
# start: starting values for the parameters in the linear predictor.
#  
# etastart: starting values for the linear predictor.
#
# mustart: starting values for the vector of means.
#  
# offset: this can be used to specify an _a priori_ known component to be included in the linear predictor during fitting.  This should be ‘NULL’ or a numeric vector of length equal to the number of cases.
#  
# method: fitting method used; currently only glm.fit is implemented (see ?glm)
#  
# glm.control: a list of parameters for controlling the fitting process.
#  
# contrasts: an optional list. See the ‘contrasts.arg’ of ‘model.matrix.default’.
# 
# init.theta: starting values for the dispersion parameter (theta).
#  
# link: a specification for the model link function. Default is log link. See ?make.link and ?family for other links.
#  
# intercept: logical. Should an intercept be included in the _null_ model?
#
#  
# @Output:
#
#     ‘glm.nb.fit’ returns an object of class inheriting from ‘"negbin" which inhetits from "glm"’ which  inherits from the class ‘"lm"’. See later in this section.  If a  non-standard ‘method’ is used, the object will also inherit from the class (if any) returned by that function.
#
 #    The function ‘summary’ (i.e., ‘summary.negbin’) can be used to obtain
 #    or print a summary of the results and the function ‘anova’ (i.e.,
 #    ‘anova.negbin’) to produce an analysis of variance table.
#
#     The generic accessor functions ‘coefficients’, ‘effects’,
#     ‘fitted.values’ and ‘residuals’ can be used to extract various
#     useful features of the value returned by ‘glm.nb.fit’.
#
#     ‘weights’ extracts a vector of weights, one for each case in the
#     fit (after subsetting and ‘na.action’).
#
#     An object of class ‘"glm"’ is a list containing at least the
#     following components:
#
# coefficients: a named vector of coefficients
#
# residuals: the _working_ residuals, that is the residuals in the final
#          iteration of the IWLS fit.  Since cases with zero weights are
#          omitted, their working residuals are ‘NA’.
#
# fitted.values: the fitted mean values, obtained by transforming the
#          linear predictors by the inverse of the link function.#
#
#    rank: the numeric rank of the fitted linear model.
#
#  family: the ‘family’ object used.
#
# linear.predictors: the linear fit on link scale.
#
# deviance: up to a constant, minus twice the maximized log-likelihood.
#          Where sensible, the constant is chosen so that a saturated
#          model has deviance zero.
#
#     aic: A version of Akaike _An Information Criterion_, minus twice
#          the maximized log-likelihood plus twice the number of
#          parameters, computed by the ‘aic’ component of the family.
#          For binomial and Poison families the dispersion is fixed at
#          one and the number of parameters is the number of
#          coefficients. For gaussian, Gamma and inverse gaussian
#          families the dispersion is estimated from the residual
#          deviance, and the number of parameters is the number of
#          coefficients plus one.  For a gaussian family the MLE of the
 #         dispersion is used so this is a valid value of AIC, but for
 #         Gamma and inverse gaussian families it is not.  For families
 #         fitted by quasi-likelihood the value is ‘NA’.
#
#null.deviance: The deviance for the null model, comparable with
#          ‘deviance’. The null model will include the offset, and an
#          intercept if there is one in the model.  Note that this will
#          be incorrect if the link function depends on the data other
#          than through the fitted mean: specify a zero offset to force
#          a correct calculation.
#
#    iter: the number of iterations of IWLS used.
#
# weights: the _working_ weights, that is the weights in the final
#          iteration of the IWLS fit.
#
#prior.weights: the weights initially supplied, a vector of ‘1’s if none
#          were.
#
#df.residual: the residual degrees of freedom.
#
# df.null: the residual degrees of freedom for the null model.
#
#       y: if requested (the default) the ‘y’ vector used. (It is a
#          vector even for a binomial model.)
#
#       x: if requested, the model matrix.
#
#
#converged: logical. Was the IWLS algorithm judged to have converged?
#
#boundary: logical. Is the fitted value on the boundary of the
#          attainable values?#
#
#  offset: the offset vector used.
#
# control: the value of the ‘control’ argument used.#
#
#  method: the name of the fitter function used, currently always
#          ‘"glm.fit"’.#
#
# contrasts: (where relevant) the contrasts used.
#
# xlevels: (where relevant) a record of the levels of the factors used
#          in fitting.#
#
# init.theta: starting values used for theta
#
# theta: estimated theta at convergence
#  
# SE.theta: Standard error of theta at convergence
#  
# twologlik: Twice the loglikelihood
  
####FUNCTION START HERE 
#negative binomial log likelihood
    loglik <- function(n, th, mu, y, w)
        sum(w*(lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
               y * log(mu + (y == 0)) - (th + y) * log(th + mu)))
 #initialise family: with init.theta=NULL a poisson is fitted, else a negative.binomial(init.theta)
    X <- x #just because I'm lazy
    Y <- y #just because I'm lazy, so I can use the V&R code right away
    link <- substitute(link)
    fam0 <- if(missing(init.theta))
        do.call("poisson", list(link = link))
    else
        do.call("negative.binomial", list(theta = init.theta, link = link))

    dots <- list(...)
    #if missing weights, all weights equal
    if(!length(w)) w <- rep(1, nrow(X))
    else if(any(w < 0)) stop("negative weights not allowed")
 
    n <- length(Y)
    #check for method and copy glm.fit if method="glm.fit" 
    #probably redundant but reused it from V&R
    if(!missing(method)) {
        if(!exists(method, mode = "function"))
            stop("unimplemented method: ", sQuote(method))
        glm.fitter <- get(method)
    } else {
        method <- "glm.fit"
        glm.fitter <- stats::glm.fit
      }
    #computation starts here
    if(control$trace > 1) message("Initial fit:")
    #glm.fitter is pretty much glm.fit with some specified parameters
    #is a poisson if init.theta=NULL
    fit <- glm.fitter(x = X, y = Y, w = w, start = start,
                      etastart = etastart, mustart = mustart,
                      offset = offset, family = fam0,
                      control = list(maxit=control$maxit,
                      epsilon = control$epsilon,
                      trace = control$trace > 1),
                      intercept = intercept)
    class(fit) <- c("glm", "lm")
    mu <- fit$fitted.values
    #calculate Maximum Likelihood theta
    th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
                             control$trace> 2))
    if(control$trace > 1)
        message("Initial value for theta:", signif(th))
    #use ML theta in family
    fam <- do.call("negative.binomial", list(theta = th, link = link))
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    g <- fam$linkfun
    Lm <- loglik(n, th, mu, y, w)
    Lm0 <- Lm + 2 * d1
    #iterate through until convergence or max.iter is reached 
    while((iter <- iter + 1) <= control$maxit &&
          (abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon) {
        eta <- g(mu)
        #calculate glm.fit
        fit <- glm.fitter(x = X, y = Y, w = w, etastart =
                          eta, offset = offset, family = fam,
                          control = list(maxit=control$maxit,
                          epsilon = control$epsilon,
                          trace = control$trace > 1),
                          intercept = intercept)
        t0 <- th
        #ML theta
        th <- theta.ml(Y, mu, sum(w), w, limit=control$maxit,
                       trace = control$trace > 2)
        fam <- do.call("negative.binomial", list(theta = th, link = link))
        mu <- fit$fitted.values
        del <- t0 - th
        Lm0 <- Lm
        Lm <- loglik(n, th, mu, Y, w)
        if(control$trace) {
           Ls <- loglik(n, th, Y, Y, w) 
            Dev <- 2 * (Ls - Lm)
            message("Theta(", iter, ") =", signif(th),
                    ", 2(Ls - Lm) =", signif(Dev))
        }
    }
    #warnings if not converged
    if(!is.null(attr(th, "warn"))) fit$th.warn <- attr(th, "warn")
    if(iter > control$maxit) {
        warning("alternation limit reached")
        fit$th.warn <- gettext("alternation limit reached")
    }
    
  # If an offset and intercept are present, iterations are needed to
  # compute the Null deviance; these are done here, unless the model
  # is NULL, in which case the computations have been done already
  # NEW
    if(length(offset) && identical(as.numeric(X[,1]),rep(1,n))) {
      #check if offset is there and if there is an intercept in X, i.e. first column is 1 vector
        null.deviance <-
            if(dim(X)[2]) #check if X is a matrix
                glm.fitter(X[,1], Y, w,
                           offset = offset, family = fam,
                           control = list(maxit=control$maxit,
                           epsilon = control$epsilon,
                           trace = control$trace > 1),
                           intercept = TRUE
                           )$deviance #fit an intercept only model
           else fit$deviance
        fit$null.deviance <- null.deviance
    }
    #wrap up
    class(fit) <- c("negbin", "glm", "lm")
    fit$init.theta <- signif(as.vector(th), 10)
    fit$link <- link 
    fit$theta <- as.vector(th)
    fit$x <- X 
    fit$y <- Y 
    fit$SE.theta <- attr(th, "SE") 
    fit$twologlik <- as.vector(2 * Lm)
    fit$aic <- -fit$twologlik + 2*fit$rank + 2
    fit$contrasts <- attr(X, "contrasts")
    fit$method <- method
    fit$control <- control
    fit$offset <- offset
    fit
}

##############################
########S4 negbin statsmodel
##############################

negbinModel <- new("StatModel",
    capabilities = new("StatModelCapabilities"),
    name = "negative binomial generalized linear regression model",
    dpp = ModelEnvFormula,
    fit = function(object, weights = NULL, ...){
        if (is.null(weights)) {
            z <- glm.nb.fit(object@get("designMatrix"),object@get("response")[,1],mustart=NULL,etastart=NULL,control=glm.control(trace=TRUE),intercept=all(object@get("designMatrix")[,1] == 1),...)
        } else {
            z <- glm.nb.fit(object@get("designMatrix"),object@get("response")[,1],
                         w = weights, mustart=NULL,etastart=NULL,control=glm.control(trace=TRUE),intercept=all(object@get("designMatrix")[,1] == 1),...)
        }
        class(z) <- c("negbinModel", "negbin", "glm", "lm")
        z$offset <- 0
        z$contrasts <- attr(object@get("designMatrix"), "contrasts")
        ## terms should be there, but still need to
	## be worked around in predictions
        z$terms <- attr(object@get("input"), "terms")
        z$predict_response <- function(newdata = NULL) {
            if (!is.null(newdata)) {
                penv <- new.env()
                object@set("input", data = newdata, env = penv)
                dm <- get("designMatrix", envir = penv, inherits = FALSE)
            } else {
                dm <- object@get("designMatrix")
            }
            pr <- z$family$linkinv(drop(dm %*% z$coef))
            return(pr)
        }
        z$addargs <- list(...)
	z$ModelEnv <- object
        z
    },
    predict = function(object, newdata = NULL, ...) 
        object$predict_response(newdata = newdata)
)

predict.negbinModel <- function(object, newdata = NULL, ...)
    object$predict_response(newdata = newdata) 

fitted.negbinModel <- function(object, ...)
    object$predict_response()

print.negbinModel <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    fam <- x$family$family
    substr(fam, 1, 1) <- toupper(substr(fam, 1, 1))
    cat(paste(fam, "GLM with coefficients:\n"))
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    invisible(x)
}


summary.negbinModel <- function(object, dispersion = 1, correlation = FALSE, ...)
{
    if(is.null(dispersion)) dispersion <- 1
    summ <- c(summary.glm(object, dispersion = dispersion,
                          correlation = correlation),
              object[c("theta", "SE.theta", "twologlik", "th.warn")])
    class(summ) <- c("summary.negbin", "summary.glm")
    summ
}

model.matrix.negbinModel <- function(object, ...)
  object$ModelEnv@get("designMatrix")

reweight.negbinModel <- function(object, weights, ...) {
    fit <- negbinModel@fit    
    do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}
