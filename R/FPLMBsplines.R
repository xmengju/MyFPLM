#' Best FPLMBsplines_fit given by a model selection criterion.
#'
#' @description Fit a FPLM model for different spline basis sizes and picks the
#'     best one according to a specified model selection criterion.
#' @param y the vector of scalar responses.
#' @param x a matrix of the functional covariates, where each row contains the
#'     functions evaluated on a (common) grid.
#' @param u the values of the explanatory variable that enters the model
#'     non-parametrically.
#' @param t the grid over which the functional covariates were evaluated.
#' @param range_freq a vector of B-spline basis sizes to try for the functional
#'     regression coefficient.
#' @param range_spl a vector of B-spline basis sizes to try for the
#'     non-parametric component.
#' @param norder the order of the B-Splines.
#' @param fLoss string specifying the loss function. 'ls' for least squares,
#'     'huang' for Huber, 'lmbrob' for MM-estimator.
#' @param criterion criterion for model selection.
#' @param trace a logical argument indicating whether partial results are
#'     printed.
#'     
#' @return A list including the following components:
#' \itemize{
#' \item{fit}{fitted parameters}
#' \item{spl}{chosen number of splines for the non-parametric component}
#' \item{freq}{chosen number of splines for the funcitonal regression coefficient}
#' }
#' @import fda robustbase
#' 
#' @export
FPLMBsplines <- function(y, x, u, t, range_freq = range_default,
                         range_spl = range_default, norder = 4,
                         fLoss = "lmrob", criterion = "bic1",
                         trace = FALSE) {

    print("running the function!")
    ## Some Setup
    opt <- spl_opt <- freq_opt <- fit_opt <- Inf
    n <- length(y)
    range_default <- floor(max(n^(1 / 5), norder)):
        floor(2 * (norder + n^(1 / 5)))

    ## Double loop
    for (spl in range_spl) {
        for (freq in range_freq) {
            fit <- FPLMBsplines_fit(y, x, u, t, freq, spl, norder, fLoss)
            val <- fit$value
            scl <- fit$scale
            crt <- goodness(n, scl, val, spl, freq, criterion)

            if (crt < opt) {
                opt <- crt
                spl_opt <- spl
                freq_opt <- freq
                fit_opt <- fit
            }
            if (trace) print(c("spl" = spl, "freq" = freq, "crit" = crt))
        }
    }
    
    print(c("optimal",   freq_opt))

    ## Best fit
    kns <- seq(min(u), max(u), length = spl_opt - norder + 2)
    base <- create.bspline.basis(
        rangeval = range(u),
        norder = norder,
        breaks = kns
    )
    #spl_uu <- getbasismatrix(u, base)
    #fit_opt$eta_est <- spl_uu %*% fit_opt$spl
    dt <- min(diff(t))
    #fit_opt$fitted <- as.vector(x %*% fit_opt$slope_fun * dt + fit_opt$eta_est)

    fit_opt$fitted <- as.vector(x %*% fit_opt$slope_fun * dt)
    
    return(list(fit = fit_opt, spl = spl_opt, freq = freq_opt, dt = dt))
}

#' @export
FPLMBsplines.predict <- function(model, newx, newy){
    
    fitted <- as.vector(newx %*% model$fit$slope_fun * model$dt)
    tmp <- list(pred= fitted)
    if(!missing(newy)){
        error <- mean( (fitted - newy)^2)
        tmp <- c(tmp, list(error = error))
    }
    return(tmp)

}

