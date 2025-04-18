#' Goodness of fit.
#'
#' @param nn sample size.
#' @param scl scale estimate.
#' @param val minimization value.
#' @param spl number of splines1.
#' @param freq number of splines2.
#' @param criterion goodness of fit criterion.
#'
#' @return Goodness of fit for the selectec criterion.
#' @export
goodness <- function(nn, scl, val, spl, freq, criterion) {
  switch(criterion,
    hic = log(scl^2 * val / nn) + spl * log(nn) / (2 * nn) + 2 * freq / nn,
    hic2 = log(scl^2 * val / nn) + freq * log(nn) / (2 * nn) + 2 * spl / nn,
    hicGR2 = log(scl^2) + val / nn + freq * log(nn) / (2 * nn) + 2 * spl / nn,
    hicGR = log(scl^2) + val / nn + spl * log(nn) / (2 * nn) + 2 * freq / nn,
    bic = log(scl^2 * val / nn) + (spl + freq) * log(nn) / (2 * nn),
    #bic1 = log(scl^2 * val / nn) + (spl + freq) * log(nn) / (nn),
    bic1 = log(scl^2 * val / nn) + (freq) * log(nn) / (nn),
    bicGR = log(scl^2) + val / nn + (spl + freq) * log(nn) / (2 * nn),
    aicCL = log(scl^2 * val / nn) + 2 * (spl + freq) / nn,
    aicCLGR = log(scl^2) + val / nn + 2 * (spl + freq) / nn
  )
}
