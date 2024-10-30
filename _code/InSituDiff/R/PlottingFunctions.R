
#' Choose colors based on a vector of non-negative, right-tailed values
#' 
#' Colors from 0-high, truncating values at a high quantile of the data so more 
#'  than just the most extreme outliers are visible. Use colorZeroCentered for vectors centered around 0.
#' @param x Vector of non-negative values
#' @param maxquant x will be thresholded above at this quantile
#' @param palette Vector of colors, which will be fed into colorRampPalette
#' @importFrom grDevices colorRampPalette
#' @export
colorRightTailed <- function(x, maxquant = 0.995, palette = c("grey80", "darkblue")) {
  if (min(x) < 0) {
    warning("negative values are present in x; these will be set to 0. Did you want to use colorZeroCentered instead?")
  }
  # define colors:
  cols <- grDevices::colorRampPalette(palette)(101)[
    1 + round(100 * pmax(pmin(x / quantile(x, maxquant), 1), 0))]
  return(cols)
  ## create legend:
  #legvals <- seq(0, quantile(x, maxquant), length.out = length(palette))
  #leg <- palette
  #names(leg) <- signif(legvals,2)
  ## output:
  #out <- list(cols = cols, legend = leg)
  #return(out)
}


#' Choose colors based on a vector of non-negative, right-tailed values
#' 
#' Colors from a vector centered around a meaningful 0, truncating values at high/low quantiles of the data so more 
#'  than just the most extreme outliers are visible
#' @param x Vector of non-negative values
#' @param maxquant x will be thresholded above at this quantile
#' @param palette Vector of colors, which will be fed into colorRampPalette
#' @importFrom grDevices colorRampPalette
#' @export
colorZeroCentered <- function(x, maxquant = 0.995, palette = c("darkblue", "blue", "grey80", "red", "darkred")) {
  maxquantile <- max(abs(quantile(x, c(maxquant, 1 - maxquant))))
  cols <- grDevices::colorRampPalette(palette)(101)[
    round(51 + 50 * pmax(pmin(x / maxquantile, 1), -1))]
  return(cols)
  ## create legend:
  #legvals <- seq(-maxquantile, maxquantile, length.out = length(palette))
  #leg <- palette
  #names(leg) <- signif(legvals,2)
  ## output:
  #out <- list(cols = cols, legend = leg)
  #return(out)
}