#' Calculate the wavelet transform coefficient.
#' @param series like plot information or time series.
#' @param mat samples corresponding to plot or time.
#' @return
#' @export
#' x, plot information; y, scale information; z, coefficient after wavelet transform.
#' @examples
#'specwt(series = 1:100, mat = matrix(runif(100*6), nrow=100))
specwt <- function(series, mat){
  midp <- function(x)  {    (x[-1] + x[-length(x)]) / 2  }
  wave <- function(lag)  {    exp(-lag ^ 2 / 2 + 2i * pi * lag) / pi ^ 4  } # Morlet
  scales = 2 ^ midp(seq(log2(2 * median(diff(series))), log2(diff(range(series))/2),
                        length = length(series) + 1))
  lmat = outer(seq(min(series), max(series), length = length(series)), series, "-")
  w = foreach(s = scales) %dopar%  {
    Conj(wave(lmat / s) / s ^ 0.5) %*% as.matrix(mat)
  }
  w = array(unlist(w), dim = c(length(series), ncol(mat), length(scales)))
  w = aperm(w, perm = c(1, 3, 2))
  return(list(x = seq(min(series), max(series), length = length(series)), y = scales, z = w))
}
