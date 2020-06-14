#' Calculate the local wavelet correlation.
#' @param data the dataset of species and locations.
#' @param w the result of species dataset calculate by 'specwt'.
#' @param scales the scale selected to analysis.
#' @return
#' @export
#' @examples
#' w.lcor(data = Fauna.data[,-12], w = specwt(1:119, Fauna.data[,-12]), scales = 57)
w.lcor <- function(data, w, scales) {   # wavelet local correlation function
  wmr.cor <- function(wmat, smoothing = 1){
    with(wmat, {
      lmat = outer(seq(min(x), max(x), length = length(x)), x, "-")
      Gauss <- function(lag) {  exp(-lag ^ 2 / 2) / sqrt(2 * pi)  }
      modrat = foreach(i = 1:length(y), .combine = cbind,
                       .export = "smoothing", .packages = "mvcwt") %dopar% {
                         kern = Gauss(lmat / y[i] / smoothing)
                         modsv = kern %*% Mod(rowSums(z, dims = 2))[,i]
                         smodv = kern %*% rowSums(Mod(z), dims = 2)[,i]
                         2*(modsv / smodv)-1
                       }
      dim(modrat) = c(length(x), length(y), 1)
      return(list(x = x, y = y, z = modrat))
    })
  }
  mat <- foreach(i = 1:dim(w$z)[3], .packages = "doParallel") %dopar% {
    foreach(j = 1:dim(w$z)[3], .combine = cbind) %dopar% {
      (wmr.cor(list(x=w$x,y=w$y,z=w$z[,,c(i,j)]),1))$z[,,1][,scales]
    }
  }
  cor.res <- foreach(i = 1:dim(w$z)[1]) %dopar% {
    res <- foreach(j = 1:dim(w$z)[3], .combine = cbind) %dopar% {
      mat[[j]][i,]
    }
    colnames(res) <- rownames(res) <- colnames(data)
    res
  }
  for(i in 1:dim(w$z)[1]){names(cor.res)[i] <- i }
  return(cor.res)
}
