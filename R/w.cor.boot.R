#' Calculate the local wavelet correlation with CI.
#'
#' @param data the dataset of species and locations.
#' @param w the result of species dataset calculate by 'specwt'.
#' @param scales the scale selected to analysis.
#' @param reps times the bootstrap were performed.
#' @return
#' @export
#' @examples
#' w.cor.boot(data = Fauna.data[,-12], w = specwt(1:119, Fauna.data[,-12]), scales = 57, reps = 500)
w.cor.boot <- function(data, w, scales, reps) {
  wmr.cor <- function(wmat, smoothing = 1){
    with(wmat, {
      lmat = outer(seq(min(x), max(x), length = length(x)), x, "-")
      Gauss <- function(lag) {  exp(-lag ^ 2 / 2) / sqrt(2 * pi)  }
      kern = Gauss(lmat / y[scales] / smoothing)
      modsv = kern %*% Mod(rowSums(z))
      smodv = kern %*% rowSums(Mod(z))
      modrat = 2*(modsv / smodv)-1
      modrat = as.vector(modrat)
      return(modrat)
    })
  }
  mat <- foreach(i = 1:dim(w$z)[3], .packages = "doParallel") %dopar% {
    foreach(j = 1:dim(w$z)[3], .combine = cbind) %dopar% {
      (wmr.cor(list(x=w$x,y=w$y,z=w$z[,scales,c(i,j)]),1))
    }
  }
  z.boot <- foreach(i = 1:dim(w$z)[3], .packages = "doParallel") %dopar% {
    foreach(j = 1:dim(w$z)[3], .combine = cbind) %dopar% {
      mr.boot = foreach(t = 1:reps, .combine=cbind, .inorder = FALSE,
                        .export=c('reps','scales'), .packages='foreach') %dopar% {
                          rphase = t(array(runif(2, -pi, pi), dim = c(2, length(w$x))))
                          zp = w$z[ , scales, c(i,j)] * complex(argument = rphase)
                          (wmr.cor(list(x=w$x,y=w$y,z=zp),1))
                        }
      res = foreach(k = 1:nrow(mr.boot),.combine = c) %dopar% {
        ecdf(mr.boot[k, ])(mat[[i]][k,j])
      }
      res
    }
  }
  mat.res <- foreach(i = 1:dim(w$z)[3] ) %dopar% {
    ifelse(z.boot[[i]]>=0.95 | z.boot[[i]]<=0.05, mat[[i]], 0 )
  }
  result <- foreach(i = 1:dim(w$z)[1])%dopar%{
    res <- foreach(j = 1:dim(w$z)[3], .combine = cbind)%dopar%{
      mat.res[[j]][i,]
    }
    colnames(res) <- rownames(res) <- colnames(data)
    diag(res) <- 0
    res
  }
  for(i in 1:dim(w$z)[1]){names(result)[i] <- i }
  return(result)
}
