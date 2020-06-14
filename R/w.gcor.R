#' Calculate the global wavelet correlation.
#' @param data the dataset of species and locations.
#' @param w the result of species dataset calculate by 'specwt'.
#' @param scales the scale selected to analysis.
#' @param reps times the bootstrap were performed.
#' @return
#' @export
#' @examples
#' w.gcor(data = Fauna.data[,-12], w = specwt(1:119, Fauna.data[,-12]), scales = 57, reps = 500)
w.gcor <- function (data, w, scales, reps) {
  require(foreach, warn.conflicts = F)
  null.model = function(w, scales, reps) {
    boot = foreach(i = 1:reps, .combine=c, .inorder = FALSE,
                   .export=c('reps','scales'), .packages='foreach') %dopar% {
                     rphase = t(array(runif(2, -pi, pi), dim = c(2, length(w$x))))
                     zp = w$z[ , scales, , drop = FALSE] * complex(argument = rphase)
                     2*(colSums(abs(apply(zp, c(1, 2), sum)))/colSums(apply(abs(zp), c(1, 2), sum)))-1
                   }
    res = c(NA,3)
    res <- c(mean(boot) + qnorm(0.025,mean=0,sd=0.5) * sd(boot),
             (2*(colSums(abs(apply(w$z, c(1, 2), sum)))/
                   colSums(apply(abs(w$z), c(1, 2), sum)))-1)[scales],
             mean(boot) + qnorm(0.975,mean=0,sd=0.5) * sd(boot) )
    return(res)
  }
  mat <- foreach(i = 1:dim(w$z)[3], .combine = rbind, .packages='doParallel') %dopar% {
    foreach(j = 1:dim(w$z)[3], .combine = rbind) %dopar% {
      null.model(list(x=w$x,y=w$y,z=w$z[,,c(i,j)]),scales, reps=reps)
    }
  }
  result <- list(lower = matrix(mat[,1],dim(w$z)[3]),
                 res = matrix(mat[,2],dim(w$z)[3]),
                 upper = matrix(mat[,3],dim(w$z)[3]) )
  colnames(result[[1]]) = rownames(result[[1]]) = colnames(result[[2]]) = colnames(data)
  rownames(result[[2]]) = colnames(result[[3]]) = rownames(result[[3]]) = colnames(data)
  diag(result[[1]]) = diag(result[[2]]) = diag(result[[3]]) = NA
  return(result)
}
