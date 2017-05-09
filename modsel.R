selectModelIc = function(HmmClass, xx, m.max, do.plot=FALSE, ...) {
  aic <- list() # TODO: should be a numeric()
  bic <- list() # TODO: should be a numeric()
  for (m in 1:m.max) {
    cat('m =', m, '\n')
    hmm <- HmmClass(m, xx, ...)
    hmm$fit(xx)
    aic[m] <- hmm$aic(xx)
    bic[m] <- hmm$bic(xx)
  }
  if (do.plot) {
    matplot(cbind(aic, bic), type='b', pch=c(1,2),
            main='', xlab='m', ylab='')
  }
  legend("topleft", bty = "n",
         legend = c('AIC', 'BIC'),
         col = c('black', 'red'),
         lty = c('solid', 'dashed'))
  c(which.min(aic), which.min(bic))
}
