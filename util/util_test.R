testIntepolation <- function(hmm, xx, x.probes, t.probes) {
  res <- matrix(0, length(xx), 7)
  colnames(res) <- c('real', 'mode', 'mean', 'median', 'd.mode', 'd.mean', 'd.median')
  t.probes <- 1:length(xx)
  res[, 1] <- xx
  avgs <- c('mode', 'mean', 'median')
  for (i in 1:length(avgs)) {
    res[, i+1] <- hmm$interpolate(xx, x.probes, t.probes, avgs[i])
    res[, 4+i] <- abs(res[, 1+i] - res[, 1])
  }
  return(res)
}