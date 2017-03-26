plotForecast <- function(fc, xx, h.max, xx.true=numeric(0)) {
  T <- length(xx)
  ylim <- range(c(xx, xx.true))
  plot(xx, type='o', pch=20, xlab='', ylab='', bty='n',
       xlim=c(1, T+h.max), ylim=ylim)
  tt.fc <- (T+1):(T+h.max)
  band.colors <- c('#B1B5CE', '#DBDBDF')
  for (j in rev(1:ncol(fc$upper))) {
    polygon(c(tt.fc, rev(tt.fc)),
            c(fc$lower[,j], rev(fc$upper[,j])),
            col=band.colors[j], border=FALSE)
  }
  par(new=TRUE)
  plot(tt.fc, fc$forecast, col='blue',
       type='o', pch=20, axes=FALSE, xlab='', ylab='',
       xlim=c(1, T+h.max), ylim=ylim)
  abline(v=T, lty=2)
  if (length(xx.true) != 0) {
    par(new=TRUE)
    plot(tt.fc, xx.true,
         xlim=c(1, T+h.max), ylim=ylim,
         xlab='', ylab='',
         type='o', pch=21, lty='dashed', axes=FALSE)
  }
  box()
  # close-up
  # plot((T.train+1):T, xx[(T.train+1):T], col='green', type='o', pch=20, ylim=ylim)
  # par(new=TRUE)
  # plot((T.train+1):T, x.probes[apply(F, 1, which.max)],
  #      type='o', pch=1, axes=FALSE, ylim=ylim)
  # plotCI((T.train+1):T, x.probes[apply(F, 1, which.max)],
  #        ui=x.probes[ui], li=x.probes[li], ylim=ylim,
  #        add=TRUE, pch=NA)
}

plotForecastDist <- function(fc.dist, xx, x.probes, h.plot=c()) {
  plot.idcs <- 1:nrow(fc.dist)
  if (length(h.plot) != 0) {
    plot.idcs <- h.plot
  }
  n.plots <- length(plot.idcs)
  par(mfrow=c(floor(sqrt(n.plots)), ceiling(n.plots/floor(sqrt(n.plots)))))
  for (h in plot.idcs) {
    plot(x.probes, fc.dist[h,], type="h", main='',
         xlab='Observations', ylab='Probabilities',
         xlim=range(xx), ylim=c(0, max(fc.dist)))
    title(paste("h =", h), line = 2)
    par(new=TRUE)
    argmax = which.max(fc.dist[h,])
    plot(x.probes[argmax], fc.dist[h, argmax], type="h", lwd=3, axes=FALSE,
         xlab='', ylab='',
         xlim=range(xx), ylim=c(0, max(fc.dist)))
    axis(3, at=round(x.probes[argmax], 2))
    # cat('h = ', h, ': F[', argmax, '] = ', fc.dist[h, argmax], '\n', sep='')
  }
  par(mfrow=c(1,1))
}

plotA <- function(hmm) {
  library('corrplot')
  col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                             "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                             "#4393C3", "#2166AC", "#053061"))
  corrplot(hmm$A,
           method='circle', is.corr=FALSE, addCoef.col = 'black',
           col=col2(100)[30:100])
}