source('../util/import_earthquakes.R')
source('../models/pois_hmm.R')
source('../models/norm_hmm.R')

xx <- earthquakes.usgs.annual$count

# One plot
# (5.05, 3.38) * 1.2
# pdf(file="~/Desktop/Rplot01.pdf", width=6.06, height=4.06, pointsize=12)
# par(mar=c(4,2,2,1)+0.1)  # bottom, left, top, right
# plot(xx, type='o', pch=20, main='Title')
# dev.off()

# Two plots side by side
# pdf(file="~/Desktop/Rplot02.pdf", width=3.03, height=4.06, pointsize=12)
# par(mar=c(4,4,0,1)+0.1)  # bottom, left, top, right
# plot(xx, type='p', pch=20, ylab='ylab', xlab='xlab', lwd=0.5)
# dev.off()

# ---

pdf(file="~/Documents/СПбГУ/Диплом/pic/eq-count-ts.pdf",
    width=6.06, height=4.06, pointsize=12)
par(mar=c(4,4,0,1)+0.1)  # bottom, left, top, right
plot(earthquakes.usgs.annual$year, xx, type='o', pch=20,
     ylab='Earthquake count', xlab='Year')
dev.off()

# ---

pdf(file="~/Documents/СПбГУ/Диплом/pic/eq-count-ts-half.pdf",
    width=3.03, height=3.38, pointsize=12)
par(mar=c(4,4,0,1)+0.1)  # bottom, left, top, right
plot(earthquakes.usgs.annual$year, xx, type='p', pch=20,
     ylab='Earthquake count', xlab='Year', lwd=0.5)
dev.off()

mllk.pois.rv <- function(lambda.wrk, xx)
  -sum(dpois(xx, exp(lambda.wrk), log=TRUE))
par.opt <- nlm(mllk.pois.rv, 1, xx)
exp(par.opt$estimate)
par.opt$minimum

lambda=mean(xx)

pdf(file="~/Documents/СПбГУ/Диплом/pic/eq-count-hist-half.pdf",
    width=3.03, height=3.38, pointsize=12)
par(mar=c(4,4,0,1)+0.1)  # bottom, left, top, right
ymax <- dpois(floor(lambda), lambda) * 1.2
hist(xx, breaks=20, freq=FALSE,
     ylab='Density', xlab='Counts', col='grey', main='',
     xlim=c(0, max(xx)), ylim=c(0, ymax))
par(new=TRUE)
plot(dpois(0:max(xx), lambda), pch=16, axes=FALSE, type='l', col='red',
     ylab='', xlab='',
     xlim=c(0, max(xx)), ylim=c(0, ymax))
legend("topright", legend = c("Pois PDF"), bty = "n",
       col = c("red"), lty = c(1))
dev.off()

# ---

source('../models/pois-mixture.R')

# pdf(file="~/Documents/СПбГУ/Диплом/pic/pois-mixture-m2.pdf",
#     width=6.06, height=4.06, pointsize=12)
pdf(file="~/Documents/СПбГУ/Диплом/pic/pois-mixture-m2-half.pdf",
    width=3.03, height=3.38, pointsize=12)
par(mar=c(4,4,0,1)+0.1)  # bottom, left, top, right
par.init <- list(lambda=c(6, 18), delta=c(1, 1)/3)
m <- length(par.init$lambda)
par.opt.wrk <- nlm(mllk, nat2wrk(par.init), xx)
par.opt.nat <- wrk2nat(par.opt.wrk$estimate)
cat("mllk =", par.opt.wrk$minimum, "\n")
par.opt.nat
eta <- mixture(dpois, par.opt.nat$lambda, par.opt.nat$delta)
h <- hist(xx, breaks=20, plot=FALSE)
ymax <- max(dpois(floor(mean(xx)), mean(xx)), max(h$density))
plot(h, freq=FALSE, col="grey", ylab='', xlab='', axes=FALSE,
     xlim=c(0, max(xx)), ylim=c(0, ymax), main='')
par(new=TRUE)
plot(eta(0:max(xx)), type="l", col='red',
     xlab='Counts', ylab='Density',
     xlim=c(0, max(xx)), ylim=c(0, ymax))
legend("topright", legend = c("Mix. PDF, m = 2"), bty = "n",
       col = c("red"), lty = c(1))
cat("Model mean =", mixture.mean(par.opt.nat$lambda, par.opt.nat$delta), "\n")
cat("Model variance =", mixture.var(par.opt.nat$lambda, par.opt.nat$delta), "\n")
cat("Data mean: ", mean(xx), "\n")
cat("Data variance: ", var(xx) * (length(xx)-1) / length(xx), "\n")
dev.off()

# ---

pdf(file="~/Documents/СПбГУ/Диплом/pic/eq-acf.pdf",
    width=6.06, height=4.06, pointsize=12)
par(mar=c(4,4,0,0)+0.1)  # bottom, left, top, right
acf(xx)
dev.off()

# ---

source('../models/pois_hmm.R')

pdf(file="~/Documents/СПбГУ/Диплом/pic/pois-hmm-m2-half.pdf",
    width=3.03, height=3.38, pointsize=12)
par(mar=c(4,4,0,0)+0.1)  # bottom, left, top, right
m <- 2
hmm <- PoisHmm(m, xx)
print(hmm)
cat('---------- Parameter conversions ----------\n')
par.wrk <- hmm$getWrkParams()
hmm$setWrkParams(par.wrk)
print(hmm)
cat('---------- alphaNorm ----------\n')
print(hmm$alphaNorm(xx))
print(hmm$alphaNorm(xx, compute.L = TRUE))
cat('---------- fit ----------\n')
nlm.out <- hmm$fit(xx)
print(nlm.out)
print(hmm)
h <- hist(xx, breaks=20, plot=FALSE)
# h$density <- h$density / length(h$breaks)
ymax <- -Inf
for (j in 1:m) {
  ymax <- max(ymax, max(hmm$pdf(xmin:xmax, hmm$pdf.params$lambda[j])))
}
ymax <- ymax * 1.2
# ymax <- max(h$density) * 2
xmax <- max(xx)
xmin <- min(xx)
# hist(xx, col='grey', freq=FALSE,
#      main=paste("m = ", m, "; -log(L) = ", round(nlm.out$minimum, 4)),
#      ylim=c(0, ymax), xlim=c(xmin, xmax))
plot(h, col='grey', freq=FALSE, ylim=c(0, ymax), xlim=c(0, xmax),
     xlab='Counts', ylab='Density', main='')
colors <- c('red', 'blue')
for (j in 1:m) {
  par(new=TRUE)
  plot(hmm$pdf(xmin:xmax, hmm$pdf.params$lambda[j]),
       type='l', col=colors[j], axes=FALSE,
       ylab='', xlab='',
       ylim=c(0, ymax), xlim=c(0, xmax))
}
legend("topright", legend = c(paste("Pois(", round(hmm$pdf.params$lambda[1], 2), ") PDF", sep=''),
                              paste("Pois(", round(hmm$pdf.params$lambda[2], 2), ") PDF", sep='')),
       bty = "n",
       col = colors, lty = c(1))
dev.off()

#---

source('../test/test_hmm.R')

pdf(file="~/Documents/СПбГУ/Диплом/pic/pois-hmm-ms.pdf",
    width=6.06, height=4.06, pointsize=12)
par(mar=c(4,4,3,0)+0.1)  # bottom, left, top, right
testHmmFit(PoisHmm, xx, 4, TRUE)
dev.off()

source('../models/hmm.R')

pdf(file="~/Documents/СПбГУ/Диплом/pic/pois-eq-aic-bic-half.pdf",
    width=3.03, height=3.38, pointsize=12)
par(mar=c(4,2,0,0)+0.1)  # bottom, left, top, right
selectModelIc(PoisHmm, xx, 6, do.plot = TRUE)
dev.off()

# ---

source('../models/pois_hmm.R')
source('../models/norm_hmm.R')

pdf(file="~/Documents/СПбГУ/Диплом/pic/norm-hmm-eq-forecast-16.pdf",
    width=6.06, height=4.06, pointsize=12)
par(mar=c(4,4,3,1)+0.1)  # bottom, left, top, right

xx <- earthquakes.usgs.annual$mag.mean
T <- length(xx)
x.probes <- seq(min(xx), max(xx), length.out = 50)
h.max <- 16
T.train <- length(xx) - h.max
xx.train <- xx[1:(T.train)]
m <- 3
hmm <- NormHmm(m, xx.train)
hmm$fit(xx.train)
F <- hmm$forecast(xx.train, x.probes, h.max)
hmm$plotForecast(F, xx, x.probes, c(1,2,3,16))
dev.off()

pdf(file="~/Documents/СПбГУ/Диплом/pic/eq-mag-split-16.pdf",
    width=6.06, height=4.06, pointsize=12)
par(mar=c(4,4,0,0)+0.1)  # bottom, left, top, right
par(xaxt='n')
plot(xx.train,
     type='o', pch=20,
     xlab='Year', ylab='Mean magnitude',
     xlim=c(1, T), ylim=range(xx))
par(xaxt='s')
axis(1, at=seq(1, length(xx), by=10),
     labels=earthquakes.usgs.annual$year[seq(1, length(xx), by=10)])
abline(v=T.train, lty=2)
par(new=TRUE)
plot((T.train+1):T, xx[(T.train+1):T],
     type='o', pch=1, xlab='', ylab='', axes=FALSE,
     xlim=c(1, T), ylim=range(xx))
dev.off()

pdf(file="~/Documents/СПбГУ/Диплом/pic/eq-mag-split-16-forecast.pdf",
    width=6.06, height=4.06, pointsize=12)
par(mar=c(4,4,0,0)+0.1)  # bottom, left, top, right
plot(xx.train,
     type='o', pch=20,
     xlab='Year', ylab='Mean magnitude',
     xaxt='n',
     xlim=c(1, T), ylim=range(xx))
axis(1, at=seq(1, length(xx), by=10),
     labels=earthquakes.usgs.annual$year[seq(1, length(xx), by=10)])
abline(v=T.train, lty=2)
par(new=TRUE)
plot((T.train+1):T, xx[(T.train+1):T],
     type='o', pch=1, xlab='', ylab='', axes=FALSE,
     xlim=c(1, T), ylim=range(xx))
par(new=TRUE)
plot((T.train+1):T, x.probes[apply(F, 1, which.max)], col='red',
     type='o', pch=20, axes=FALSE, xlab='', ylab='',
     xlim=c(1, T), ylim=range(xx))
dev.off()


pdf(file="~/Documents/СПбГУ/Диплом/pic/eq-forecast-16-closeup.pdf",
    width=6.06, height=4.06, pointsize=12)
par(mar=c(4,4,0,0)+0.1)  # bottom, left, top, right
# close-up
plot(earthquakes.usgs.annual$year[(T.train+1):T], xx[(T.train+1):T],
     type='o', pch=21, ylim=range(xx),
     ylab='Mean magnitude', xlab='Year')
par(new=TRUE)
plot((T.train+1):T, x.probes[apply(F, 1, which.max)], col='red',
     type='o', pch=20, axes=FALSE, ylim=range(xx), xlab='', ylab='')

conf.level <- 0.9
ui <- numeric(h.max)
li <- numeric(h.max)
for (h in 1:h.max) {
  print(h)
  am <- which.max(F[h,])
  acc <- F[h, am]
  ui[h] <- am
  li[h] <- am
  while ((li[h] > 1) || (ui[h] < length(F[h,]))) {
    if (li[h] > 1)  {
      if ((F[h, li[h] - 1] + acc) > conf.level) {
        break
      }
      li[h] <- li[h] - 1
      acc <- acc + F[h, li[h]]
    }
    if (ui[h] < length(F[h,])) {
      if ((F[1, ui[h] + 1] + acc) > conf.level) {
        break
      }
      ui[h] <- ui[h] + 1
      acc <- acc + F[h, ui[h]]
    }
    print(c(li[h], ui[h], acc))
  }
}
plotCI((T.train+1):T, x.probes[apply(F, 1, which.max)],
       ui=x.probes[ui], li=x.probes[li], ylim=range(xx),
       slty='dashed',
       add=TRUE, pch=NA)
legend("topright", legend = c("Observed values", "HMM forecast", "90% conf. int."),
       bty = "n",
       col = c('black', "red", 'black'),
       lty = c(1, 1, 2))
dev.off()