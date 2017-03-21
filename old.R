source('import_earthquakes.R')
source('pois_hmm.R')
source('norm_hmm.R')

library('plotrix')

cumulativeForecast <- function(xx, h.max) {
  T <- length(xx)
  xx.range <- range(xx)
  plot(xx, type='o', pch=20,
       xlim=c(0, T + h.max), ylim=xx.range)
  hmm <- PoisHmm(3, xx)
  for (h in 1:h.max) {
    hmm$fit(xx)
    F <- hmm$forecastDist(xx, 0:50, 1)
    print(which.max(F))
    xx <- c(xx, which.max(F[1,]))
  }
  par(new=TRUE)
  plot((T+1):(T+h.max), xx[(T+1):length(xx)],
       type='o', pch=20, col='red',
       xlim=c(0, T + h.max), ylim=xx.range)
  abline(v=T, lty=2)
}

cumulativeForecast(earthquakes.zuc.annual$count, 16)

trainForecast <- function(HmmClass, m, xx, x.probes, h.max, conf.level = 0.90) {
  T <- length(xx)
  T.train <- T-h.max
  xx.train <- xx[1:(T.train)]
  hmm <- HmmClass(m, xx.train)
  hmm$fit(xx.train)
  F <- hmm$forecastDist(xx.train, x.probes, h.max)
  hmm$plotForecast(F, xx, x.probes)
  
  plot(xx.train, type='o', pch=20,
       xlim=c(1, T), ylim=range(xx))
  par(new=TRUE)
  
  plot((T.train+1):T, x.probes[apply(F, 1, which.max)], col='red',
       type='o', pch=20, axes=FALSE,
       xlim=c(1, T), ylim=range(xx))
  abline(v=T.train, lty=2)
  par(new=TRUE)
  plot((T.train+1):T, xx[(T.train+1):T],
       type='o', pch=1,
       xlim=c(1, T), ylim=range(xx))
  
  # close-up
  plot((T.train+1):T, xx[(T.train+1):T], col='green', type='o', pch=20, ylim=range(xx))
  par(new=TRUE)
  plot((T.train+1):T, x.probes[apply(F, 1, which.max)],
       type='o', pch=1, axes=FALSE, ylim=range(xx))
  # , xaxt = 'n')
  # axis(1, at=1:h.max, labels=(length(xx)-h.max+1):length(xx))
  
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
         add=TRUE, pch=NA)
}

xx <- earthquakes.zuc.annual$count
x.probes <- 0:50
trainForecast(PoisHmm, 4, xx, x.probes, 16, 0.9)

xx <- earthquakes.usgs.annual$mag.mean
x.probes <- seq(min(xx), max(xx), length.out = 50)
trainForecast(NormHmm, 3, xx, x.probes, 16, 0.9)

xx <- earthquakes.usgs.annual$mag.mean
h.max <- 32
T <- length(xx)
T.train <- T-h.max
xx.train <- xx[1:(T.train)]
m.opt <- min(selectModelIc(NormHmm, xx, 6, TRUE))

hmm <- NormHmm(m.opt, xx.train)
hmm$fit(xx.train)
plot(xx.train, type='o', pch=20,
     xlim=c(1, T), ylim=range(xx))
par(new=TRUE)
plot((T.train+1):T, hmm$genSample(h.max),
     type='o', pch=1, axes=FALSE,
     xlim=c(1, T), ylim=range(xx))
abline(v=T.train, lty=2)
par(new=TRUE)
plot((T.train+1):T, xx[(T.train+1):T], col='green',
     type='o', pch=20,
     xlim=c(1, T), ylim=range(xx))

x.probes <- seq(min(xx), max(xx), length.out = 50)
trainForecast(NormHmm, m.opt, xx, x.probes, h.max, 0.9)