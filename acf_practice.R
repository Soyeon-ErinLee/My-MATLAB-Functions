library(astsa)
acf = ARMAacf(ar=c(1.5,-0.75),ma=0,50)
plot(acf,type="h",xlab="lag")
abline(h=0)


set.seed(5)
ar2 <- arima.sim(list(order=c(2,0,0), ar = c(1.5, -0.75)), n=144)
plot.ts(ar2, axes=F); box(); axis(2)
axis(1,seq(0,144,24))
abline(v=seq(0,144,12),lty="dotted")
Periodogram <- abs(fft(ar2)/144)**2
frequency = 2*pi*c(0:143)/144
plot(frequency, Periodogram,type="o")
arma.spec( ar = c(1.5, -0.75), log = "no", main = "Autoregressive")
