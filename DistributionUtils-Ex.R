pkgname <- "DistributionUtils"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('DistributionUtils')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("besselRatio")
### * besselRatio

flush(stderr()); flush(stdout())

### Name: Bessel K Ratio
### Title: Ratio of Bessel K Functions
### Aliases: besselRatio
### Keywords: math

### ** Examples

nus <- c(0:5, 10, 20)
x <- seq(1, 4, length.out = 11)
k <- 3

raw <- matrix(nrow = length(nus), ncol = length(x))
scaled <- matrix(nrow = length(nus), ncol = length(x))
compare <- matrix(nrow = length(nus), ncol = length(x))

for (i in 1:length(nus)){
    for (j in 1:length(x)) {
        raw[i,j] <- besselRatio(x[j], nus[i],
                                orderDiff = k)
        scaled[i,j] <- besselRatio(x[j], nus[i],
                                orderDiff = k, useExpScaled = 1)
        compare[i,j] <- raw[i,j]/scaled[i,j]
    }
}
raw
scaled
compare




cleanEx()
nameEx("distCalcRange")
### * distCalcRange

flush(stderr()); flush(stdout())

### Name: distCalcRange
### Title: Range of a Unimodal Distribution
### Aliases: distCalcRange
### Keywords: distribution univar

### ** Examples

normRange <- distCalcRange("norm", tol = 10^(-7), mean = 4, sd = 1)
normRange
tRange <- distCalcRange("t", tol = 10^(-5), df = 4)
tRange



cleanEx()
nameEx("distIneqMassart")
### * distIneqMassart

flush(stderr()); flush(stdout())

### Name: distIneqMassart
### Title: Massart Inequality for Distributions
### Aliases: distIneqMassart
### Keywords: distribution univariate

### ** Examples

## Normal distribution is the default
distIneqMassart()
## Specify parameter values
distIneqMassart(mean = 1, sd = 2)
## Gamma distribution has no default value for shape
distIneqMassart("gamma", shape = 1)



cleanEx()
nameEx("distIneqMassartPlot")
### * distIneqMassartPlot

flush(stderr()); flush(stdout())

### Name: distIneqMassartPlot
### Title: Massart Inequality Plot Function
### Aliases: distIneqMassartPlot
### Keywords: distribution univar

### ** Examples

### The Massart Inequality plot for standard Normal Distribution
distIneqMassartPlot()

### The Massart Inequality plot for Gamma Distribution
distIneqMassartPlot("gamma", shape = 1)



cleanEx()
nameEx("distMode")
### * distMode

flush(stderr()); flush(stdout())

### Name: distMode
### Title: Mode of a Unimodal Distribution
### Aliases: distMode
### Keywords: distribution univar

### ** Examples

normRange <- distCalcRange("norm", tol = 10^(-7), mean = 4, sd = 1)
curve(dnorm(x, mean = 4, sd = 1), normRange[1], normRange[2])
abline(v = distMode("norm", mean = 4, sd = 1), col = "blue")



cleanEx()
nameEx("distStepSize")
### * distStepSize

flush(stderr()); flush(stdout())

### Name: distStepSize
### Title: Step Size for Calculating the Range of a Unimodal Distribution
### Aliases: distStepSize
### Keywords: distribution univar

### ** Examples

normRange <- distCalcRange("norm", tol = 10^(-7), mean = 4, sd = 1)
normRange
tRange <- distCalcRange("t", tol = 10^(-5), df = 4)
tRange



cleanEx()
nameEx("incompleteBesselK")
### * incompleteBesselK

flush(stderr()); flush(stdout())

### Name: incompleteBesselK
### Title: The Incomplete Bessel K Function
### Aliases: incompleteBesselK incompleteBesselKR SSFcoef combinatorial
###   GDENOM GNUM
### Keywords: math distribution

### ** Examples

### Harris (2008) gives accurate values (16 figures) for
### x = 0.01, y = 4, and nu = 0:9
### nu = 0, Harris value is 2.22531 07612 66469
options(digits = 16)
incompleteBesselK(0.01, 4, 0)
### nu = 9, Harris value is 0.00324 67980 03149
incompleteBesselK(0.01, 4, 9)

### Other values given in Harris (2008)
### x = 4.95, y = 5.00, nu = 2
incompleteBesselK(4.95, 5, 2) ## 0.00001 22499 87981
### x = 10, y = 2, nu = 6
### Slevinsky and Safouhi (2010) suggest Harris (2008) value
### is incorrect, give value 0.00000 04150 01064 21228
incompleteBesselK(10, 2, 6)
### x = 3.1, y = 2.6, nu = 5
incompleteBesselK(3.1, 2.6, 5) ## 0.00052 85043 25244

### Check values when x > y using numeric integration
(numIBF <- sapply(0:9, incompleteBesselK, x = 4, y = 0.01))

besselFn <- function(t, x, y, nu) {
  (t^(-nu - 1))*exp(-x*t - y/t)
}

(intIBF <- sapply(0:9, integrate, f = besselFn, lower = 1, upper = Inf,
                 x = 4, y = 0.01))
intIBF <- as.numeric(intIBF[1, ])
numIBF - intIBF
max(abs(numIBF - intIBF)) ## 1.256649992398273e-11

options(digits = 7)



cleanEx()
nameEx("integrateDens")
### * integrateDens

flush(stderr()); flush(stdout())

### Name: integrateDens
### Title: Integrates a Density Function
### Aliases: integrateDens
### Keywords: distribution univar

### ** Examples

integrateDens("norm", mean = 1, sd = 1)
integrateDens("t", df = 4)
integrateDens("exp", rate = 2)
integrateDens("weibull", shape = 1)



cleanEx()
nameEx("inversionTests")
### * inversionTests

flush(stderr()); flush(stdout())

### Name: inversionTests
### Title: Inversion Tests for Distributions
### Aliases: inversionTestpq inversionTestqp
### Keywords: distribution univariate

### ** Examples

## Default distribution is normal
inversionTestpq()
inversionTestqp()
## Supply parameters
inversionTestpq(mean = 1, sd = 2)
inversionTestqp(mean = 1, sd = 2)
## Gamma distribution, must specify shape
inversionTestpq("gamma", shape = 1)
inversionTestqp("gamma", shape = 1)



cleanEx()
nameEx("is.wholenumber")
### * is.wholenumber

flush(stderr()); flush(stdout())

### Name: is.wholenumber
### Title: Is Object Numeric and Whole Numbers
### Aliases: is.wholenumber
### Keywords: classes

### ** Examples

is.wholenumber(-3:5)                           # TRUE
is.wholenumber(c(0,0.1,1.3,5))                 # FALSE
is.wholenumber(-3:5 + .Machine$double.eps)     # TRUE
is.wholenumber(-3:5 + .Machine$double.eps^0.5) # FALSE
is.wholenumber(c(2L,3L))                       # TRUE
is.wholenumber(c("2L","3L"))                   # FALSE
is.wholenumber(0i ^ (-3:3))                    # FALSE
is.wholenumber(matrix(1:6, nrow = 3))          # TRUE
is.wholenumber(list(-1:3,2:6))                 # FALSE
is.numeric(list(-1:3,2:6))                     # FALSE
is.wholenumber(unlist(list(-1:3,2:6)))         # TRUE



cleanEx()
nameEx("logHist")
### * logHist

flush(stderr()); flush(stdout())

### Name: logHist
### Title: Plot Log-Histogram
### Aliases: logHist
### Keywords: hplot distribution

### ** Examples

x <- rnorm(200)
hist(x)
### default
logHist(x)
### log histogram only
logHist(x, htype = "h")
### points only, some options
logHist(x, htype = "p", pch = 20, cex = 2, col = "steelblue")



cleanEx()
nameEx("momChangeAbout")
### * momChangeAbout

flush(stderr()); flush(stdout())

### Name: momChangeAbout
### Title: Obtain Moments About a New Location
### Aliases: momChangeAbout
### Keywords: distribution univar

### ** Examples

### Gamma distribution
k <- 4
shape <- 2
old <- 0
new <- 1
sampSize <- 1000000

### Calculate 1st to 4th raw moments 
m <- numeric(k)
for (i in 1:k){
   m[i] <- gamma(shape + i)/gamma(shape)
}
m

### Calculate 4th moment about new 
momChangeAbout(k, m, old, new)
### Calculate 3rd about new
momChangeAbout(3, m, old, new)

### Calculate 1st to 4th moments about new
momChangeAbout(oldMom = m, oldAbout = old, newAbout = new)
momChangeAbout(order = "all", m, old, new)
  
### Approximate kth moment about new using sampling
x <- rgamma(sampSize, shape)
mean((x - new)^k) 



cleanEx()
nameEx("momIntegrated")
### * momIntegrated

flush(stderr()); flush(stdout())

### Name: momIntegrated
### Title: Moments Using Integration
### Aliases: momIntegrated
### Keywords: distribution univar

### ** Examples

require(GeneralizedHyperbolic)
### Calculate the mean of a generalized hyperbolic distribution
### Compare the use of integration and the formula for the mean
m1 <- momIntegrated("ghyp", param = c(0, 1, 3, 1, 1 / 2), order = 1, about = 0)
m1
ghypMean(param = c(0, 1, 3, 1, 1 / 2))
### The first moment about the mean should be zero
momIntegrated("ghyp", order = 1, param = c(0, 1, 3, 1, 1 / 2), about = m1)
### The variance can be calculated from the raw moments
m2 <- momIntegrated("ghyp", order = 2, param = c(0, 1, 3, 1, 1 / 2), about = 0)
m2
m2 - m1^2
### Compare with direct calculation using integration
momIntegrated("ghyp", order = 2, param = c(0, 1, 3, 1, 1 / 2), about = m1)
momIntegrated("ghyp", param = c(0, 1, 3, 1, 1 / 2), order = 2,
              about = m1)
### Compare with use of the formula for the variance
ghypVar(param = c(0, 1, 3, 1, 1 / 2))



cleanEx()
nameEx("momSE")
### * momSE

flush(stderr()); flush(stdout())

### Name: momSE
### Title: Standard Errors of Sample Moments
### Aliases: momSE
### Keywords: distribution univar

### ** Examples

### Moments of the normal distribution, mean 1, variance 4
mu <- 1
sigma <- 2
mom <- c(0,sigma^2,0,3*sigma^4,0,15*sigma^6,0,105*sigma^8)
### standard error of sample variance
momSE(2, 100, mom[1:4])
### should be
sqrt(2*sigma^4)/10
### standard error of sample central third moment
momSE(3, 100, mom[1:6])
### should be
sqrt(6*sigma^6)/10
### standard error of sample central fourth moment
momSE(4, 100, mom)
### should be
sqrt(96*sigma^8)/10



cleanEx()
nameEx("moranTest")
### * moranTest

flush(stderr()); flush(stdout())

### Name: moranTest
### Title: Moran's Log Spacings Test
### Aliases: moranTest
### Keywords: distribution univariate

### ** Examples


### Normal Distribution
x <- rnorm(100, mean = 0, sd = 1)
muhat <- mean(x)
sigmahat <- sqrt(var(x)*(100 - 1)/100)
result <- moranTest(x, "norm", mean = muhat, sd = sigmahat)
result

### Exponential Distribution
y <- rexp(200, rate = 3)
lambdahat <- 1/mean(y)
result <- moranTest(y, "exp", rate = lambdahat)
result



cleanEx()
nameEx("pDist")
### * pDist

flush(stderr()); flush(stdout())

### Name: pDist
### Title: Distribution and Quantile Functions for Unimodal Distributions
### Aliases: pDist qDist
### Keywords: distribution univar

### ** Examples

pDist("norm", q = 2, mean = 1, sd = 1)
pDist("t", q = 0.5, df = 4)
require(GeneralizedHyperbolic)
pDist("ghyp", q = 0.1)
require(SkewHyperbolic)
qDist("skewhyp", p = 0.4, param = c(0, 1, 0, 10))
qDist("t", p = 0.2, df = 4)



cleanEx()
nameEx("safeIntegrate")
### * safeIntegrate

flush(stderr()); flush(stdout())

### Name: safeIntegrate
### Title: Safe Integration of One-Dimensional Functions
### Aliases: safeIntegrate print.integrate
### Keywords: math utilities

### ** Examples

integrate(dnorm, -1.96, 1.96)
safeIntegrate(dnorm, -1.96, 1.96)  # Same as for integrate()
integrate(dnorm, -Inf, Inf)
safeIntegrate(dnorm, -Inf, Inf)    # Same as for integrate()
integrate(dnorm, 1.96, 1.96)       # OK here but can give an error
safeIntegrate(dnorm, 1.96, 1.96)
integrate(dnorm, -Inf, -Inf)
safeIntegrate(dnorm, -Inf, -Inf)   # Avoids nonsense answer
integrate(dnorm, Inf, Inf)
safeIntegrate(dnorm, Inf, Inf)     # Avoids nonsense answer



cleanEx()
nameEx("sampleMoments")
### * sampleMoments

flush(stderr()); flush(stdout())

### Name: Sample Moments
### Title: Sample Skewness and Kurtosis
### Aliases: skewness kurtosis
### Keywords: univar

### ** Examples

x <- rnorm(100)
skewness(x)
kurtosis(x)



cleanEx()
nameEx("tailPlot")
### * tailPlot

flush(stderr()); flush(stdout())

### Name: tailPlot
### Title: Tail Plot Functions
### Aliases: tailPlot tailPlotLine normTailPlotLine tTailPlotLine
###   gammaTailPlotLine
### Keywords: distribution univar

### ** Examples

### Draw tail plot of some data
x <- rnorm(100, 1, 2)
tailPlot(x)
### Add normal distribution line
normTailPlotLine(x, mean = 1, sd = 2)
### Add t distribution line
tTailPlotLine(x, df = 5, lty = 2)
### Use fitted values
normTailPlotLine(x, mean = mean(x), sd = sd(x), lty = 3)

### Gamma distribution
x <- rgamma(100, shape = 1, scale = 1)
tailPlot(x)
### Add gamma distribution line
gammaTailPlotLine(x, shape = 1, scale = 1)
### Left tail example
tailPlot(x, side = "l")
### Add gamma distribution line
gammaTailPlotLine(x, shape = 1, scale = 1, side = "l")
### Log scale on both axes
tailPlot(x, side = "l", log = "xy")
### Add gamma distribution line
gammaTailPlotLine(x, shape = 1, scale = 1, side = "l")

### Add line from a standard distribution with default parameters
x <- rlnorm(100)
tailPlot(x)
tailPlotLine(x, distrFn = "lnorm")

### Add line from a distribution with 'param' argument
require(VarianceGamma)
param <- c(0,0.5,0,0.5)
x <- rvg(100, param = param)
tailPlot(x)
tailPlotLine(x, distrFn = "vg", param = param)




cleanEx()
nameEx("tsHessian")
### * tsHessian

flush(stderr()); flush(stdout())

### Name: tsHessian
### Title: Calculate Two-Sided Hessian Approximation
### Aliases: tsHessian
### Keywords: math

### ** Examples

### Consider Hessian of log(1 + x + 2y)
### Example from Lang: A Second Course in Calculus, p.74
fun <- function(param){
  x <- param[1]
  y <- param[2]
  return(log(1 + x + 2*y))
}

### True value of Hessian at (0,0)
trueHessian <- matrix( c(-1,-2,
                         -2,-4), byrow = 2, nrow = 2)
trueHessian

### Value from tsHessian
approxHessian <- tsHessian(c(0,0), fun = fun)
approxHessian
maxDiff <- max(abs(trueHessian - approxHessian))
### Should be approximately 0.045
maxDiff




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
