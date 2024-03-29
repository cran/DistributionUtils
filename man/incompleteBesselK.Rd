\name{incompleteBesselK}
\title{The Incomplete Bessel K Function}
\alias{incompleteBesselK}
\alias{incompleteBesselKR}
\alias{SSFcoef}
\alias{combinatorial}
\alias{GDENOM}
\alias{GNUM}
\description{
  Calculates the incomplete Bessel K function using the algorithm and code
  provided by Slavinsky and Safouhi (2009).
}
\usage{
incompleteBesselK (x, y, nu, tol = .Machine$double.eps^0.85, nmax = 120)
incompleteBesselKR(x, y, nu, tol = .Machine$double.eps^0.85, nmax = 120)
SSFcoef(nmax, nu)
combinatorial(nu)
GDENOM(n, x, y, nu, An, nmax, Cnp)
GNUM(n, x, y, nu, Am, An, nmax, Cnp, GM, GN)
}
\arguments{
  \item{x}{Numeric. Value of the first argument of the incomplete Bessel K
    function.}
  \item{y}{Numeric. Value of the  second argument of the incomplete Bessel K
    function.}
  \item{nu}{Numeric. The order of the incomplete Bessel K function.}
  \item{tol}{Numeric. The tolerance for the difference between
    successive approximations of the incomplete Bessel K function.}
  \item{nmax}{Integer. The maximum order allowed for the approximation
    of the incomplete Bessel K function.}
  \item{n}{Integer. Current order of the approximation. Not required to be
    specified by users.}
  \item{An}{Matrix of coefficients. Not required to be specified by
    users.}
  \item{Am}{Matrix of coefficients. Not required to be specified by
    users.}
  \item{Cnp}{Vector of elements of Pascal's triangle.  Not required to
    be specified by users.}
  \item{GN}{Vector of denominators used for approximation.  Not required
    to be specified by users.}
  \item{GM}{Vector of numerators used for approximation.  Not required
    to be specified by users.}
}
\details{
  The function \code{incompleteBesselK} implements the algorithm
  proposed by Slavinsky and Safouhi (2010) and uses code provided by
  them.

  The incomplete Bessel K function is defined by

  \deqn{K_\nu(x,y)=\int_1^\infty t^{-nu-1}\exp(-xt-y/t)\,dt}{%
    K_nu(x,y)=integral_1^infinity t^{-nu-1} exp(-xt-y/t) dt}

  see Slavinsky and Safouhi (2010), or Harris (2008).

  \code{incompleteBesselK} calls a Fortran routine to carry out the
  calculations. \code{incompleteBesselKR()} is a pure R version of the
  routine for computing the incomplete Bessel K function.

  The functions \code{SSFcoef}, \code{combinatorial}, \code{GDENOM}, and
  \code{GNUM} are \dQuote{subroutines}, i.e., auxiliary functions used in
  \code{incompleteBesselKR()}.  They are not expected to be called by the
  user.

  The approximation to the incomplete Bessel K function returned by
  \code{incompleteBesselK} is highly accurate.  The default value of
  \code{tol} is about 10^(-14) on a 32-bit computer.  It appears that
  even higher accuracy is possible when \code{x > y}.  Then the tolerance
  can be taken as \code{.Machine$double.eps} and the number of correct
  figures essentially coincides with the number of figures representable
  in the machine being used.

  \code{incompleteBesselKR} is very slow compared to the Fortran version
  and is only included for those who wish to see the algorithm in \R
  rather than Fortran.
}

\value{
  \code{incompleteBesselK} and \code{incompleteBesselKR} both return an
  approximation to the incomplete Bessel K function as defined above.
}

\references{
  Harris, Frank E. (2008)
  Incomplete Bessel, generalized incomplete gamma, or leaky aquifer functions.
  \emph{J. Comp. Appl. Math.} \bold{215}, 260--269.

  Slevinsky, Richard M., and Safouhi, Hassan (2009)
  New formulae for higher order derivatives and applications.
  \emph{J. Comp. Appl. Math.} \bold{233}, 405--419.

  Slevinsky, Richard M., and Safouhi, Hassan (2010)
  A recursive algorithm for the G transformation and accurate
  computation of incomplete Bessel functions.
  \emph{Applied Numerical Mathematics} \bold{60}(12), 1411--1417.
}
\author{
  David Scott \email{d.scott@auckland.ac.nz}, Thomas Tran, Richard
  Slevinsky, Hassan Safouhi.
}
\note{
  The problem of calculation of the incomplete Bessel K function is
  equivalent to the problem of calculation of the cumulative
  distribution function of the generalized inverse Gaussian
  distribution. See \code{\link[GeneralizedHyperbolic:dgig]{Generalized
  Inverse Gaussian}.}
}
\seealso{
\code{\link{besselK}}
}
\examples{
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
}
\keyword{math}
\keyword{distribution}

