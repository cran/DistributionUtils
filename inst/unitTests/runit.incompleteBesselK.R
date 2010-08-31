### Unit tests of function incompleteBesselK

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make

test.incompleteBesselK <- function()
{
  ## Purpose: Level 1 test of incompleteBesselK
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Local, Date: 31 August, 2010

  ## Data from Harris (2008)
  HarrisCase1 <- c(2.225310761266469,
                   0.213894166822940,
                   0.054503469799701,
                   0.023253121507708,
                   0.013042750996080,
                   0.008567534990649,
                   0.006208676806601,
                   0.004801085238177,
                   0.003884072049627,
                   0.003246798003149)

  HarrisCase2 <- 0.000012249987981

  HarrisCase3 <- 0.000000415004594

  HarrisCase4 <- 0.000528504325244
                   
  ## Calculations for Harris Case 1
  numIBF <- numeric(10)
  for (i in 0:9){
    ibf <- incompleteBesselK(0.01, 4, i, nmax = 100)
    numIBF[i + 1] <- ibf
  }
  maxDiff <- max(abs(numIBF - HarrisCase1))
  checkTrue(maxDiff <= 10^(-13), msg = "Harris Case 1")

  ## Calculations for Harris Case 2
  ibf <- incompleteBesselK(4.95, 5, 2)
  maxDiff <- max(abs(ibf - HarrisCase2))
  checkTrue(maxDiff <= 10^(-14), msg = "Harris Case 2")
  
  ## Calculations for Harris Case 3
  ibf <- incompleteBesselK(10, 2, 6)
  maxDiff <- max(abs(ibf - HarrisCase3))
  checkTrue(maxDiff <= 10^(-14), msg = "Harris Case 3")

  ## Calculations for Harris Case 4
  ibf <- incompleteBesselK(3.1, 2.6, 5)
  maxDiff <- max(abs(ibf - HarrisCase4))
  checkTrue(maxDiff <= 10^(-13), msg = "Harris Case 4")
  

  return()
}

