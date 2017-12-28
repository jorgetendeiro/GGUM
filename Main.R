# TABLE OF CONTENTS
# 0. Prepare environment
# 1. Set up parameters
# 2. Generate data
# 3. Run GUM
# 4. Run GGUM
# 5. Compare generated and estimated parameters
# 6. Test GGUM2004 related functions
#

# 0. Prepare environment ----
rm(list=ls())
if (!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])
library(psych) 
library(abind)
library(viridis) # new color palettes

library(devtools)
library(roxygen2)
library(knitr)
load_all()
# END SECTION

# 1. Set up parameters ----
# 
set.seed(738)
N <- 500
I <- 5
C <- sample(2:6, I, replace = TRUE)
# END SECTION

# 2. Generate data ----
data.gen <- GenData.GGUM(N, I, C, model = "GGUM", seed = 4367)
data     <- data.gen$data
# Data with missing values, to test the code:
pos.NA   <- matrix(rbinom(N * I, 1, .80), nrow = N)
pos.NA[pos.NA == 0] <- NA
data.NA  <- data * pos.NA
# END SECTION

# 3. Run GUM ----
Model3.IP    <- GUM            (data, C)
Model3.Th    <- Theta.EAP         (Model3.IP, TRUE)
Model3.Th.SE <- Theta.SE          (data, C, Model3.IP, Model3.Th)
Model3.IP.SE <- ItemParamModel3.SE(data, C, Model3.IP)
# Plots: 
# Response category curves per item:
plot.GGUM(C, Model3.IP, items = NULL, x.lim = 4, ThetaminDelta = FALSE)
# Test characteristic curve:
plot.TestCharacteristicCurve.GGUM(data, C, Model3.IP, Model3.Th)
# Item characteristic curves:
plot.ItemCharacteristicCurve.GGUM(data, C, Model3.IP, Model3.Th, items = NULL)
# Test information:
plot.TestInf(data, C, Model3.IP, Model3.Th)
# Item information:
plot.ItemInf(data, C, Model3.IP, Model3.Th, items = NULL)
# END SECTION

# 4. Run GGUM ----
Model8.res <- GGUM      (data, C)
Th.GGUM    <- Theta.EAP   (data, C, Model8.res)
Th.SE.GGUM <- Theta.SE    (data, C, Model8.res, Th.GGUM)
It.SE.GGUM <- ItemParamModel8.SE(data, C, Model8.res)
# Plots: 
# Response category curves per item:
plot.GGUM(C, Model8.res, items = NULL, x.lim = 4, ThetaminDelta = FALSE)
# Test characteristic curve:
plot.TestCharacteristicCurve.GGUM(data, C, Model8.res, Th.GGUM)
# Item characteristic curves:
plot.ItemCharacteristicCurve.GGUM(data, C, Model8.res, Th.GGUM, items = NULL)
# Test information:
plot.TestInf(data, C, Model8.res, Th.GGUM)
# Item information:
plot.ItemInf(data, C, Model8.res, Th.GGUM, items = NULL)
# END SECTION

# 5. Compare generated and estimated parameters ----
IP.est <- Model8.res
Th.est <- Th.GGUM$Th.full
# 
MAD.alpha  <- round(sum(abs(IP.est$alpha - data.gen$alpha.gen)) / I, 4)
BIAS.alpha <- round(sum(    IP.est$alpha - data.gen$alpha.gen)  / I, 4)
cor.alpha  <- round(cor(    IP.est$alpha,  data.gen$alpha.gen)     , 4)
# 
cor.delta  <- round(cor(    IP.est$delta,  data.gen$delta.gen)     , 4)
MAD.delta  <- round(sum(abs(IP.est$delta - sign(cor.delta)*data.gen$delta.gen)) / I, 4)
BIAS.delta <- round(sum(    IP.est$delta - sign(cor.delta)*data.gen$delta.gen)  / I, 4)
# 
MAD.taus   <- round(sum(abs(IP.est$taus[, 1:max(C)] -   data.gen$taus.gen[, 1:max(C)])) / sum(C), 4)
BIAS.taus  <- round(sum(    IP.est$taus[, 1:max(C)] -   data.gen$taus.gen[, 1:max(C)])  / sum(C), 4)
cor.taus   <- round(cor(  c(IP.est  $taus    [, 1:max(C)][IP.est  $taus    [, 1:max(C)] != 0]), 
                          c(data.gen$taus.gen[, 1:max(C)][data.gen$taus.gen[, 1:max(C)] != 0])) , 4)
# 
cor.th     <- round(cor(    Th.est,  data.gen$theta.gen)     , 4)
MAD.th     <- round(sum(abs(Th.est - sign(cor.th)*data.gen$theta.gen)) / N, 4)
BIAS.th    <- round(sum(    Th.est - sign(cor.th)*data.gen$theta.gen)  / N, 4)

# 6. Test GGUM2004 related functions ----
# Export data to GGUM2004:
export.GGUM2004(data.NA, file.name = "C:/GGUM2004/DataTest")
write.GGUM2004("inputTest", "C:/GGUM2004/DataTest.txt", I, C,
               cutoff = 2, model = "GUM"  )

run.GGUM2004("inputTest", I, C, N, model = "GUM" )


# END SECTION



