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
# devtools::install_github("jorgetendeiro/GGUM")
library(GGUM)
# END SECTION

# 1. Set up parameters ----
# 
set.seed(739)
N <- 98
I <- 4
C <- 3 # sample(2:3, I, replace = TRUE)
# END SECTION

# 2. Generate data ----
data.gen <- GenData.GGUM(N, I, C, model = "GUM", seed = 4367)
data     <- data.gen$data
# Data with missing values, to test the code:
pos.NA   <- matrix(rbinom(N * I, 1, .80), nrow = N)
pos.NA[pos.NA == 0] <- NA
data.NA  <- data * pos.NA
data <- data.NA
# END SECTION

# 3. Run GUM ----
if (max(C) - min(C) == 0)
{
  Model3.res <- GUM       (data, C)
  Th.GUM     <- Theta.EAP (Model3.res)
  # Plots: 
  # Category response curves per item:
  plotCRC(Model3.res, items = c(1, 3))
  # Test characteristic curve:
  plotTCC(Model3.res, Th.GUM)
  # Item characteristic curves:
  plotICC(Model3.res, Th.GUM, items = 3)
  # Test information:
  plotTIF(Model3.res, Th.GUM)
  # Item information:
  plotIIF(Model3.res, Th.GUM, items = 3)
}
# END SECTION

# 4. Run GGUM ----
Model8.res <- GGUM      (data, C)
Th.GGUM    <- Theta.EAP (Model8.res)
# Plots: 
# Category response curves per item:
plotCRC(Model8.res, items = c(1, 3))
# Test characteristic curve:
plotTCC(Model8.res, Th.GGUM)
# Item characteristic curves:
plotICC(Model8.res, Th.GGUM, items = c(1, 3))
# Test information:
plotTIF(Model8.res, Th.GGUM)
# Item information:
plotIIF(Model8.res, Th.GGUM, items = c(1, 3))
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
export.GGUM2004(data)
write.GGUM2004(I, C, model = "GUM")
run.GGUM2004("cmd.txt", I, C, N, model = "GGUM" )


# END SECTION



