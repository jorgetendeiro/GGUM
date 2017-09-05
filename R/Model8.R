Model8 <- function(data, C, method = "MML", 
                   N.nodes = 30, max.outer = 60, max.inner = 60, tol = .001, precision = 4, verbose = "intermediate")
{
  data.condensed  <- GGUM.data.condense(data)$data.condensed
  I               <- ncol(data.condensed) - 1
  S               <- nrow(data.condensed)
  if (length(C) == 1) {C <- rep(C, I)}
  C.max           <- max(C)
  M               <- 2 * C + 1
  
  # Initial values, Step 1 of 2 - derived from Model 3:
  Model3.res     <- Model3(data, C, verbose = "intermediate")
  delta.old      <- Model3.res$delta
  taus.old       <- Model3.res$taus
  rm(Model3.res)
  # Initial values, Step 2 of 2 - derived from Model 4, based on the estimates from Model 3 as initial values:
  Model4.res     <- Model4(data, C, verbose = "intermediate", Init.vals = list(delta = delta.old, taus = taus.old))
  delta.old      <- Model4.res$delta
  taus.old       <- Model4.res$taus
  rm(Model4.res)
  # 
  alpha.old      <- rep(1, I)
  alphadelta.old <- rbind(alpha.old, delta.old)
  
  # Nodes and weights:
  nodes   <- seq(-4, 4, length.out = N.nodes)
  weights <- dnorm(nodes) / sum(dnorm(nodes))
  # 
  rs.arr      <- array(rep(data.condensed[, I + 1], N.nodes * I * (C.max + 1)), dim = c(S, N.nodes, I, C.max + 1)) # S x N.nodes x I x (C+1)
  weights.arr <- array(rep(weights, S * I * (C.max + 1)), dim = c(N.nodes, S, I, C.max + 1))
  weights.arr <- aperm(weights.arr, c(2, 1, 3, 4))                                                                 # S x N.nodes x I x (C+1)
  H.siz       <- array(0, dim = c(S, I, C.max + 1))
  for (s in 1:S) {
    for (i in 1:I) {
      H.siz[s, i, data.condensed[s, i] + 1] <- 1
    }
  }
  rm(s, i)
  H.siz  <- array(rep(H.siz, N.nodes), dim = c(S, I, C.max + 1, N.nodes))
  H.siz  <- aperm(H.siz, c(1, 4, 2, 3))                                                                            # S x N.nodes x I x (C+1)
  # 
  
  iter.outer <- 0
  iter.inner <- 10 # just to pass the first 'while' test below.
  
  if (verbose == "intermediate") {cat("\nStep 3 of 3: Estimating item parameters for GGUM... \n")}
  
  while ((iter.outer <= (max.outer - 1)) && (iter.inner > 1))
  {
    if (verbose == "yes") {cat("\niter.outer = ", iter.outer + 1, "\n")}
    if (verbose == "intermediate")
    {
      cat("\r", "|", rep("-", iter.outer+1), rep(" ", max.outer - iter.outer-1), "|", sep = "")
      flush.console()
    }
    
    iter.inner <- 0
    curr.tol   <- 1
    
    # E stage: Compute r.bar.izf and N.bar.if, 
    #             given the current alpha, delta, and taus.
    # Ls:
    Ls.mat <- Ls(data.condensed, alpha.old, delta.old, taus.old, nodes, C)
    
#     C.tab             <- table(C)
#     C.tab.vals        <- as.numeric(names(C.tab))
#     C.tab.frqs        <- as.vector(C.tab)
#     Ls.mat            <- sapply(1:length(C.tab), function(c)
#     {
#       c.use  <- C.tab.vals[c]
#       it.use <- (1:I)[which(C == c.use)]
#       Ls(data.condensed[, c(it.use, I + 1), drop = FALSE], rep(NA, C.tab.frqs[c]), delta.old[it.use], nodes, 
#          taus.old[it.use, , drop = FALSE][, (C.max + 1 - c.use):(C.max + 1 + c.use)], c.use, model)
#     })
#     dim(Ls.mat)       <- c(S, N.nodes, length(C.tab))
#     Ls.mat            <- apply(Ls.mat, 1:2, prod)                                                 # S x N.nodes
    Ls.arr            <- array(rep(Ls.mat, I * (C.max + 1)), dim = c(S, N.nodes, I, C.max + 1))   # S x N.nodes x I x (C+1)
    # 
    P.tilde.s.arr     <- array(rep(P.tilde.s.vec(Ls.mat, weights), N.nodes * I * (C.max + 1)), dim = c(S, N.nodes, I, C.max + 1)) # S x N.nodes x I x (C+1)
    r.bar.izf         <- apply(H.siz * rs.arr * Ls.arr * weights.arr / P.tilde.s.arr, 2:4, sum)             # N.nodes x I x (C+1)
    r.bar.izf.taus    <- array(rep(r.bar.izf, C.max), dim = c(N.nodes, I, C.max + 1, C.max))                                   # N.nodes x I x (C+1) x C
    N.bar.if          <- apply(r.bar.izf, 1:2, sum)
    N.bar.if.arr      <- array(rep(N.bar.if, (C.max + 1)), c(N.nodes, I, C.max + 1))                                         # N.nodes x I x (C+1)
    N.bar.if.arr.taus <- array(rep(N.bar.if, (C.max + 1) *  C.max *  C.max), c(N.nodes, I, C.max + 1,  C.max,  C.max))                           # N.nodes x I x (C+1) x C x C
    #
    
    while ((iter.inner <= (max.inner - 1)) && (max(curr.tol) > tol))
    {
      if (verbose == "yes")
      {
        cat("\r", "|", rep("-", iter.inner+1), rep(" ", max.inner - iter.inner-1), "|", sep = "")
        flush.console()
      }
      
      # M stage, part 1 of 2:
      #     Update taus, for fixed alphas and deltas.
      P.izf.arr      <- P.izf(alpha.old, delta.old, taus.old, nodes, C)                                     # N.nodes x I x (C+1)
      P.izf.arr.taus <- array(rep(P.izf.arr, C.max), dim = c(N.nodes, I, C.max+1, C.max))                                    # N.nodes x I x (C+1) x C
      
      dP             <- dP.phi(alpha.old, delta.old, taus.old, nodes, C, param = "taus")
      D1             <- DlogL.dphi (param = "taus", dP, r.bar.izf.taus, P.izf.arr.taus)
      DlogL.taus     <- D1$taus
      # 
      P.izf.arr.taus.taus <- array(rep(P.izf.arr.taus, C.max), dim = c(N.nodes, I, C.max+1, C.max, C.max))                       # N.nodes x I x (C+1) x C x C
      dP.taus.taus <- array(NA, c(N.nodes, I, C.max+1, C.max, C.max))                                                        # N.nodes x I x (C+1) x C x C
      for (f in 1:N.nodes) {
        for (i in 1:I) {
          for (z in 0:C.max) {
            dP.taus.taus[f, i, z+1, , ] <- dP$taus[f, i, z+1, ] %*% t(dP$taus[f, i, z+1, ])
          }
        }
      }
      rm(f, i, z)
      Inf.arr <- apply(N.bar.if.arr.taus * dP.taus.taus / P.izf.arr.taus.taus, c(2, 4, 5), sum, na.rm = TRUE) #!#
      # 
      taus.new <- matrix(0, nrow = I, ncol = C.max) 
      for (i in 1:I)
      {
        c.use  <- C[i]
        taus.new[i, (C.max - c.use + 1):C.max] <- 
          c(taus.old[i, (C.max - c.use + 1):C.max]) + c(solve(Inf.arr[i, 1:c.use, 1:c.use]) %*% DlogL.taus[i, 1:c.use])
      }
      # Extra control (needed in weird cases):
      taus.new <- -abs(taus.new)
      taus.new[taus.new < -10]  <- -10
      # 
      taus.new <- cbind(taus.new, 0, -taus.new[, C.max:1])
      tol.taus <- max(abs(taus.old[, 1:C.max] - taus.new[, 1:C.max]))
      taus.old <- taus.new
      
      # tol.taus <- tol
      
      # M stage, part 2 of 2:
      #     Update alphas and deltas, for fixed taus.
      P.izf.arr        <- P.izf(alpha.old, delta.old, taus.old, nodes, C)                     # N.nodes x I x (C+1)
      dP               <- dP.phi(alpha.old, delta.old, taus.old, nodes, C, param = "alphadelta")
      D1               <- DlogL.dphi (param = "alphadelta", dP, r.bar.izf, P.izf.arr)
      DlogL.alphadelta <- rbind(D1$alpha, D1$delta)
      # 
      dP.alpha.alpha   <- (dP$alpha)^2                                                        # N.nodes x I x (C+1)
      dP.alpha.delta   <- dP$alpha * dP$delta                                                 # N.nodes x I x (C+1)
      dP.delta.delta   <- (dP$delta)^2                                                        # N.nodes x I x (C+1)
      Inf.arr <- array(NA, c(I, 2, 2))
      Inf.arr[, 1, 1] <- apply(N.bar.if.arr * dP.alpha.alpha / P.izf.arr, 2, sum, na.rm = TRUE)
      Inf.arr[, 1, 2] <- apply(N.bar.if.arr * dP.alpha.delta / P.izf.arr, 2, sum, na.rm = TRUE)
      Inf.arr[, 2, 1] <- Inf.arr[, 1, 2]
      Inf.arr[, 2, 2] <- apply(N.bar.if.arr * dP.delta.delta / P.izf.arr, 2, sum, na.rm = TRUE)
      # 
      alphadelta.new <- matrix(NA, nrow = 2, ncol = I)
      for (i in 1:I) { 
        alphadelta.new[, i] <- c(alphadelta.old[, i] + solve(Inf.arr[i, , ]) %*% DlogL.alphadelta[, i])
      }
      # Extra control (needed in weird cases):
      alphadelta.new[1, ][alphadelta.new[1, ] < .1] <- .1 
      alphadelta.new[1, ][alphadelta.new[1, ] >  10] <-  10
      alphadelta.new[2, ][alphadelta.new[2, ] < -10] <- -10
      alphadelta.new[2, ][alphadelta.new[2, ] >  10] <-  10
      # 
      tol.alphadelta <-  max(abs(alphadelta.old - alphadelta.new))
      alphadelta.old <- alphadelta.new
      alpha.old <- alphadelta.old[1, ]
      delta.old <- alphadelta.old[2, ]
      
      # curr.tol <- c(tol.delta)
      # curr.tol <- c(tol.taus, tol.alpha, tol.delta)
      curr.tol <- c(tol.taus, tol.alphadelta)
      iter.inner <- iter.inner + 1
    }
    #
    iter.outer <- iter.outer + 1
  }
  
  # Information criteria:
  N      <- nrow(data)
  s.vec  <- data.condensed[, I+1]
  Ls.mat <- Ls(data.condensed, alpha.old, delta.old, taus.old, nodes, C)
  P.Xs   <- P.tilde.s.vec(Ls.mat, weights)
  # log.L  <- -sum(log(factorial(s.vec))) + sum(s.vec * log(P.Xs))
  log.L  <- sum(s.vec * log(P.Xs))
  k      <- I + I + sum(C)
  AIC    <- -2 * log.L + k * 2
  BIC    <- -2 * log.L + k * log(N)
  CAIC   <- -2 * log.L + k * (log(N) + 1)
  Inf.df <- data.frame(log.L, N.param = k, AIC, BIC, CAIC)
  
  return(list(alpha = round(alpha.old, precision), delta = round(delta.old, precision), taus = round(taus.old, precision), 
              N.nodes = N.nodes, tol = max(curr.tol), iter.inner = iter.inner, model = "GGUM", InformationCrit = Inf.df))
}