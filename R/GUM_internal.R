# This version of GUM is to get initial values for GGUM().
# It generalizes GUM() by allowing varying C values across items.
# This is a stripped-down version of GUM().
GUM.internal <- function(data, C, 
                N.nodes = 30, max.outer = 60, max.inner = 60, tol = .001)
{
  I <- ncol(data)
  
  tmp            <- GGUM.data.condense(data)
  data.condensed <- tmp$data.condensed
  S              <- nrow(data.condensed)
  rm(tmp)
  
  # Because models 4 and GGUM (based on possibly different C per item) need 
  # function GUM() to estimate starting values, the code of GUM() needs to be 
  # versatile and allow for a vector of varying C values.
  # Observe that GUM assumes C fixed for all items.
  # Below, we therefore accommodate the code for vector of C values for GGUM 
  # purposes only.
  
  if (length(C) == 1) {C <- rep(C, I)}
  C.max          <- max(C)
  M              <- 2 * C + 1
  
  # Initial values derived using the procedure in the appendix of Roberts and 
  # Laughlin (1996):
  Initprm        <- GUM.initprm(data.condensed, C)
  delta.old      <- Initprm$delta.ini
  taus.old       <- Initprm$taus.ini
  rm(Initprm)
  alpha          <- rep(1, I) # Constant throughout
  
  # Nodes and weights:
  nodes   <- seq(-4, 4, length.out = N.nodes)
  weights <- dnorm(nodes) / sum(dnorm(nodes))
  # 
  rs.arr      <- array(rep(data.condensed[, I + 1], N.nodes * I * (C.max + 1)), 
                       dim = c(S, N.nodes, I, C.max + 1))
  weights.arr <- array(rep(weights, S * I * (C.max + 1)), 
                       dim = c(N.nodes, S, I, C.max + 1))
  weights.arr <- aperm(weights.arr, c(2, 1, 3, 4))
  H.siz       <- array(0, dim = c(S, I, C.max + 1))
  for (s in 1:S) {
    for (i in 1:I) {
      H.siz[s, i, data.condensed[s, i] + 1] <- 1
    }
  }
  rm(s, i)
  H.siz  <- array(rep(H.siz, N.nodes), dim = c(S, I, C.max + 1, N.nodes))
  H.siz  <- aperm(H.siz, c(1, 4, 2, 3))
  # 
  
  iter.outer <- 0
  iter.inner <- 10 # just to pass the first 'while' test below.
  
  while ((iter.outer <= (max.outer - 1)) && (iter.inner > 1))
  {
    cat("\r", "|", rep("-", iter.outer+1), rep(" ", max.outer - iter.outer-1), "|", sep = "")
    flush.console()
    
    
    iter.inner <- 0
    curr.tol   <- 1
    
    # E stage: Compute r.bar.izf and N.bar.if, given the current delta, and taus.
    # Ls:
    Ls.mat <- Ls(data.condensed, alpha, delta.old, taus.old, nodes, C)
    Ls.arr <- array(rep(Ls.mat, I * (C.max + 1)), dim = c(S, N.nodes, I, C.max + 1))
    # 
    P.tilde.s.arr     <- array(rep(P.tilde.s.vec(Ls.mat, weights), N.nodes * I * (C.max + 1)), dim = c(S, N.nodes, I, C.max + 1))
    r.bar.izf         <- apply(H.siz * rs.arr * Ls.arr * weights.arr / P.tilde.s.arr, 2:4, sum)
    r.bar.izf.taus    <- array(rep(r.bar.izf, C.max), dim = c(N.nodes, I, C.max + 1, C.max))
    N.bar.if          <- apply(r.bar.izf, 1:2, sum)
    N.bar.if.arr      <- array(rep(N.bar.if, (C.max + 1)), c(N.nodes, I, C.max + 1))
    N.bar.if.arr.taus <- array(rep(N.bar.if, (C.max + 1) *  C.max * C.max), c(N.nodes, I, C.max + 1, C.max, C.max))
    #
    
    while ((iter.inner <= (max.inner - 1)) && (max(curr.tol) > tol))
    {
      
      # M stage, part 1 of 2:
      #     Update taus, for fixed deltas.
      # P.izf.arr:
      P.izf.arr      <- P.izf(alpha, delta.old, taus.old, nodes, C)
      P.izf.arr.taus <- array(rep(P.izf.arr, C.max), dim = c(N.nodes, I, C.max + 1, C.max))
      # dP:
      dP             <- dP.phi(alpha, delta.old, taus.old, nodes, C, param = "taus")
      D1             <- DlogL.dphi(param = "taus", dP, r.bar.izf.taus, P.izf.arr.taus)
      
      # Dirty fix for R 4.1:
      tmp.dims   <- dim(D1$taus)
      DlogL.taus <- colSums(matrix(unlist(D1$taus), tmp.dims))
      
      # 
      P.izf.arr.taus.taus <- array(rep(P.izf.arr.taus, C.max), dim = c(N.nodes, I, C.max + 1, C.max, C.max))
      dP.taus.taus        <- array(NA, c(N.nodes, I, C.max + 1, C.max, C.max))
      for (f in 1:N.nodes) {
        for (i in 1:I) {
          for (z in 0:C.max) {
            dP.taus.taus[f, i, z+1, , ] <- dP$taus[f, i, z+1, ] %*% t(dP$taus[f, i, z+1, ])
          }
        }
      }
      rm(f, i, z)
      Inf.arr <- apply(N.bar.if.arr.taus * dP.taus.taus / P.izf.arr.taus.taus, c(4, 5), sum, na.rm = TRUE)
      # 
      taus.new <- matrix(0, nrow = I, ncol = C.max)
      for (i in 1:I)
      {
        c.use  <- C[i]
        taus.new[i, (C.max - c.use + 1):C.max] <- 
          c(taus.old[i, (C.max - c.use + 1):C.max]) + c(solve(Inf.arr[1:c.use, 1:c.use]) %*% DlogL.taus[1:c.use])
      }
      # Extra control (needed in weird cases):
      taus.new <- -abs(taus.new)
      taus.new[taus.new < -10]  <- -10
      # 
      taus.new <- cbind(taus.new, 0, -taus.new[, C.max:1])
      tol.taus <- max(abs(taus.old[, 1:C.max] - taus.new[, 1:C.max]))
      taus.old <- taus.new
      
      # M stage, part 2 of 2:
      #     Update deltas, for fixed taus.
      # P.izf.arr:
      P.izf.arr <- P.izf(alpha, delta.old, taus.old, nodes, C)
      # dP:
      dP             <- dP.phi(alpha, delta.old, taus.old, nodes, C, param = "delta")
      D1             <- DlogL.dphi(param = "delta", dP, r.bar.izf, P.izf.arr)
      DlogL.delta    <- D1$delta
      # 
      dP.delta.delta <- (dP$delta)^2
      Inf.arr        <- apply(N.bar.if.arr * dP.delta.delta / P.izf.arr, 2, sum, na.rm = TRUE)
      # 
      # 
      # Dirty fix for R 4.1:
      delta.new <- delta.old + (1/Inf.arr) * unlist(DlogL.delta)
      # Extra control (needed in weird cases):
      delta.new[delta.new < -10] <- -10
      delta.new[delta.new >  10] <-  10
      #
      tol.delta <-  max(abs(delta.old - delta.new))
      delta.old <- delta.new
      
      curr.tol <- c(tol.taus, tol.delta)
      iter.inner <- iter.inner + 1
    }
    #
    iter.outer <- iter.outer + 1
  }
  
  res <- list(
    delta           = delta.old, 
    taus            = taus.old)
  return(res)
}

