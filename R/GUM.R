#' @title Fit the graded unfolding model (GUM)
#'   
#' @description \code{GUM} estimates all item parameters for the GUM.
#'   
#' @param data The \eqn{N\times I}{NxI} data matrix. The item scores are coded 
#'   \eqn{0, 1, \ldots, C}{0, 1, ..., C} for an item with \eqn{(C+1)} observable
#'   response categories.
#' @param C \eqn{C} is the number of observable response categories minus 1
#'   (i.e., the item scores will be in the set \eqn{\{0, 1, ..., C\}}{{0, 1,
#'   ..., C}}). It should be a scalar since the GUM expects all items to be
#'   based on the same number of observable response categories.
#' @param SE Logical value: Estimate the standard errors of the item parameter 
#'   estimates? Default is \code{TRUE}.
#' @param precision Number of decimal places of the results (default = 4).
#' @param N.nodes Number of nodes for numerical integration (default = 30).
#' @param max.outer Maximum number of outer iterations (default = 60).
#' @param max.inner Maximum number of inner iterations (default = 60).
#' @param tol Convergence tolerance (default = .001).
#'   
#' @return The function returns a list (an object of class \code{GGUM}) with 12
#'   elements: \item{data}{Data matrix.} \item{C}{Vector \eqn{C}.} 
#'   \item{alpha}{In case of the GUM this is simply a vector of 1s.} 
#'   \item{delta}{The estimated difficulty parameters.} \item{taus}{The
#'   estimated threshold parameters.} \item{SE}{The standard errors of the item
#'   parameters estimates.} \item{rows.rm}{Indices of rows removed from the data
#'   before fitting the model, due to complete disagreement.} 
#'   \item{N.nodes}{Number of nodes for numerical integration.} 
#'   \item{tol.conv}{Loss function value at convergence (it is smaller than 
#'   \code{tol} upon convergence).} \item{iter.inner}{Number of inner iterations
#'   (it is equal to 1 upon convergence).} \item{model}{Model fitted.} 
#'   \item{InformationCrit}{Loglikelihood, number of model parameters, AIC, BIC,
#'   CAIC.}
#'   
#' @section Details: The graded unfolding model (GUM; Roberts & Laughlin, 1996)
#'   is a constrained version of the GGUM (Roberts et al., 2000; see
#'   \code{\link[GGUM]{GGUM}}). GUM is constrained in two ways: All
#'   discrimination parameters are fixed to unity and the threshold parameters
#'   are shared across items. In particular, the last constraint implies that
#'   only data with the same response categories across items should be used
#'   (i.e., \eqn{C} is constant for all items).
#'   
#'   Estimated GUM parameters are used as the second step of fitting the more 
#'   general GGUM. Since under the GGUM data may include items with different 
#'   number of response categories, the code to fitting the GUM was internally 
#'   extended to accommodate for this.
#'   
#'   The marginal maximum likelihood algorithm of Roberts et al. (2000) was 
#'   implemented.
#'   
#' @references \insertRef{RobertsLaughlin1996}{GGUM}
#' 
#' \insertRef{Robertsetal2000}{GGUM}
#' 
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' # Generate data:
#' gen <- GenData.GGUM(400, 5, 3, "GUM", seed = 139)
#' # Fit the GUM:
#' fit <- GUM(gen$data, 3)
#' # Compare true and estimated item parameters:
#' cbind(gen$delta, fit$delta)
#' cbind(c(gen$taus[, 5:7]), c(fit$taus[,  5:7]))
#' 
#' @export
GUM <- function(data, C, SE = TRUE, precision = 4, 
                N.nodes = 30, max.outer = 60, max.inner = 60, tol = .001)
{
  I <- ncol(data)
  
  # Sanity check - data:
  Sanity.data(data)
  # Sanity check - C:
  Sanity.C(C, I)
  
  # Discard response patterns due to complete disagreement;
  # ceiling(C/2) gives for each item either the first agree category
  # or the middle/neutral category (if one exists) --
  # if a row has all responses as either missing or disagree,
  # we want to remove that row.
  # Note we use colSums(t(data) ...) rather than rowSums(data ...)
  # for the first comparison so that vectorization over C
  # coincides with the columns of the matrix
  rows.rm <- which(colSums(t(data) < ceiling(C/2), na.rm = TRUE)
                   + rowSums(is.na(data)) == I)
  # If there are any such rows, we need to remove them;
  # if there are no such rows, it's imperative we don't try to --
  # if we do it will remove all such rows
  data.sv <- data
  if ( length(rows.rm) > 0 ) {
    data    <- data[-rows.rm, ]
  }
  
  tmp            <- GGUM.data.condense(data)
  data.condensed <- tmp$data.condensed
  S              <- nrow(data.condensed)
  rm(tmp)
  
  # Sanity check - Cfixed:
  Sanity.Cfixed(C)
  
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
    cat("\r", "|", rep("-", iter.outer+1), 
        rep(" ", max.outer - iter.outer-1), "|", sep = "")
    flush.console()
    
    
    iter.inner <- 0
    curr.tol   <- 1
    
    # E stage: Compute r.bar.izf and N.bar.if, given the current delta, and taus.
    # Ls:
    Ls.mat <- Ls(data.condensed, alpha, delta.old, taus.old, nodes, C)
    Ls.arr <- array(rep(Ls.mat, I * (C.max + 1)), 
                    dim = c(S, N.nodes, I, C.max + 1))
    # 
    P.tilde.s.arr     <- array(rep(P.tilde.s.vec(Ls.mat, weights), 
                                   N.nodes * I * (C.max + 1)), 
                               dim = c(S, N.nodes, I, C.max + 1))
    r.bar.izf         <- apply(H.siz * rs.arr * Ls.arr * weights.arr / P.tilde.s.arr, 2:4, sum)
    r.bar.izf.taus    <- array(rep(r.bar.izf, C.max), 
                               dim = c(N.nodes, I, C.max + 1, C.max))
    N.bar.if          <- apply(r.bar.izf, 1:2, sum)
    N.bar.if.arr      <- array(rep(N.bar.if, (C.max + 1)), 
                               c(N.nodes, I, C.max + 1))
    N.bar.if.arr.taus <- array(rep(N.bar.if, (C.max + 1) *  C.max * C.max), 
                               c(N.nodes, I, C.max + 1, C.max, C.max))
    #
    
    while ((iter.inner <= (max.inner - 1)) && (max(curr.tol) > tol))
    {
      
      # M stage, part 1 of 2:
      #     Update taus, for fixed deltas.
      # P.izf.arr:
      P.izf.arr      <- P.izf(alpha, delta.old, taus.old, nodes, C)
      P.izf.arr.taus <- array(rep(P.izf.arr, C.max), 
                              dim = c(N.nodes, I, C.max + 1, C.max))
      # dP:
      dP             <- dP.phi(alpha, delta.old, taus.old, nodes, C, 
                               param = "taus")
      D1             <- DlogL.dphi(param = "taus", dP, r.bar.izf.taus, 
                                   P.izf.arr.taus)
      
      # Dirty fix for R 4.1:
      tmp.dims   <- dim(D1$taus)
      DlogL.taus <- colSums(matrix(unlist(D1$taus), tmp.dims))
      
      # 
      P.izf.arr.taus.taus <- array(rep(P.izf.arr.taus, C.max), 
                                   dim = c(N.nodes, I, C.max + 1, C.max, C.max))
      dP.taus.taus        <- array(NA, c(N.nodes, I, C.max + 1, C.max, C.max))
      for (f in 1:N.nodes) {
        for (i in 1:I) {
          for (z in 0:C.max) {
            dP.taus.taus[f, i, z+1, , ] <- dP$taus[f, i, z+1, ] %*% t(dP$taus[f, i, z+1, ])
          }
        }
      }
      rm(f, i, z)
      Inf.arr <- apply(N.bar.if.arr.taus * dP.taus.taus / P.izf.arr.taus.taus, 
                       c(4, 5), sum, na.rm = TRUE)
      # 
      taus.new <- matrix(0, nrow = I, ncol = C.max)
      for (i in 1:I)
      {
        c.use  <- C[i]
        taus.new[i, (C.max - c.use + 1):C.max] <- 
          c(taus.old[i, (C.max - c.use + 1):C.max]) + 
          c(solve(Inf.arr[1:c.use, 1:c.use]) %*% DlogL.taus[1:c.use])
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
      dP             <- dP.phi(alpha, delta.old, taus.old, nodes, C, 
                               param = "delta")
      D1             <- DlogL.dphi(param = "delta", dP, r.bar.izf, P.izf.arr)
      DlogL.delta    <- D1$delta
      # 
      dP.delta.delta <- (dP$delta)^2
      Inf.arr        <- apply(N.bar.if.arr * dP.delta.delta / P.izf.arr, 2, 
                              sum, na.rm = TRUE)
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
  cat("\n\n")
  
  # Information criteria:
  N      <- nrow(data)
  s.vec  <- data.condensed[, I+1]
  Ls.mat <- Ls(data.condensed, alpha, delta.old, taus.old, nodes, C)
  P.Xs   <- P.tilde.s.vec(Ls.mat, weights)
  log.L  <- sum(s.vec * log(P.Xs))
  k      <- I + C[1]
  AIC    <- -2 * log.L + k * 2
  BIC    <- -2 * log.L + k * log(N)
  CAIC   <- -2 * log.L + k * (log(N) + 1)
  Inf.df <- data.frame(log.L, N.param = k, AIC, BIC, CAIC)
  
  # SE:
  if (SE)
  {
    Ls.arr    <- array(rep(Ls.mat, I * (C.max + 1)), 
                       dim = c(S, N.nodes, I, C.max + 1))
    P.tilde.s <- P.tilde.s.vec(Ls.mat, weights)
    P.OBS.s   <- data.condensed[, I + 1] / N
    P.izf.arr <- array(rep(P.izf(alpha, delta.old, taus.old, nodes, C), S), 
                       dim = c(N.nodes, I, C.max + 1, S))
    P.izf.arr <- aperm(P.izf.arr, c(4, 1, 2, 3))
    
    # 
    tmp      <- dP.phi(alpha, delta.old, taus.old, nodes, C, param = "delta")
    dP.delta <- aperm(array(rep(tmp$delta, S), 
                            dim = c(N.nodes, I, C.max + 1, S)), c(4, 1, 2, 3))
    tmp      <- dP.phi(alpha, delta.old, taus.old, nodes, C, param = "taus")
    dP.taus  <- aperm(array(rep(tmp$taus, S), 
                            dim = c(N.nodes, I, C.max + 1, C.max, S)), 
                      c(5, 1, 2, 3, 4))
    rm(tmp)
    
    #
    dP.tilde.s.vec.delta <- apply(Ls.arr * weights.arr * dP.delta * H.siz / 
                                    P.izf.arr, c(1, 3), sum, na.rm = TRUE)
    rm(dP.delta)
    Ls.arr.taus          <- array(rep(Ls.arr, C.max), 
                                  dim = c(S, N.nodes, I, C.max + 1, C.max))
    rm(Ls.arr)
    weights.arr.taus     <- array(rep(weights.arr, C.max), 
                                  dim = c(S, N.nodes, I, C.max + 1, C.max))
    rm(weights.arr)
    H.siz.taus           <- array(rep(H.siz, C.max), 
                                  dim = c(S, N.nodes, I, C.max + 1, C.max))
    P.izf.arr.taus       <- array(rep(P.izf.arr, C.max), 
                                  dim = c(S, N.nodes, I, C.max + 1, C.max))
    rm(P.izf.arr)
    dP.tilde.s.vec.taus  <- apply(Ls.arr.taus * weights.arr.taus * dP.taus * 
                                    H.siz.taus / P.izf.arr.taus, c(1, 5), sum, 
                                  na.rm = TRUE)
    rm(dP.taus, Ls.arr.taus, weights.arr.taus, H.siz.taus, P.izf.arr.taus)
    
    #
    Inf.mat <- matrix(0, nrow = I + C.max, ncol = I + C.max)
    for (i1 in 1:I) {
      dP.tilde      <- dP.tilde.s.vec.delta[, i1] * dP.tilde.s.vec.delta[, i1]
      Inf.mat[i1, i1] <- N * sum(P.OBS.s * dP.tilde / (P.tilde.s^2))
    }
    for (i1 in 1:C.max) {
      for (i2 in i1:C.max) {
        dP.tilde        <- dP.tilde.s.vec.taus[, i1, drop = FALSE] * 
          dP.tilde.s.vec.taus[, i2, drop = FALSE]
        Inf.mat[I + i1, I + i2] <- N * sum(P.OBS.s * dP.tilde / 
                                             (P.tilde.s^2))
      }
    }
    for (i1 in 1:I) {
      for (i2 in 1:C.max) {
        dP.tilde        <- dP.tilde.s.vec.delta[, i1] * 
          dP.tilde.s.vec.taus[, i2, drop = FALSE]
        Inf.mat[i1, I + i2] <- N * sum(P.OBS.s * dP.tilde / (P.tilde.s^2))
      }
    }
    Inf.mat <- Inf.mat + t(Inf.mat) - diag(diag(Inf.mat))
    SE.mat  <- matrix(sqrt(diag(solve(Inf.mat))), ncol = 1)
    rownames(SE.mat) <- c(paste0("SE.delta", 1:I), paste0("SE.tau", C.max:1))
    
    SE.out <- list(
      Delta = round(SE.mat[1:I, ], precision), 
      Taus = round(SE.mat[(I + C.max):(I+1), ], precision))
  } else {SE.out <- NULL}
  
  res <- list(
    data            = data.sv, 
    C               = C, 
    alpha           = alpha, 
    delta           = round(delta.old, precision), 
    taus            = round(taus.old, precision), 
    SE              = SE.out, 
    rows.rm         = rows.rm, 
    N.nodes         = N.nodes, 
    tol.conv        = max(curr.tol), 
    iter.inner      = iter.inner, 
    model           = "GUM", 
    InformationCrit = Inf.df)
  class(res) <- "GGUM"
  return(res)
}
