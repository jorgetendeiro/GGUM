#' @title Fit the generalized graded unfolding model (GGUM)
#'   
#' @description \code{GGUM} estimates all item parameters for the GGUM.
#'   
#' @param data The \eqn{N\times I}{NxI} data matrix. The item scores are coded 
#'   \eqn{0, 1, \ldots, C}{0, 1, ..., C} for an item with \eqn{(C+1)} observable
#'   response categories.
#' @param C \eqn{C} is the number of observable response categories minus 1
#'   (i.e., the item scores will be in the set \eqn{\{0, 1, ..., C\}}{{0, 1,
#'   ..., C}}). It should either be a vector of \eqn{I} elements or a scalar. In
#'   the latter case, it is assumed that \eqn{C} applies to all items.
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
#'   \item{alpha}{The estimated discrimination parameters for the GGUM.} 
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
#' @section Details: The generalized graded unfolding model (GGUM; Roberts &
#'   Laughlin, 1996; Roberts et al., 2000) is given by \deqn{P(Z_i=z|\theta_n) =
#'   \frac{f(z) + f(M-z)}{\sum_{w=0}^C\left[f(w)+f(M-w)\right]}, }{P(Z_i =
#'   z|t_n) = ( f(z) + f(M-z) ) / (sum( f(w) + f(M - w); w = 0, ..., C )),}
#'   
#'   \deqn{f(w) = exp\left\{\alpha_i\left[w(\theta_n-\delta_i)- 
#'   \sum_{k=0}^w\tau_{ik}\right]\right\}, }{f(w) = exp( alpha_i ( w(t_n -
#'   delta_i) - sum( tau_ik; k = 0, ..., w) ) ),}
#'   
#'   where: \itemize{ \item The subscripts \eqn{i} and \eqn{n} identify the item
#'   and person, respectively. \item \eqn{z=0,\ldots,C}{z = 0, ..., C} denotes
#'   the observed answer response. \item \eqn{M = 2C + 1} is the number of
#'   subjective response options minus 1. \item \eqn{\theta_n}{t_n} is the
#'   latent trait score for person \eqn{n}. \item \eqn{\alpha_i}{alpha_i} is the
#'   item slope (discrimination). \item \eqn{\delta_i}{delta_i} is the item
#'   location. \item \eqn{\tau_{ik}}{tau_ik} (\eqn{k=1,\ldots,M}{k = 1, ..., M}
#'   ) are the threshold parameters. }
#'   
#'   Parameter \eqn{\tau_{i0}}{tau_i0} is arbitrarily constrained to zero and
#'   the threshold parameters are constrained to symmetry around zero, that is, 
#'   \eqn{\tau_{i(C+1)}=0}{tau_{i(C+1)} = 0} and 
#'   \eqn{\tau_{iz}=-\tau_{i(M-z+1)}}{tau_{iz} = -tau_{i(M-z+1)}} for 
#'   \eqn{z\not= 0}{z != 0}.
#'   
#'   The marginal maximum likelihood algorithm of Roberts et al. (2000) was 
#'   implemented.
#'   
#' @references \insertRef{RobertsLaughlin1996}{GGUM}
#' 
#' \insertRef{Robertsetal2000}{GGUM}
#' 
#' @author Jorge N. Tendeiro, \email{j.n.tendeiro@rug.nl}
#'   
#' @examples
#' \dontrun{
#' # Example 1 - Same value C across items:
#' # Generate data:
#' gen1 <- GenData.GGUM(2000, 10, 2, seed = 125)
#' # Fit the GGUM:
#' fit1 <- GGUM(gen1$data, 2)
#' # Compare true and estimated item parameters:
#' cbind(gen1$alpha, fit1$alpha)
#' cbind(gen1$delta, fit1$delta)
#' cbind(c(gen1$taus[, 4:5]), c(fit1$taus[, 4:5]))
#' 
#' # Example 2 - Different C across items:
#' # Generate data:
#' set.seed(1); C <- sample(3:5, 10, replace = TRUE)
#' gen2 <- GenData.GGUM(2000, 10, C, seed = 125)
#' # Fit the GGUM:
#' fit2 <- GGUM(gen2$data, C)
#' # Compare true and estimated item parameters:
#' cbind(gen2$alpha, fit2$alpha)
#' cbind(gen2$delta, fit2$delta)
#' cbind(c(gen2$taus[, 7:11]), c(fit2$taus[, 7:11]))
#' }
#' @export
GGUM <- function(data, C, SE = TRUE, precision = 4, 
                 N.nodes = 30, max.outer = 60, max.inner = 60, tol = .001)
{
  I <- ncol(data)
  
  # Sanity check - data:
  Sanity.data(data)
  # Sanity check - C:
  Sanity.C(C, I)
  
  # Discard response patterns due to complete disagreement:
  rows.rm <- which(rowSums(data<2, na.rm = TRUE) + rowSums(is.na(data)) == I)
  data.sv <- data
  data    <- data[-rows.rm, ]
  
  tmp            <- GGUM.data.condense(data)
  data.condensed <- tmp$data.condensed
  S              <- nrow(data.condensed)
  rm(tmp)
  
  if (length(C) == 1) {C <- rep(C, I)}
  C.max           <- max(C)
  M               <- 2 * C + 1
  
  # Initial values, Step 1 of 2 - derived from GUM:
  cat("\nStep 1 of 3: Calibrating initial parameters by means of GUM... \n")
  GUM.res        <- GUM.internal(data, C)
  delta.old      <- GUM.res$delta
  taus.old       <- GUM.res$taus
  rm(GUM.res)
  # Initial values, Step 2 of 2 - derived from Model 4, based on the estimates 
  #   from GUM as initial values:
  cat("\nStep 2 of 3: Calibrating initial parameters by means of GUM extended... \n")
  Model4.res     <- Model4(data, C, 
                           Init.vals = list(delta = delta.old, taus = taus.old))
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
  
  cat("\nStep 3 of 3: Estimating item parameters for GGUM... \n")
  
  while ((iter.outer <= (max.outer - 1)) && (iter.inner > 1))
  {
    
    cat("\r", "|", rep("-", iter.outer+1), 
        rep(" ", max.outer - iter.outer-1), "|", sep = "")
    flush.console()
    
    
    iter.inner <- 0
    curr.tol   <- 1
    
    # E stage: Compute r.bar.izf and N.bar.if, 
    #             given the current alpha, delta, and taus.
    # Ls:
    Ls.mat <- Ls(data.condensed, alpha.old, delta.old, taus.old, nodes, C)
    Ls.arr            <- array(rep(Ls.mat, I * (C.max + 1)), 
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
    N.bar.if.arr.taus <- array(rep(N.bar.if, (C.max + 1) *  C.max *  C.max), 
                               c(N.nodes, I, C.max + 1,  C.max,  C.max))
    #
    
    while ((iter.inner <= (max.inner - 1)) && (max(curr.tol) > tol))
    {
      
      
      # M stage, part 1 of 2:
      #     Update taus, for fixed alphas and deltas.
      P.izf.arr      <- P.izf(alpha.old, delta.old, taus.old, nodes, C)
      P.izf.arr.taus <- array(rep(P.izf.arr, C.max), 
                              dim = c(N.nodes, I, C.max+1, C.max))
      
      dP             <- dP.phi(alpha.old, delta.old, taus.old, nodes, C, 
                               param = "taus")
      D1             <- DlogL.dphi (param = "taus", dP, r.bar.izf.taus, 
                                    P.izf.arr.taus)
      DlogL.taus     <- D1$taus
      # 
      P.izf.arr.taus.taus <- array(rep(P.izf.arr.taus, C.max), 
                                   dim = c(N.nodes, I, C.max+1, C.max, C.max))
      dP.taus.taus <- array(NA, c(N.nodes, I, C.max+1, C.max, C.max))
      for (f in 1:N.nodes) {
        for (i in 1:I) {
          for (z in 0:C.max) {
            dP.taus.taus[f, i, z+1, , ] <- dP$taus[f, i, z+1, ] %*% t(dP$taus[f, i, z+1, ])
          }
        }
      }
      rm(f, i, z)
      Inf.arr <- apply(N.bar.if.arr.taus * dP.taus.taus / P.izf.arr.taus.taus, 
                       c(2, 4, 5), sum, na.rm = TRUE)
      # 
      taus.new <- matrix(0, nrow = I, ncol = C.max) 
      for (i in 1:I)
      {
        c.use  <- C[i]
        taus.new[i, (C.max - c.use + 1):C.max] <- 
          c(taus.old[i, (C.max - c.use + 1):C.max]) + 
          c(solve(Inf.arr[i, 1:c.use, 1:c.use]) %*% DlogL.taus[i, 1:c.use])
      }
      # Extra control (needed in cases of difficult convergence):
      taus.new <- -abs(taus.new)
      taus.new[taus.new < -10]  <- -10
      # 
      taus.new <- cbind(taus.new, 0, -taus.new[, C.max:1])
      tol.taus <- max(abs(taus.old[, 1:C.max] - taus.new[, 1:C.max]))
      taus.old <- taus.new
      
      # M stage, part 2 of 2:
      #     Update alphas and deltas, for fixed taus.
      P.izf.arr        <- P.izf(alpha.old, delta.old, taus.old, nodes, C)
      dP               <- dP.phi(alpha.old, delta.old, taus.old, nodes, C, 
                                 param = "alphadelta")
      D1               <- DlogL.dphi (param = "alphadelta", dP, r.bar.izf, 
                                      P.izf.arr)
      DlogL.alphadelta <- rbind(D1$alpha, D1$delta)
      # 
      dP.alpha.alpha   <- (dP$alpha)^2
      dP.alpha.delta   <- dP$alpha * dP$delta
      dP.delta.delta   <- (dP$delta)^2
      Inf.arr <- array(NA, c(I, 2, 2))
      Inf.arr[, 1, 1] <- apply(N.bar.if.arr * dP.alpha.alpha / P.izf.arr, 2, 
                               sum, na.rm = TRUE)
      Inf.arr[, 1, 2] <- apply(N.bar.if.arr * dP.alpha.delta / P.izf.arr, 2, 
                               sum, na.rm = TRUE)
      Inf.arr[, 2, 1] <- Inf.arr[, 1, 2]
      Inf.arr[, 2, 2] <- apply(N.bar.if.arr * dP.delta.delta / P.izf.arr, 2, 
                               sum, na.rm = TRUE)
      # 
      alphadelta.new <- matrix(NA, nrow = 2, ncol = I)
      for (i in 1:I) { 
        alphadelta.new[, i] <- c(alphadelta.old[, i] + solve(Inf.arr[i, , ]) %*% DlogL.alphadelta[, i])
      }
      # Extra control (needed in cases of difficult convergence):
      alphadelta.new[1, ][alphadelta.new[1, ] < .1] <- .1 
      alphadelta.new[1, ][alphadelta.new[1, ] >  10] <-  10
      alphadelta.new[2, ][alphadelta.new[2, ] < -10] <- -10
      alphadelta.new[2, ][alphadelta.new[2, ] >  10] <-  10
      # 
      tol.alphadelta <-  max(abs(alphadelta.old - alphadelta.new))
      alphadelta.old <- alphadelta.new
      alpha.old <- alphadelta.old[1, ]
      delta.old <- alphadelta.old[2, ]
      
      curr.tol <- c(tol.taus, tol.alphadelta)
      iter.inner <- iter.inner + 1
    }
    #
    iter.outer <- iter.outer + 1
  }
  cat("\n\n")
  
  # Information criteria:
  N      <- nrow(data)
  s.vec  <- data.condensed[, I+1]
  Ls.mat <- Ls(data.condensed, alpha.old, delta.old, taus.old, nodes, C)
  P.Xs   <- P.tilde.s.vec(Ls.mat, weights)
  log.L  <- sum(s.vec * log(P.Xs))
  k      <- I + I + sum(C)
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
    P.izf.arr <- array(rep(P.izf(alpha.old, delta.old, taus.old, nodes, C), S), 
                       dim = c(N.nodes, I, C.max + 1, S))
    P.izf.arr <- aperm(P.izf.arr, c(4, 1, 2, 3))
    
    # 
    tmp      <- dP.phi(alpha.old, delta.old, taus.old, nodes, C, 
                       param = "alphadelta")
    dP.alpha <- aperm(array(rep(tmp$alpha, S), 
                            dim = c(N.nodes, I, C.max + 1, S)), c(4, 1, 2, 3))
    dP.delta <- aperm(array(rep(tmp$delta, S), 
                            dim = c(N.nodes, I, C.max + 1, S)), c(4, 1, 2, 3))
    tmp      <- dP.phi(alpha.old, delta.old, taus.old, nodes, C, 
                       param = "taus")
    dP.taus  <- tmp$taus
    rm(tmp)
    
    #
    dP.tilde.s.vec.alpha <- apply(Ls.arr * weights.arr * dP.alpha * H.siz / 
                                    P.izf.arr, c(1, 3), sum, na.rm = TRUE)
    dP.tilde.s.vec.delta <- apply(Ls.arr * weights.arr * dP.delta * H.siz / 
                                    P.izf.arr, c(1, 3), sum, na.rm = TRUE)
    rm(dP.alpha, dP.delta)
    dP.tilde.s.vec.taus <- array(NA, dim = c(S, I, C.max))
    for (i in 1:C.max)
    {
      dP.taus.use <- aperm(array(rep(dP.taus[ , , , i], S), 
                                 dim = c(N.nodes, I, C.max + 1, S)), 
                           c(4, 1, 2, 3))
      dP.tilde.s.vec.taus[, , i] <- apply(Ls.arr * weights.arr * dP.taus.use * 
                                            H.siz / P.izf.arr, c(1, 3), sum, 
                                          na.rm = TRUE)
    }
    rm(dP.taus, Ls.arr, weights.arr, H.siz, P.izf.arr, i, dP.taus.use)
    #
    SE.mat <- matrix(0, nrow = I, ncol = C.max + 2)
    for (i in 1:I) {
      dP.tilde      <- cbind(dP.tilde.s.vec.alpha[, i], 
                             dP.tilde.s.vec.delta[, i], 
                             dP.tilde.s.vec.taus[, i, ])
      P.OBS.s.mat   <- matrix(rep(P.OBS.s  , C.max + 2), nrow = S, 
                              byrow = FALSE)
      P.tilde.s.mat <- matrix(rep(P.tilde.s, C.max + 2), nrow = S, 
                              byrow = FALSE)
      sum.arg       <- sqrt(P.OBS.s.mat) * dP.tilde / P.tilde.s.mat
      Inf.i         <- N * (t(sum.arg) %*% sum.arg)
      tmp           <- sqrt(diag(solve(Inf.i[1:(C[i] + 2), 1:(C[i] + 2)])))
      SE.mat[i, 1:(C[i] + 2)]   <- c(tmp[1:2], rev(tmp[-(1:2)]))
    }
    
    colnames(SE.mat) <- c("SE.alpha", "SE.delta", paste0("SE.tau", 1:C.max))
    SE.out <- round(SE.mat, precision)
  } else {SE.out <- NULL}
  
  res <- list(
    data            = data.sv, 
    C               = C, 
    alpha           = round(alpha.old, precision), 
    delta           = round(delta.old, precision), 
    taus            = round(taus.old, precision), 
    SE              = SE.out, 
    rows.rm         = rows.rm,
    N.nodes         = N.nodes, 
    tol.conv        = max(curr.tol), 
    iter.inner      = iter.inner, 
    model           = "GGUM", 
    InformationCrit = Inf.df)
  class(res) <- "GGUM"
  return(res)
}
