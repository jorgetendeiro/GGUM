# GPCM (base for GUM and GGUM) ----
#    y     : Either scalar (whose value is then replicated I times), or vector 
#            of length I
#    alpha : Vector of length I
#    delta : Vector of length I
#    taus  : Either vector of length M (which is then replicated I times), or 
#           (generalized) matrix I x max(M)
#    theta : Vector of length N
#    M     : Either scalar (whose value is then replicated I times), or vector 
#            of length I
# 
# GPCM applies to GGUM in its most general form.

# P.GPCM ----
P.GPCM <- function(y, alpha, delta, taus, theta, M)
{
  N         <- length(theta)
  I         <- length(delta)
  if (length(y) == 1) y <- rep(y, I)
  if (is.vector(taus)) taus <- matrix(rep(taus, I), nrow = I, byrow = TRUE)
  if (length(M) == 1) M <- rep(M, I)
  taus.zero <- cbind(0, taus)
  taus.cum  <- t(apply(taus.zero, 1, cumsum))
  part      <- function(i, w) 
  {
    if ((0 <= w) && (w <= M[i])) {
      exp(alpha[i] * (w * (theta - delta[i]) - taus.cum[i, ((max(M) - M[i])/2) + w + 1]))
    } else rep(0, N)
  }
  num       <- sapply(1:I,      function(i) part(i, y[i]), simplify = "array")
  tmp       <- sapply(0:max(M), function(w) sapply(1:I, function(i) part(i, w)))
  den       <- matrix(rowSums(tmp, na.rm = TRUE), ncol = I, byrow = FALSE)
  return(num / den)
}

# P.GGUM ----
P.GGUM <- function(z, alpha, delta, taus, theta, C)
{
  N         <- length(theta)
  I         <- length(delta)
  if (length(z) == 1) z <- rep(z, I)
  if (is.vector(taus)) taus <- matrix(rep(taus, I), nrow = I, byrow = TRUE)
  if (length(C) == 1) C <- rep(C, I)
  M         <- 2 * C + 1
  mat.ind   <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
  return( mat.ind * (P.GPCM(z, alpha, delta, taus, theta, M) + P.GPCM(M - z, alpha, delta, taus, theta, M)) )
}

# probs.GGUM ----
#' @title Compute model probabilities for the GGUM
#'   
#' @description \code{probs.GGUM} computes model probabilities for the GGUM (and
#'   the GUM) for given item and person parameters.
#'   
#' @param alpha A vector of length \eqn{I} with the discrimination parameters.
#' @param delta A vector of length \eqn{I} with the difficulty parameters.
#' @param taus An \eqn{I\times M}{IxM} matrix with the threshold parameters 
#'   (\eqn{M = 2\times\max{C}+1}{M = 2*max(C)+1}).
#' @param theta A vector of length \eqn{N} with the person parameters.
#' @param C \eqn{C} is the number of observable response categories minus 1
#'   (i.e., the item scores will be in the set \eqn{\{0, 1, ..., C\}}). It
#'   should either be a vector of \eqn{I} elements or a scalar. In the latter
#'   case, it is assumed that \eqn{C} applies to all items.
#'   
#' @return The function returns an \eqn{N\times I\times K}{NxIxK} array with the
#'   GGUM probabilities, with \eqn{K=\max{C}+1}{K=max(C)+1}. To retrieve the
#'   GUM-based probabilities just constrain alpha to a unit vector of length {I}
#'   (i.e., \code{alpha = rep(1, I)}). In this case, make sure \code{C} is
#'   constant across items.
#'   
#' @section Details: This function computes the GGUM-based probabilities for all
#'   (person, item, response category) combinations. For the GGUM formula see
#'   the help for function \code{GGUM} (\code{\link[GGUM]{GGUM}}).
#'   
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' C <- c(3, 3, 3, 5, 5)
#' gen <- GenData.GGUM(10, 5, C, seed = 456)
#' gen.alpha <- gen$alpha.gen
#' gen.delta <- gen$delta.gen
#' gen.taus  <- gen$taus.gen
#' gen.theta <- gen$theta.gen
#'  
#' # Compute model probabilities for the parameters above:
#' Ps <- probs.GGUM(gen.alpha, gen.delta, gen.taus, gen.theta, C)
#' Ps
#' # In particular, the sum of the probabilities across all response options 
#' # (i.e., the third dimension) should be 1 for all (person, item) combinations:
#' apply(Ps, 1:2, sum)
#' @export
probs.GGUM <- function(alpha, delta, taus, theta, C)
{
  # Sanity check - parameters:
  Sanity.params(alpha, delta, taus, theta, C)
  
  N     <- length(theta)
  I     <- length(alpha)
  C.max <- max(C)
  res <- array(0, dim = c(N, I, C.max + 1))
  for (c in 0:C.max) res[, , c+1] <- P.GGUM(c, alpha, delta, taus, theta, C)
  dimnames(res)[[1]] <- paste0("N", 1:N)
  dimnames(res)[[2]] <- paste0("I", 1:I)
  dimnames(res)[[3]] <- paste0("C=", 0:C.max)
  return(res)
}

# P.GRM ----
P.GRM <- function(C, IP, theta)
{
  N         <- length(theta)
  I         <- nrow(IP)
  alpha     <- IP[, ncol(IP)]
  betas     <- IP[, -ncol(IP)]
  res.cum   <- array(NA, c(N, I, C))
  for (i in 1:I)
  {
    for (c in 1:C)
    {
      arg              <- alpha[i] * (theta - betas[i, c])
      res.cum[ , i, c] <- exp(arg) / (1 + exp(arg))
    }
  }
  res.cum <- array(c(matrix(1, N, I), res.cum, matrix(0, N, I)), dim = c(N, I, C + 2))
  res     <- array(NA, c(N, I, C + 1))
  for (c in 1:(C + 1)) res[, , c] <- res.cum[, , c] - res.cum[, , c + 1]
  return(res)
}
