# GPCM (base for Models 3, 4, 7, and 8) ----
#    y     : Either scalar (whose value is then replicated I times), or vector of length I
#    alpha : Vector of length I
#    delta : Vector of length I
#    taus  : Either vector of length M (which is then replicated I times), or (generalized) matrix I x max(M)
#    theta : Vector of length N
#    M     : Either scalar (whose value is then replicated I times), or vector of length I
# 
# GPCM applies to Model 8 in its most general form.
#    When all alphas are equal to 1, we get Model 4.
#    When taus are constant across items (so taus is essentially a vector and M a scalar), we get Model 7.
#    Finally, when taus are constant across items (so taus is essentially a vector and M a scalar), 
#      and all alphas are equal to 1, we get Model 7.
# So, to decide which model to use, one needs to define carefully alpha (all 1 versus not all 1) 
#    and taus (one vector for all items versus a matrix of varying taus across items).

P.GPCM <- function(y, alpha, delta, taus, theta, M)
{
  N         <- length(theta)
  I         <- length(delta)
  if (length(y) == 1) {y <- rep(y, I)}
  if (is.vector(taus)) {taus <- matrix(rep(taus, I), nrow = I, byrow = TRUE)}
  if (length(M) == 1) {M <- rep(M, I)}
  taus.zero <- cbind(0, taus)
  taus.cum  <- t(apply(taus.zero, 1, cumsum))
  part      <- function(i, w) 
  {
    if ((0 <= w) && (w <= M[i])) {
      exp(alpha[i] * (w * (theta - delta[i]) - taus.cum[i, ((max(M) - M[i])/2) + w + 1]))
    } else {rep(0, N)}
  }
  num       <- sapply(1:I,      function(i) part(i, y[i]), simplify = "array")
  tmp       <- sapply(0:max(M), function(w) sapply(1:I, function(i) part(i, w)))
  den       <- matrix(rowSums(tmp, na.rm = TRUE), ncol = I, byrow = FALSE)
  return(num / den)
}

P.Model.3478 <- function(z, alpha, delta, taus, theta, C)
{
  N         <- length(theta)
  I         <- length(delta)
  if (length(z) == 1) {z <- rep(z, I)}
  if (is.vector(taus)) {taus <- matrix(rep(taus, I), nrow = I, byrow = TRUE)}
  if (length(C) == 1) {C <- rep(C, I)}
  M         <- 2 * C + 1
  mat.ind   <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
  return( mat.ind * (P.GPCM(z, alpha, delta, taus, theta, M) + P.GPCM(M - z, alpha, delta, taus, theta, M)) )
}

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
      arg            <- alpha[i] * (theta - betas[i, c])
      res.cum[, i, c] <- exp(arg) / (1 + exp(arg))
    }
  }
  res.cum <- array(c(matrix(1, N, I), res.cum, matrix(0, N, I)), dim = c(N, I, C + 2))
  res <- array(NA, c(N, I, C + 1))
  for (c in 1:(C + 1))
  {
    res[, , c] <- res.cum[, , c] - res.cum[, , c + 1]
  }
  return(res)
}
















# 
# 
# # GPCM (base for Models 4 and 8) ----
# P.GPCM <- function(y, alpha, delta, theta, taus, M)
# {
#   N         <- length(theta)
#   I         <- length(alpha)
#   if (length(y) == 1) {
#     y <- rep(y, I)
#     M <- rep(M, I)
#   }
#   if (is.vector(taus)) {taus <- matrix(rep(taus, I), nrow = I, byrow = TRUE)}
#   taus.zero <- cbind(0, taus)
#   taus.cum  <- t(apply(taus.zero, 1, cumsum))
#   part      <- function(i, w) 
#   {
#     if ((0 <= w) && (w <= M[i])) {
#       exp(alpha[i] * (w * (theta - delta[i]) - taus.cum[i, ((max(M) - M[i])/2) + w + 1]))
#     } else {rep(0, N)}
#   }
#   num       <- matrix(sapply(1:I, function(i) part(i, y[i])), nrow = N, byrow = FALSE)
#   tmp       <- sapply(0:max(M), function(w) matrix(sapply(1:I, function(i) part(i, w)), nrow = N, byrow = FALSE))
#   den       <- matrix(rowSums(tmp, na.rm = TRUE), ncol = I, byrow = FALSE)
#   return(num / den)
# }
# 
# # Model4 ----
# P.Model4 <- function(z, delta, theta, taus, C)
# {
#   I       <- length(delta)
#   M       <- 2 * C + 1
#   if (length(C) == 1) {C <- rep(C, I)}
#   N       <- length(theta)
#   mat.ind <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
#   I       <- length(delta)
#   return( mat.ind * (P.GPCM(z, rep(1, I), delta, theta, taus, M) + P.GPCM(M - z, rep(1, I), delta, theta, taus, M)) )
# }
# 
# # Model8 ----
# P.Model8 <- function(z, alpha, delta, theta, taus, C)
# {
#   if (length(C) == 1) {C <- rep(C, I)}
#   M       <- 2 * C + 1
#   N       <- length(theta)
#   mat.ind <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
#   return( mat.ind * (P.GPCM(z, alpha, delta, theta, taus, M) + P.GPCM(M - z, alpha, delta, theta, taus, M)) )
# }
# 
# #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# 
# # GRSM (base for Models 3 and 7) ----
# P.GRSM <- function(y, alpha, delta, theta, taus, M)
# {
#   N         <- length(theta)
#   I         <- length(alpha)
#   if (length(y) == 1) {
#     y <- rep(y, I)
#     M <- rep(M, I)
#   }
#   if (is.vector(taus)) {taus <- matrix(rep(taus, I), nrow = I, byrow = TRUE)}
#   taus.zero <- cbind(0, taus)
#   taus.cum  <- t(apply(taus.zero, 1, cumsum))
#   part      <- function(i, w) 
#   {
#     if ((0 <= w) && (w <= M[i])) {
#       exp(alpha[i] * (w * (theta - delta[i]) - taus.cum[i, ((max(M) - M[i])/2) + w + 1]))
#     } else {rep(0, N)}
#   }
#   num       <- matrix(sapply(1:I, function(i) part(i, y[i])), nrow = N, byrow = FALSE)
#   tmp       <- sapply(0:max(M), function(w) matrix(sapply(1:I, function(i) part(i, w)), nrow = N, byrow = FALSE))
#   den       <- matrix(rowSums(tmp, na.rm = TRUE), ncol = I, byrow = FALSE)
#   return(num / den)
# }
# 
# # Model3 ----
# P.Model3 <- function(z, delta, theta, taus, C)
# {
#   I       <- length(delta)
#   M       <- 2 * C + 1
#   if (length(C) == 1) {C <- rep(C, I)}
#   N       <- length(theta)
#   mat.ind <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
#   return( mat.ind * (P.GRSM(z, rep(1, I), delta, theta, taus, M) + P.GRSM(M - z, rep(1, I), delta, theta, taus, M)) )
# }
# 
# # # Model7 ----
# # P.Model7 <- function(z, alpha, delta, theta, taus, C)
# # {
# #   if (length(C) == 1) {C <- rep(C, I)}
# #   M       <- 2 * C + 1
# #   N       <- length(theta)
# #   mat.ind <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
# #   return( mat.ind * (P.GRSM(z, alpha, delta, theta, taus, M) + P.GRSM(M - z, alpha, delta, theta, taus, M)) )
# # }
# # 
# # #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# # #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# # #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# # 
# # # Multiple unit model (base for Models 2 and 6) ----
# # P.GMultipleUnit <- function(y, alpha, delta, theta, lambda, M)
# # {
# #   N        <- length(theta)
# #   I        <- length(alpha)
# #   if (length(y) == 1) {y <- rep(y, I)}
# #   if (length(M) == 1) {M <- rep(M, I)}
# #   part     <- function(i, w) 
# #   {
# #     if ((0 <= w) && (w <= M[i])) {
# #       exp(alpha[i] * (w * (theta - delta[i]) + w * (M[i] - w) * lambda[i]))
# #     } else {rep(0, N)}
# #   }
# #   num      <- matrix(sapply(1:I, function(i) part(i, y[i])), nrow = N, byrow = FALSE)
# #   tmp      <- sapply(0:max(M), function(w) matrix(sapply(1:I, function(i) part(i, w)), nrow = N, byrow = FALSE))
# #   if (I == 1) {tmp <- matrix(tmp, nrow = 1)}
# #   den      <- matrix(rowSums(tmp, na.rm = TRUE), ncol = I, byrow = FALSE)
# #   return(num / den)
# # }
# # 
# # # Model2 ----
# # P.Model2 <- function(z, delta, theta, lambda, C)
# # {
# #   if (length(C) == 1) {C <- rep(C, I)}
# #   M       <- 2 * C + 1
# #   N       <- length(theta)
# #   mat.ind <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
# #   I       <- length(delta)
# #   return( mat.ind * (P.GMultipleUnit(z, rep(1, I), delta, theta, lambda, M) + P.GMultipleUnit(M - z, rep(1, I), delta, theta, lambda, M)) )
# # }
# # 
# # # Model6 ----
# # P.Model6 <- function(z, alpha, delta, theta, lambda, C)
# # {
# #   if (length(C) == 1) {C <- rep(C, I)}
# #   M       <- 2 * C + 1
# #   N       <- length(theta)
# #   mat.ind <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
# #   return( mat.ind * (P.GMultipleUnit(z, alpha, delta, theta, lambda, M) + P.GMultipleUnit(M - z, alpha, delta, theta, lambda, M)) )
# # }
# # 
# # #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# # #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# # #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# # 
# # # Constant unit model (base for Models 1 and 5) ----
# # P.GConstantUnit <- function(y, alpha, delta, theta, lambda, M)
# # {
# #   N        <- length(theta)
# #   I        <- length(alpha)
# #   if (length(y) == 1) {y <- rep(y, I)}
# #   part     <- function(i, w) 
# #   {
# #     if ((0 <= w) && (w <= M)) {
# #       exp(alpha[i] * (w * (theta - delta[i]) + w * (M - w) * lambda))
# #     } else {rep(0, N)}
# #   }
# #   num      <- matrix(sapply(1:I, function(i) part(i, y)), nrow = N, byrow = FALSE)
# #   tmp      <- sapply(0:M, function(w) matrix(sapply(1:I, function(i) part(i, w)), nrow = N, byrow = FALSE))
# #   if (I == 1) {tmp <- matrix(tmp, nrow = 1)}
# #   den      <- matrix(rowSums(tmp, na.rm = TRUE), ncol = I, byrow = FALSE)
# #   return(num / den)
# # }
# # 
# # # Model1 ----
# # P.Model1 <- function(z, delta, theta, lambda, C)
# # {
# #   if (length(C) == 1) {C <- rep(C, I)}
# #   M       <- 2 * C + 1
# #   N       <- length(theta)
# #   mat.ind <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
# #   I       <- length(delta)
# #   return( mat.ind * (P.GConstantUnit(z, rep(1, I), delta, theta, lambda, M) + P.GConstantUnit(M - z, rep(1, I), delta, theta, lambda, M)) )
# # }
# # 
# # # Model5 ----
# # P.Model5 <- function(z, alpha, delta, theta, lambda, C)
# # {
# #   if (length(C) == 1) {C <- rep(C, I)}
# #   M       <- 2 * C + 1
# #   N       <- length(theta)
# #   mat.ind <- matrix(rep(z <= C, N), nrow = N, byrow = TRUE)
# #   return( mat.ind * (P.GConstantUnit(z, alpha, delta, theta, lambda, M) + P.GConstantUnit(M - z, alpha, delta, theta, lambda, M)) )
# # }
