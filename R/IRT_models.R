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
