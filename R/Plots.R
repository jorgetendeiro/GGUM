# Plot IRF ----
plot.GGUM <- function(C, IP.res, items = NULL, x.lim = 4, ThetaminDelta = TRUE) 
{
  if (length(C) == 1) {C <- rep(C, length(IP.res$alpha))}
  M     <- 2 * C + 1
  C.max <- max(C)
  
  I         <- length(IP.res$alpha)
  alpha     <- IP.res$alpha
  delta     <- IP.res$delta
  taus      <- IP.res$taus
  if (ThetaminDelta == TRUE) {
    th.lims   <- cbind(delta - x.lim, delta + x.lim)
    th.vals   <- t(apply(th.lims, 1, function(vec) {seq(vec[1], vec[2], 
                                                        length.out = 100)}))
    x         <- seq(-x.lim, x.lim, length.out = 100)
  } else {
    th.lims   <- t(sapply(delta, function(d) {c(-max(ceiling(d), 4), 
                                                max(4, ceiling(d)))}))
    th.vals   <- t(apply(th.lims, 1, function(vec) {seq(vec[1], vec[2], 
                                                        length.out = 100)}))
    x         <- seq(-x.lim, x.lim, length.out = 100)
  }
  #
  if (is.null(items)) {I.plot <- 1:I} else {I.plot <- items}
  for (i in 1:length(I.plot)) {
    it <- I.plot[i]
    invisible(readline(prompt="Press [enter] to continue"))
    #
    plot(x, P.GGUM(0, alpha[it], delta[it], 
                   taus[it, (C.max - C[it] + 1):(C.max - C[it] + M[it])], 
                   th.vals[it, ], C[it]),
         type = "l", lty = 1, col = plasma(C[it] + 1)[1], lwd = 2,
         xlim = c(-x.lim, x.lim), ylim = 0:1,
         xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = paste("Item", it))
    for (c in 1:C[it]) {
      points(x, P.GGUM(c, alpha[it], delta[it], 
                       taus[it, (C.max - C[it] + 1):(C.max - C[it] + M[it])], 
                       th.vals[it, ], C[it]),
             type = "l", lty = 1, col = plasma(C[it] + 1)[c + 1], lwd = 2,
             xlim = c(-x.lim, x.lim), ylim = 0:1)
    }
    axis(1, at = seq(-x.lim, x.lim, 1))
    axis(2, at = seq(0, 1, .20), las = 1)
    if (ThetaminDelta == TRUE) {mtext(expression(paste(theta, " - ", delta)), 
                                      side = 1, line = 2.5, cex = 1.2)} else {
                                        mtext(expression(paste("Item location ", delta)), side = 1, 
                                              line = 2.5, cex = 1.2)
                                      }
    mtext(expression(paste("P(X=1|", theta, ")")), side = 2, line = 2.5, 
          cex = 1.2)
    #
    legend("top", paste0("C = ", 0:C[it]), col = plasma(C[it] + 1),
           lty = 1, lwd = 2, inset = .01, cex = .8, ncol = 2)
  }
}

# Plot test characteristic curve ----
plot.TestCharacteristicCurve.GGUM <- function(data, C, IP.res, Th.res) {
  if (length(C) == 1) {C <- rep(C, length(IP.res$alpha))}
  M     <- 2 * C + 1
  C.max <- max(C)
  #
  alpha     <- IP.res$alpha
  delta     <- IP.res$delta
  taus      <- IP.res$taus
  theta     <- Th.res[[2]]
  I         <- length(IP.res$alpha)
  N         <- length(theta)
  # Test characteristic curve:
  OBS.scores      <- data
  theta.exp       <- seq(min(-4, min(theta)), max(4, max(theta)), 
                         length.out = 1000)
  N.groups        <- 100
  int.lims        <- quantile(theta.exp, probs = seq(0, 1, 1 / (N.groups - 1)))
  names(int.lims) <- NULL
  #
  scores.mat <- matrix(0, nrow = I, ncol = C.max + 1)
  for (i in 1:I) {
    scores.mat[i, 1:(C[i] + 1)] <- 0:C[i]
  }
  scores.arr <- aperm(array(rep(scores.mat, N), 
                            dim = c(I, C.max + 1, N.groups)), c(3, 1, 2))
  EXP.scores <- apply(P.izf(alpha, delta, taus, int.lims, C) * 
                        scores.arr, 1, sum)
  res        <- cbind(th =  int.lims, observed = rep(NA, N.groups), 
                      expected = EXP.scores)
  for (int in 1:N.groups) {
    pos.obs     <- which( (theta >= int.lims[int]) & (theta < int.lims[int + 1]))
    if (length(pos.obs) > 0) {
      res[int, 2] <- mean(rowSums(OBS.scores[pos.obs, , drop = FALSE], na.rm = TRUE))
    }
  }
  #
  plot(res[, 1], res[, 3], type = "l", lty = 1, lwd = 2,
       xlim = c(theta.exp[1], theta.exp[1000]),
       ylim = c(0, ceiling(max(res[, 3], na.rm = TRUE))),
       xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = "Test Characteristic Curve")
  points(res[, 1], res[, 2], type = "p", pch = 21, bg = "grey")
  axis(1, at = seq(floor(theta.exp[1]), ceiling(theta.exp[1000]), 1))
  axis(2, at = seq(0, ceiling(max(res[, 3], na.rm = TRUE)), 10), las = 1)
  mtext(expression(theta), side = 1, line = 2.5, cex = 1.2)
  mtext("Total Score", side = 2, line = 2.5, cex = 1.2)
  #
  cor.OBS.EXP.means <- cor(res[, 2], res[, 3], use = "pairwise.complete.obs")
  scores.arr.N      <- aperm(array(rep(scores.mat, N), dim = c(I, C.max + 1, N)), c(3, 1, 2))
  EXP.scores.N      <- apply(P.izf(alpha, delta, taus, theta, C) * scores.arr.N, 1:2, sum)
  cor.OBS.EXP.raw   <- cor(c(OBS.scores), c(EXP.scores.N), use = "pairwise.complete.obs")
  #
  res.out           <- res[, c(1, 3)]
  colnames(res.out) <- c("Theta", "TotalScore")
  return(list(res.out, cor.OBS.EXP.raw = round(cor.OBS.EXP.raw, 4), cor.OBS.EXP.means = round(cor.OBS.EXP.means, 4)))
}

# Plot item characteristic curves ----
plot.ItemCharacteristicCurve.GGUM <- function(data, C, IP.res, Th.res, items = NULL) {
  if (length(C) == 1) {C <- rep(C, length(IP.res$alpha))}
  M     <- 2 * C + 1
  C.max <- max(C)
  #
  alpha     <- IP.res$alpha
  delta     <- IP.res$delta
  taus      <- IP.res$taus
  theta     <- Th.res[[2]]
  I         <- length(IP.res$alpha)
  N         <- length(theta)
  # Item characteristic curves:
  OBS.scores      <- data
  theta.exp       <- seq(min(-4, min(theta)), max(4, max(theta)), length.out = 1000)
  N.groups        <- 50
  int.lims        <- quantile(theta.exp, probs = seq(0, 1, 1 / (N.groups - 1)))
  names(int.lims) <- NULL
  #
  scores.mat <- matrix(0, nrow = I, ncol = C.max + 1)
  for (i in 1:I) {
    scores.mat[i, 1:(C[i] + 1)] <- 0:C[i]
  }
  scores.arr <- aperm(array(rep(scores.mat, N.groups), dim = c(I, C.max + 1, N.groups)), c(3, 1, 2))
  EXP.scores <- apply(P.izf(alpha, delta, taus, int.lims, C) * scores.arr, 1:2, sum)
  res        <- list(th =  int.lims, observed = matrix(NA, N.groups, I), expected = EXP.scores)
  for (int in 1:N.groups) {
    pos.obs     <- which( (theta >= int.lims[int]) & (theta < int.lims[int + 1]))
    if (length(pos.obs) > 0) {
      res[[2]][int, ] <- colMeans(OBS.scores[pos.obs, , drop = FALSE], na.rm = TRUE)
    }
  }
  #
  if (is.null(items)) {I.plot <- 1:I} else {I.plot <- items}
  for (i in 1:length(I.plot)) {
    it <- I.plot[i]
    invisible(readline(prompt="Press [enter] to continue"))
    #
    plot(res[[1]], res[[3]][, it], type = "l", lty = 1, lwd = 2,
         xlim = c(theta.exp[1], theta.exp[1000]),
         ylim = c(0, ceiling(max(res[[3]][, it], na.rm = TRUE))),
         xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = paste0("Item Characteristic Curve (Item", it, ")"))
    points(res[[1]], res[[2]][, it], type = "p", pch = 21, bg = "grey")
    axis(1, at = seq(floor(theta.exp[1]), ceiling(theta.exp[1000]), 1))
    axis(2, at = seq(0, ceiling(max(res[[3]][, it], na.rm = TRUE)), 1), las = 1)
    mtext(expression(theta), side = 1, line = 2.5, cex = 1.2)
    mtext("Item Score", side = 2, line = 2.5, cex = 1.2)
    legend("topleft", c("Expected", "Observed"), col = "black",
           lty = c(1, NA), pch = c(NA, 21), pt.bg = "grey", lwd = c(2, 1), inset = .01, cex = .8, pt.cex = 1)
  }
  #
  scores.arr.N <- aperm(array(rep(scores.mat, N), dim = c(I, C.max + 1, N)), c(3, 1, 2))
  EXP.scores.N <- apply(P.izf(alpha, delta, taus, theta, C) * scores.arr.N, 1:2, sum)
  return(cor.OBS.EXP.items = round(diag(cor(OBS.scores, EXP.scores.N, use = "pairwise.complete.obs")), 4))
}

# Plot test information ----
plot.TestInf <- function(data, C, IP.res, Th.res) {
  if (length(C) == 1) {C <- rep(C, length(IP.res$alpha))}
  alpha     <- IP.res$alpha
  delta     <- IP.res$delta
  taus      <- IP.res$taus
  theta     <- Th.res[[2]]
  #
  d2res         <- d2logP.dtheta2.arr(data, alpha, delta, taus, theta, C)
  # Nodes:
  N.nodes <- d2res$N.nodes
  nodes   <- d2res$nodes
  #
  probs         <- P.izf(alpha, delta, taus, nodes, C)
  res           <- -apply(probs * (d2res$d2logP.dtheta2), 1, sum, na.rm = TRUE)
  res           <- cbind(nodes, res)
  colnames(res) <- c("Theta", "Inf")
  #
  plot(res[, 1], res[, 2], type = "l", lty = 1, lwd = 2,
       xlim = c(nodes[1], nodes[N.nodes]),
       ylim = c(0, ceiling(max(res[, 2], na.rm = TRUE)) + 1),
       xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = "Test Information Function")
  axis(1, at = seq(floor(nodes[1]), ceiling(nodes[N.nodes]), 1))
  axis(2, at = seq(0, ceiling(max(res[, 2], na.rm = TRUE) + 1), 10), las = 1)
  mtext(expression(theta), side = 1, line = 2.5, cex = 1.2)
  mtext("Information", side = 2, line = 2.5, cex = 1.2)
  #
  return(TestInformation = res)
}

# Plot item information ----
plot.ItemInf <- function(data, C, IP.res, Th.res, items = NULL) {
  if (length(C) == 1) {C <- rep(C, length(IP.res$alpha))}
  alpha     <- IP.res$alpha
  delta     <- IP.res$delta
  taus      <- IP.res$taus
  theta     <- Th.res[[2]]
  I         <- length(alpha)
  #
  d2res         <- d2logP.dtheta2.arr(data, alpha, delta, taus, theta, C)
  # Nodes:
  N.nodes <- d2res$N.nodes
  nodes   <- d2res$nodes
  #
  probs         <- P.izf(alpha, delta, taus, nodes, C)
  res           <- -apply(probs * (d2res$d2logP.dtheta2), 1:2, sum, na.rm = TRUE)
  res           <- cbind(nodes, res)
  colnames(res) <- c("Theta", paste0("Inf", 1:I))
  #
  if (is.null(items)) {I.plot <- 1:I} else {I.plot <- items}
  for (i in 1:length(I.plot)) {
    it <- I.plot[i]
    invisible(readline(prompt="Press [enter] to continue"))
    #
    plot(res[, 1], res[, it + 1], type = "l", lty = 1, lwd = 2,
         xlim = c(nodes[1], nodes[N.nodes]),
         ylim = c(0, ceiling(max(res[, 2:(I + 1)], na.rm = TRUE)) + 1),
         xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = paste0("Item Information Function (Item", it, ")"))
    axis(1, at = seq(floor(nodes[1]), ceiling(nodes[N.nodes]), 1))
    axis(2, at = seq(0, ceiling(max(res[, 2:(I + 1)], na.rm = TRUE) + 1), 2), las = 1)
    mtext(expression(theta), side = 1, line = 2.5, cex = 1.2)
    mtext("Information", side = 2, line = 2.5, cex = 1.2)
  }
  #
  return(ItemInformation = res[, c(1, I.plot + 1)])
}
