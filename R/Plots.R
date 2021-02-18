# Plot CRCs ----
#' @title Plot item category response curves (CRCs)
#'   
#' @description \code{plot.CRC} plots item CRCs for the GUM and the GGUM.
#'   
#' @param IP Object of class \code{GGUM}.
#' @param items Vector indicating the items for which the CRCs are to be 
#'   plotted. Default is all items.
#' @param x.lim Controls the limits of the x-axis. Default is -4 through +4.
#' @param ThetaminDelta Logical; if \code{TRUE}, plot the CRCs centered at 0, 
#'   otherwise plot the CRCs centered at \eqn{\delta}{delta} (item's
#'   difficulty). Default is \code{TRUE}.
#' @param quiet Render all plots for \code{items} at once? Default is 
#' \code{FALSE}.
#'   
#' @return The function returns a three-dimensional array with the probabilities
#'   associated to each item's CRC. These are the values shown in the plot.
#'   
#' @section Details: This function plots the item category response curves
#'   (CRCs) for the requested items.
#'   
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' # For GUM:
#' # Generate data:
#' gen1 <- GenData.GGUM(400, 5, 3, "GUM", seed = 139)
#' # Fit the GUM:
#' fit1 <- GUM(gen1$data, 3)
#' # Plot CRCs:
#' plotCRC(fit1, items = 1, quiet = TRUE)
#' \dontrun{
#' # For GGUM:
#' # Generate data:
#' set.seed(1); C <- sample(3:5, 10, replace = TRUE)
#' gen2 <- GenData.GGUM(2000, 10, C, "GGUM", seed = 156)
#' # Fit the GGUM:
#' fit2 <- GGUM(gen2$data, C)
#' # Plot CRCs:
#' plotCRC(fit2, items = 1, quiet = TRUE)
#' }
#' 
#' @export
plotCRC <- function(IP, items = NULL, x.lim = 4, ThetaminDelta = TRUE, 
                    quiet = FALSE) 
{
  # Sanity check - class:
  Sanity.class(IP)
  
  C <- IP$C
  M     <- 2 * C + 1
  C.max <- max(C)
  
  I         <- ncol(IP$data)
  alpha     <- IP$alpha
  delta     <- IP$delta
  taus      <- IP$taus
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
  
  probs.arr <- array(0, c(100, C.max + 2, length(I.plot)))
  for (i in 1:length(I.plot)) probs.arr[, 1, i] <- x
  dimnames(probs.arr)[[2]] <- c("Theta", paste0("C=", 0:C.max))
  dimnames(probs.arr)[[3]] <- paste0("Item", I.plot)
  
  for (i in 1:length(I.plot)) {
    it <- I.plot[i]
    if (!quiet) invisible(readline(prompt="Press [enter] to continue"))
    #
    probs.arr[, 2, i] <- P.GGUM(0, alpha[it], delta[it], 
                                taus[it, (C.max - C[it] + 1):(C.max - C[it] + M[it])], 
                                th.vals[it, ], C[it])
    # 
    if (C[it] <= 3) par(mar = c(6, 4.5, 1.5, .5)) else 
    par(mar = c(7, 4.5, 1.5, .5))
    plot(x, probs.arr[, 2, i],
         type = "l", lty = 1, col = viridis(C[it] + 1)[1], lwd = 2,
         xlim = c(-x.lim, x.lim), ylim = 0:1,
         xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
         main = paste0("CRC (Item", it, ")"))
    for (c in 1:C[it]) {
      probs.arr[, c+2, i] <- P.GGUM(c, alpha[it], delta[it], 
                                  taus[it, (C.max - C[it] + 1):(C.max - C[it] + M[it])], 
                                  th.vals[it, ], C[it])
      points(x, probs.arr[, c+2, i],
             type = "l", lty = 1, col = viridis(C[it] + 1)[c + 1], lwd = 2,
             xlim = c(-x.lim, x.lim), ylim = 0:1)
    }
    axis(1, at = seq(-x.lim, x.lim, 1))
    axis(2, at = seq(0, 1, .20), las = 1)
    if (ThetaminDelta == TRUE) 
    {mtext(expression(paste(theta, " - ", delta)), 
           side = 1, line = 2.5, cex = 1.2)} else {
             mtext(expression(paste("Location ", delta)), side = 1, 
                   line = 2.5, cex = 1.2)
           }
    mtext(expression(paste("P(X=c|", theta, ")")), side = 2, line = 2.5, 
          cex = 1.2)
    #
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 4.5, 0, .5), 
        new = TRUE)
    plot(0, type="n", axes=F, xlab="", ylab="")
    if (C[it] <= 3)
    {
      legend("bottom", paste0("C = ", 0:C[it]), col = viridis(C[it] + 1),
             lty = 1, lwd = 2, inset = .01, cex = .8, horiz = TRUE, 
             bg = "gray95")
    } else {
      legend("bottom", paste0("C = ", 0:C[it]), col = viridis(C[it] + 1),
             lty = 1, lwd = 2, inset = .01, cex = .8, 
             ncol = (C[it]+1) %/% 2 + (C[it]+1) %% 2, 
             bg = "gray95")
    }
  }
  invisible(probs.arr)
}

# Plot TCC ----
#' @title Plot test characteristic curve (TCC)
#'   
#' @description \code{plot.TCC} plots the TCC for the GUM and the GGUM.
#'   
#' @param IP Object of class \code{GGUM}.
#' @param Th Theta estimates from function \code{Theta.EAP()}.
#'   
#' @return The function returns a list with three elements: \item{coords}{(x, y)
#'   coordinates of the TCC.} \item{cor.OBS.EXP}{Correlation between observed
#'   and expected test scores (missing values pairwise removed).} 
#'   \item{cor.OBS.EXP.means}{Correlation between observed and expected mean
#'   test scores (missing values pairwise removed). The \eqn{\theta}{theta}
#'   interval between \eqn{-4}{-4} through \eqn{+4}{+4} is divided in 100
#'   subintervals of equal length. The observed and expected mean scores are
#'   computed for each subinterval.}
#'   
#' @section Details: This function plots the test characteristic curve (TCC).
#'   
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' \dontrun{
#' # For GUM:
#' # Generate data
#' #   (toy example: Too few items (due to computation time constraints) for 
#' #   accurate estimation of person parameters; larger number of items is 
#' #   required in practice):
#' gen1 <- GenData.GGUM(400, 5, 3, "GUM", seed = 139)
#' # Fit the GUM:
#' fit1 <- GUM(gen1$data, 3)
#' th1  <- Theta.EAP(fit1)
#' # Plot TCC:
#' plotTCC(fit1, th1)
#' }
#' \dontrun{
#' # For GGUM:
#' # Generate data:
#' set.seed(1); C <- sample(3:5, 10, replace = TRUE)
#' gen2 <- GenData.GGUM(2000, 10, C, "GGUM", seed = 156)
#' # Fit the GGUM:
#' fit2 <- GGUM(gen2$data, C)
#' th2  <- Theta.EAP(fit2)
#' # Plot TCC:
#' plotTCC(fit2, th2)
#' }
#' 
#' @export
plotTCC <- function(IP, Th) 
{
  # Sanity check - class:
  Sanity.class(IP)
  
  C <- IP$C
  M     <- 2 * C + 1
  C.max <- max(C)
  
  I         <- ncol(IP$data)
  N         <- nrow(IP$data)
  alpha     <- IP$alpha
  delta     <- IP$delta
  taus      <- IP$taus
  theta     <- if (class(Th) == "matrix") Th[, "Theta"] else Th
  
  # Test characteristic curve:
  OBS.scores      <- IP$data
  theta.exp       <- seq(min(-4, min(theta, na.rm = TRUE)), 
                         max(4, max(theta, na.rm = TRUE)), 
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
    pos.obs <- which( (theta >= int.lims[int]) & (theta < int.lims[int + 1]))
    if (length(pos.obs) > 0) {
      res[int, 2] <- mean(rowSums(OBS.scores[pos.obs, , drop = FALSE], 
                                  na.rm = TRUE))
    }
  }
  #
  par(mar = c(6, 4.5, 1.5, .5))
  plot(res[, 1], res[, 3], type = "l", lty = 1, lwd = 2,
       xlim = c(theta.exp[1], theta.exp[1000]),
       ylim = c(0, ceiling(max(res[, 3], na.rm = TRUE))),
       xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
       main = "Test Characteristic Curve")
  points(res[, 1], res[, 2], type = "p", pch = 21, bg = "grey")
  axis(1, at = seq(floor(theta.exp[1]), ceiling(theta.exp[1000]), 1))
  axis(2, at = seq(0, ceiling(max(res[, 3], na.rm = TRUE)), length.out = 5), 
       las = 1)
  mtext(expression(theta), side = 1, line = 2.5, cex = 1.2)
  mtext("Total Score", side = 2, line = 3, cex = 1.2)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 4.5, 0, .5), 
      new = TRUE)
  plot(0, type="n", axes=F, xlab="", ylab="")
  legend("bottom", c("Expected", "Observed"), col = "black",
         lty = c(1, NA), pch = c(NA, 21), pt.bg = "grey", lwd = c(2, 1), 
         inset = .01, cex = .8, pt.cex = 1, horiz = TRUE, bg = "gray95")
  
  #
  cor.OBS.EXP.means <- cor(res[, 2], res[, 3], use = "pairwise.complete.obs")
  scores.arr.N      <- aperm(array(rep(scores.mat, N), 
                                   dim = c(I, C.max + 1, N)), c(3, 1, 2))
  EXP.scores.N      <- apply(P.izf(alpha, delta, taus, theta, C) * scores.arr.N,
                             1:2, sum)
  cor.OBS.EXP.raw   <- cor(c(OBS.scores), c(EXP.scores.N), 
                           use = "pairwise.complete.obs")
  #
  colnames(res) <- c("Theta", "TCC.obs", "TCC.exp")
  invisible(list(coords = res, cor.OBS.EXP = round(cor.OBS.EXP.raw, 4), 
              cor.OBS.EXP.means = round(cor.OBS.EXP.means, 4)))
}

# Plot ICCs ----
#' @title Plot item characteristic curves (ICCs)
#'   
#' @description \code{plot.ICC} plots the ICCs for the GUM and the GGUM.
#'   
#' @param IP Object of class \code{GGUM}.
#' @param Th Theta estimates from function \code{Theta.EAP()}.
#' @param items Vector indicating the items for which the ICCs are to be 
#'   plotted. Default is all items.
#' @param quiet Render all plots for \code{items} at once? Default is 
#' \code{FALSE}.
#'   
#' @return The function returns the correlation between observed and expected 
#'   item scores (missing values pairwise removed).
#'   
#' @section Details: This function plots the item characteristic curves (ICCs).
#'   
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' \dontrun{
#' # For GUM:
#' # Generate data
#' #   (toy example: Too few items (due to computation time constraints) for 
#' #   accurate estimation of person parameters; larger number of items is 
#' #   required in practice):
#' gen1 <- GenData.GGUM(400, 5, 3, "GUM", seed = 139)
#' # Fit the GUM:
#' fit1 <- GUM(gen1$data, 3)
#' th1  <- Theta.EAP(fit1)
#' # Plot ICCs:
#' plotICC(fit1, th1, items = 1, quiet = TRUE)
#' }
#' \dontrun{
#' # For GGUM:
#' # Generate data:
#' set.seed(1); C <- sample(3:5, 10, replace = TRUE)
#' gen2 <- GenData.GGUM(2000, 10, C, "GGUM", seed = 156)
#' # Fit the GGUM:
#' fit2 <- GGUM(gen2$data, C)
#' th2  <- Theta.EAP(fit2)
#' # Plot ICCs:
#' plotICC(fit2, th2, items = 1, quiet = TRUE)
#' }
#' 
#' @export
plotICC <- function(IP, Th, items = NULL, quiet = FALSE) 
{
  # Sanity check - class:
  Sanity.class(IP)
  
  C <- IP$C
  M     <- 2 * C + 1
  C.max <- max(C)
  
  I         <- ncol(IP$data)
  N         <- nrow(IP$data)
  alpha     <- IP$alpha
  delta     <- IP$delta
  taus      <- IP$taus
  theta     <- if (class(Th) == "matrix") Th[, "Theta"] else Th
  
  # Item characteristic curves:
  OBS.scores      <- IP$data
  theta.exp       <- seq(min(-4, min(theta, na.rm = TRUE)), 
                         max(4, max(theta, na.rm = TRUE)), 
                         length.out = 1000)
  N.groups        <- 50
  int.lims        <- quantile(theta.exp, probs = seq(0, 1, 1 / (N.groups - 1)))
  names(int.lims) <- NULL
  #
  scores.mat <- matrix(0, nrow = I, ncol = C.max + 1)
  for (i in 1:I) {
    scores.mat[i, 1:(C[i] + 1)] <- 0:C[i]
  }
  scores.arr <- aperm(array(rep(scores.mat, N.groups), 
                            dim = c(I, C.max + 1, N.groups)), c(3, 1, 2))
  EXP.scores <- apply(P.izf(alpha, delta, taus, int.lims, C) * scores.arr, 1:2, 
                      sum)
  res        <- list(th =  int.lims, observed = matrix(NA, N.groups, I), 
                     expected = EXP.scores)
  for (int in 1:N.groups) {
    pos.obs <- which( (theta >= int.lims[int]) & (theta < int.lims[int + 1]))
    if (length(pos.obs) > 0) {
      res[[2]][int, ] <- colMeans(OBS.scores[pos.obs, , drop = FALSE], 
                                  na.rm = TRUE)
    }
  }
  #
  if (is.null(items)) {I.plot <- 1:I} else {I.plot <- items}
  for (i in 1:length(I.plot)) {
    it <- I.plot[i]
    if (!quiet) invisible(readline(prompt="Press [enter] to continue"))
    #
    par(mar = c(6, 4.5, 1.5, .5))
    plot(res[[1]], res[[3]][, it], type = "l", lty = 1, lwd = 2,
         xlim = c(theta.exp[1], theta.exp[1000]),
         ylim = c(0, ceiling(max(res[[3]][, it], na.rm = TRUE))),
         xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
         main = paste0("ICC (Item", it, ")"))
    points(res[[1]], res[[2]][, it], type = "p", pch = 21, bg = "grey")
    axis(1, at = seq(floor(theta.exp[1]), ceiling(theta.exp[1000]), 1))
    axis(2, at = seq(0, ceiling(max(res[[3]][, it], na.rm = TRUE)), 1), las = 1)
    mtext(expression(theta), side = 1, line = 2.5, cex = 1.2)
    mtext("Item Score", side = 2, line = 2.5, cex = 1.2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 4.5, 0, .5), 
        new = TRUE)
    plot(0, type="n", axes=F, xlab="", ylab="")
    legend("bottom", c("Expected", "Observed"), col = "black",
           lty = c(1, NA), pch = c(NA, 21), pt.bg = "grey", lwd = c(2, 1), 
           inset = .01, cex = .8, pt.cex = 1, horiz = TRUE, bg = "gray95")
  }
  #
  scores.arr.N <- aperm(array(rep(scores.mat, N), dim = c(I, C.max + 1, N)), c(3, 1, 2))
  EXP.scores.N <- apply(P.izf(alpha, delta, taus, theta, C) * scores.arr.N, 1:2,
                        sum)
  cor.OBS.EXP  <- round(diag(cor(OBS.scores, EXP.scores.N, use = "pairwise.complete.obs")), 4)
  names(cor.OBS.EXP) <- paste0("Item", 1:I)
  cat("\n")
  invisible(cbind(Theta = res[[1]], 
                  ICC.obs = res[[2]][, I.plot], 
                  ICC.exp = res[[3]][, I.plot]))
}

# Plot TIF ----
#' @title Plot test information function (TIF)
#'   
#' @description \code{plot.TIF} plots the TIF for the GUM and the GGUM.
#'   
#' @param IP Object of class \code{GGUM}.
#' @param Th Theta estimates from function \code{Theta.EAP()}.
#'   
#' @return The function returns the (x, y) coordinates of the TIF.
#'   
#' @section Details: This function plots the test information function (TIF).
#'   
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' \dontrun{
#' # For GUM:
#' # Generate data
#' #   (toy example: Too few items (due to computation time constraints) for 
#' #   accurate estimation of person parameters; larger number of items is 
#' #   required in practice):
#' gen1 <- GenData.GGUM(400, 5, 3, "GUM", seed = 139)
#' # Fit the GUM:
#' fit1 <- GUM(gen1$data, 3)
#' th1  <- Theta.EAP(fit1)
#' # Plot TIF:
#' plotTIF(fit1, th1)
#' }
#' \dontrun{
#' # For GGUM:
#' # Generate data:
#' set.seed(1); C <- sample(3:5, 10, replace = TRUE)
#' gen2 <- GenData.GGUM(2000, 10, C, "GGUM", seed = 156)
#' # Fit the GGUM:
#' fit2 <- GGUM(gen2$data, C)
#' th2  <- Theta.EAP(fit2)
#' # Plot TIF:
#' plotTIF(fit2, th2)
#' }
#' 
#' @export
plotTIF <- function(IP, Th) 
{
  # Sanity check - class:
  Sanity.class(IP)
  
  data  <- IP$data
  C     <- IP$C
  alpha <- IP$alpha
  delta <- IP$delta
  taus  <- IP$taus
  theta <- if (class(Th) == "matrix") Th[, "Theta"] else Th
  
  d2res         <- d2logP.dtheta2.arr(data, alpha, delta, taus, theta, C)
  # Nodes:
  N.nodes <- d2res$N.nodes
  nodes   <- d2res$nodes
  #
  probs         <- P.izf(alpha, delta, taus, nodes, C)
  res           <- -apply(probs * (d2res$d2logP.dtheta2), 1, sum, na.rm = TRUE)
  res           <- cbind(nodes, res)
  colnames(res) <- c("Theta", "TIF")
  #
  par(mar = c(3.5, 5, 1.5, .5))
  plot(res[, 1], res[, 2], type = "l", lty = 1, lwd = 2,
       xlim = c(nodes[1], nodes[N.nodes]),
       ylim = c(0, ceiling(max(res[, 2], na.rm = TRUE)) + 1),
       xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
       main = "Test Information Function")
  axis(1, at = seq(floor(nodes[1]), ceiling(nodes[N.nodes]), 1))
  axis(2, at = seq(0, ceiling(max(res[, 2], na.rm = TRUE) + 1), length.out = 5),
       las = 1)
  mtext(expression(theta), side = 1, line = 2.5, cex = 1.2)
  mtext("Information", side = 2, line = 2.5, cex = 1.2)
  #
  invisible(res)
}

# Plot IIFs ----
#' @title Plot item information functions (IIFs)
#'   
#' @description \code{plot.IIF} plots the IIFs for the GUM and the GGUM.
#'   
#' @param IP Object of class \code{GGUM}.
#' @param Th Theta estimates from function \code{Theta.EAP()}.
#' @param items Vector indicating the items for which the ICCs are to be 
#'   plotted. Default is all items.
#' @param quiet Render all plots for \code{items} at once? Default is 
#' \code{FALSE}.
#'   
#' @return The function returns the (x, y) coordinates of the IIFs.
#'   
#' @section Details: This function plots the item information functions (IIFs).
#'   
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' \dontrun{
#' # For GUM:
#' # Generate data
#' #   (toy example: Too few items (due to computation time constraints) for 
#' #   accurate estimation of person parameters; larger number of items is 
#' #   required in practice):
#' gen1 <- GenData.GGUM(400, 5, 3, "GUM", seed = 139)
#' # Fit the GUM:
#' fit1 <- GUM(gen1$data, 3)
#' th1  <- Theta.EAP(fit1)
#' # Plot IIFs:
#' plotIIF(fit1, th1, items = 1, quiet = TRUE)
#' }
#' \dontrun{
#' # For GGUM:
#' # Generate data:
#' set.seed(1); C <- sample(3:5, 10, replace = TRUE)
#' gen2 <- GenData.GGUM(2000, 10, C, "GGUM", seed = 156)
#' # Fit the GGUM:
#' fit2 <- GGUM(gen2$data, C)
#' th2  <- Theta.EAP(fit2)
#' # Plot IIFs:
#' plotIIF(fit2, th2, items = 1, quiet = TRUE)
#' }
#' 
#' @export
plotIIF <- function(IP, Th, items = NULL, quiet = FALSE) 
{
  # Sanity check - class:
  Sanity.class(IP)
  
  data  <- IP$data
  I     <- ncol(data)
  C     <- IP$C
  alpha <- IP$alpha
  delta <- IP$delta
  taus  <- IP$taus
  theta <- if (class(Th) == "matrix") Th[, "Theta"] else Th
  
  #
  d2res <- d2logP.dtheta2.arr(data, alpha, delta, taus, theta, C)
  # Nodes:
  N.nodes <- d2res$N.nodes
  nodes   <- d2res$nodes
  #
  probs         <- P.izf(alpha, delta, taus, nodes, C)
  res           <- -apply(probs * (d2res$d2logP.dtheta2), 1:2, sum, 
                          na.rm = TRUE)
  res           <- cbind(nodes, res)
  colnames(res) <- c("Theta", paste0("IIF", 1:I))
  #
  if (is.null(items)) {I.plot <- 1:I} else {I.plot <- items}
  for (i in 1:length(I.plot)) {
    it <- I.plot[i]
    if (!quiet) invisible(readline(prompt="Press [enter] to continue"))
    #
    par(mar = c(3.5, 5, 1.5, .5))
    plot(res[, 1], res[, it + 1], type = "l", lty = 1, lwd = 2,
         xlim = c(nodes[1], nodes[N.nodes]),
         ylim = c(0, ceiling(max(res[, 2:(I + 1)], na.rm = TRUE))),
         xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
         main = paste0("IIF (Item", it, ")"))
    axis(1, at = seq(floor(nodes[1]), ceiling(nodes[N.nodes]), 1))
    axis(2, at = seq(0, ceiling(max(res[, 2:(I + 1)], na.rm = TRUE)), 
                     length.out = 5), las = 1)
    mtext(expression(theta), side = 1, line = 2.5, cex = 1.2)
    mtext("Information", side = 2, line = 2.5, cex = 1.2)
  }
  #
  invisible(res[, c(1, I.plot + 1)])
}
