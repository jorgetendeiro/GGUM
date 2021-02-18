# Computing chisq/df ratios for single, pairs, and triples of items 
# (Drasgow et al., 1995; LaHuis, Clark, & O'Brien, 2011).
# Run Stark's MODFIT Excel program and look into file chisqr.dbg, it is very 
# informative.
# For the adjusted chisq to N = 3000, see LaHuis et al. (2011, p. 14).
# 

#' @title MODFIT for the GGUM
#'   
#' @description \code{MODFIT} computes the adjusted \eqn{\chi^2}{chi-square} 
#'   degrees of freedom ratios (\eqn{\chi^2/df}{chisq/df}) introduced by Drasgow
#'   et al. (1995) for the GGUM.
#'   
#' @param IP Object of class \code{GGUM}.
#' @param precision Number of decimal places of the results (default = 4).
#'   
#' @return A list (an object of class \code{MODFIT}) with four elements: The 
#'   results for singlets, doublets, triples, and a summary result.
#'   
#' @section Details: This function computes the adjusted
#'   \eqn{\chi^2}{chi-square} degrees of freedom ratios
#'   (\eqn{\chi^2/df}{chisq/df}) introduced by Drasgow et al. (1995). These
#'   \eqn{\chi^2}{chi-square} statistics are based on expected frequencies that
#'   depend on the estimated item parameters and the distribution of
#'   \eqn{\theta}{theta}. The \emph{unadjusted} statistic for item \eqn{i} is
#'   given by
#'   
#'   \deqn{\chi^2_i = \sum_{z=0}^C \frac{(O_{iz} - E_{iz})^2}{E_{iz}}, } 
#'   {chisq_i = sum( (O_iz - E_iz)^2 / E_iz; z = 0, ..., C ), }
#'   
#'   with
#'   
#'   \deqn{E_{iz} = N\int P_{iz}(\theta)\varphi(\theta)d\theta.}{E_iz = N 
#'   int(P_iz(th)phi(th)dth).}
#'   
#'   \eqn{O_{iz}}{O_iz} is the observed frequency of choosing answer \eqn{z} for
#'   item \eqn{i} and \eqn{\varphi(\theta)}{phi(th)} is the standard normal 
#'   density. The equation above applies to single items ('singlets'). The 
#'   formula is easily extendible to pairs and triples of items. For a large
#'   number of items, the function selects suitable subsets of doublets and
#'   triples to perform the computations since its total number increases
#'   quickly with test length (Drasgow et al., 1995).
#'   
#'   The formula is adjusted to a sample size of 3,000, as follows (see also 
#'   LaHuis et al., 2011):
#'   
#'   \deqn{\chi^2_i/df = \frac{3,000(\chi^2_i - df)}{N}+df,}{chisq/df = 
#'   3,000(chsqr-df)/N + df,}
#'   
#'   where \eqn{df} is a number of degrees of freedom that depends on the number
#'   of singlets, doublets, and triplets.
#'   
#'   As an heuristic, values of \eqn{\chi^2/df}{chisq/df} larger than 3 are 
#'   indicative of model misfit.
#'   
#'   This function produces the same numerical results as the MODFIT program 
#'   (Stark, 2001) for the GGUM.
#'   
#' @references \insertRef{Drasgowetal1995}{GGUM}
#' 
#' \insertRef{LaHuisetal2011}{GGUM}
#' 
#' \insertRef{MODFITsoftware}{GGUM}
#' 
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' # For GUM:
#' # Generate data:
#' gen1 <- GenData.GGUM(400, 5, 3, "GUM", seed = 139)
#' # Fit the GUM:
#' fit1 <- GUM(gen1$data, 3)
#' # Compute the adjusted chi square degrees of freedom ratios:
#' modfit.res1 <- MODFIT(fit1)
#' modfit.res1$Singlets
#' modfit.res1$Doublets
#' modfit.res1$Triplets
#' modfit.res1$Summary
#' \dontrun{
#' # For GGUM:
#' # Generate data:
#' set.seed(1); C <- sample(3:5, 10, replace = TRUE)
#' gen2 <- GenData.GGUM(2000, 10, C, "GGUM", seed = 156)
#' # Fit the GGUM:
#' fit2 <- GGUM(gen2$data, C)
#' # Compute the adjusted chi square degrees of freedom ratios:
#' modfit.res2 <- MODFIT(fit1)
#' modfit.res2$Singlets
#' modfit.res2$Doublets
#' modfit.res2$Triplets
#' modfit.res2$Summary
#' }
#' 
#' @export
MODFIT <- function(IP, precision = 4)
{
  data  <- IP$data
  C     <- IP$C
  model <- IP$model
  
  N     <- nrow(data)
  I     <- ncol(data)
  N.NAs <- N - colSums(is.na(data))
  
  if (I <= 10)
  {
    doublets <- t(combn(I, 2))
    triplets <- t(combn(I, 3))
  } else {
    # Find packets:
    obs.props   <- colMeans(data > 0, na.rm = TRUE)
    group.low   <- sort(order(obs.props)[1:ceiling(I / 3)])
    group.med   <- sort(order(obs.props)[(ceiling(I / 3) + 1) : ceiling(2*I/3)])
    group.high  <- sort(order(obs.props)[(ceiling(2*I / 3) + 1) : I])
    groups      <- cbind(group.low = group.low[1:floor(I / 3)], 
                         group.med[1:floor(I / 3)], 
                         group.high[1:floor(I / 3)])
    packets     <- lapply(seq_len(nrow(groups)), function(row) sort(groups[row, ]))
    if ((I %% 3) == 1) {packets[[1]][4] <- group.low[ceiling(I / 3)]}
    if ((I %% 3) == 2) {
      packets[[1]][4] <- group.low[ceiling(I / 3)]
      packets[[2]][4] <- group.med[ceiling(I / 3)]
    }
    
    # 
    # singlets    <- 1:I
    doublets    <- matrix(unlist(lapply(packets, function(x) combn(x,2))), 
                          ncol = 2, byrow = TRUE)
    triplets    <- matrix(unlist(lapply(packets, function(x) combn(x,3))), 
                          ncol = 3, byrow = TRUE)
  }
  
  # NAs for doublets and triplets:
  N.NAs.doublets <- apply(doublets, 1, 
                          function(vec) N - sum(rowSums(is.na(data[, vec])) > 0))
  N.NAs.triplets <- apply(triplets, 1, 
                          function(vec) N - sum(rowSums(is.na(data[, vec])) > 0))
  
  # Nodes and weights:
  nodes.chi   <- seq(-3, 3, length.out = 61)
  N.nodes.chi <- length(nodes.chi)
  weights     <- dnorm(nodes.chi) / sum(dnorm(nodes.chi))
  
  # Singlets:
  probs.array.aber.drasgow <- array(NA, dim = c(N.nodes.chi, I, max(C) + 1))
  if (model == "GGUM")
  {
    if (length(C) == 1)
    {
      for (z in 0:C) 
      {
        probs.array.aber.drasgow[, , z + 1] <- 
          P.GGUM(z, IP$alpha, IP$delta, IP$taus, nodes.chi, C)
      }
    } else
    {
      for (i in 1:I)
      {
        for (z in 0:C[i]) 
        {
          probs.array.aber.drasgow[, i, z + 1] <- 
          P.GGUM(z, IP$alpha[i], IP$delta[i], 
                 IP$taus[i, (max(C)-C[[i]]+1):(2*max(C)+1-(max(C)-C[i]))], 
                 nodes.chi, C[i])
        }
      }
    }
  }
  if (model == "GRM")
  {
    probs.array.aber.drasgow <- P.GRM(C, IP, nodes.chi)
  }
  weights.arr              <- array(rep(weights, I * (max(C) + 1)), 
                                    c(N.nodes.chi, I, max(C) + 1))
  N.NAs.mat                <- matrix(rep(N.NAs, max(C) + 1), nrow = I, 
                                     byrow = FALSE)
  if (length(C) > 1) for (i in 1:I) N.NAs.mat[i, (C[i] + 1):(max(C) + 1)] <- NA
  expected.mat.drasgow     <- N.NAs.mat * apply((probs.array.aber.drasgow * weights.arr), 2:3, sum)
  if (length(C) == 1) 
  {
    observed.mat.drasgow     <- t(apply(data, 2, 
                                        function(vec) table(factor(vec, levels = 0:C))))
  } else 
  {
    observed.mat.drasgow <- matrix(NA, nrow = I, ncol = max(C) + 1)
    for (i in 1:I) 
    {
      observed.mat.drasgow[i, 1:(C[i] + 1)] <- table(factor(data[, i], levels = 0:C[i]))
    }
  }
  # Merge cells with expected frequencies < 5:
  expected.order                 <- t(apply(expected.mat.drasgow, 1, order))
  expected.mat.drasgow           <- t(sapply(1:I, function(it) expected.mat.drasgow[it, expected.order[it, ]]))
  observed.mat.drasgow           <- t(sapply(1:I, function(it) observed.mat.drasgow[it, expected.order[it, ]]))
  expected.mat.drasgow.less5     <- rowSums(expected.mat.drasgow < 5, na.rm = TRUE)
  N.expected.mat.drasgow.less5   <- sum(expected.mat.drasgow.less5 > 0)
  pos.expected.mat.drasgow.less5 <- which(expected.mat.drasgow.less5 > 0)
  df                             <- if (length(C) == 1) rep(C, I) else C
  if (N.expected.mat.drasgow.less5 > 0)
  {
    sapply(1:N.expected.mat.drasgow.less5, function(it) 
    {
      item    <- pos.expected.mat.drasgow.less5[it]
      pos.sum <- expected.mat.drasgow.less5[item]
      if (sum(expected.mat.drasgow[item, 1:pos.sum]) < 5) {pos.sum <- pos.sum + 1}
      expected.mat.drasgow[item, pos.sum]         <- sum(expected.mat.drasgow[item, 1:pos.sum])
      expected.mat.drasgow[item, 1:(pos.sum - 1)] <- 1
      observed.mat.drasgow[item, pos.sum]         <- sum(observed.mat.drasgow[item, 1:pos.sum])
      observed.mat.drasgow[item, 1:(pos.sum - 1)] <- 1
      df[item] <- if (length(C) == 1) (C + 1) - pos.sum else (C[item] + 1) - pos.sum
    })
  }
  # Compute (adjusted) chi squares (/df):
  chisq        <- rowSums(((observed.mat.drasgow - expected.mat.drasgow)^2) / expected.mat.drasgow, na.rm = TRUE)
  chisq.df     <- chisq / df
  chisq.adj    <- sapply(1:I, function(it) max(0, 3000 * (chisq[it] - df[it]) / N.NAs[it] + df[it]))
  chisq.adj.df <- chisq.adj / df
  singlets.res <- cbind(Item = 1:I, N.NAs, df, chisq, chisq.df, chisq.adj, chisq.adj.df)
  
  # Doublets:
  doublets.NAs <- cbind(doublets, N.NAs.doublets)
  doublets.res <- t(apply(doublets.NAs, 1, function(vec)
  {
    item1   <- vec[1]
    item2   <- vec[2]
    N.NAs.d <- vec[3]
    probs.array.aber.drasgow.item1     <- probs.array.aber.drasgow[, item1, ]
    probs.array.aber.drasgow.item2     <- probs.array.aber.drasgow[, item2, ]
    probs.array.aber.drasgow.item1.arr <- array(rep(probs.array.aber.drasgow.item1, max(C) + 1), c(N.nodes.chi, max(C) + 1, max(C) + 1))
    probs.array.aber.drasgow.item2.arr <- array(rep(probs.array.aber.drasgow.item2, max(C) + 1), c(N.nodes.chi, max(C) + 1, max(C) + 1))
    probs.array.aber.drasgow.item2.arr <- aperm(probs.array.aber.drasgow.item2.arr, c(1, 3, 2))
    weights.arr2                       <- array(rep(weights.arr, (max(C) + 1) * (max(C) + 1)), c(N.nodes.chi, max(C) + 1, max(C) + 1)) # For doublets
    expected.mat.drasgow.it1.it2       <- N.NAs.d * apply(probs.array.aber.drasgow.item1.arr * probs.array.aber.drasgow.item2.arr * weights.arr2, 2:3, sum)
    if (length(C) > 1)
    {
      if (C[item1] < max(C)) expected.mat.drasgow.it1.it2[(C[item1] + 2):(max(C) + 1), ] <- NA
      if (C[item2] < max(C)) expected.mat.drasgow.it1.it2[, (C[item2] + 2):(max(C) + 1)] <- NA
    }
    if (length(C) == 1) 
    {
      observed.mat.drasgow.it1.it2 <- table(factor(data[, item1], levels = 0:C), factor(data[, item2], levels = 0:C)) # table(data[, item1], data[, item2])
    } else 
    {
      observed.mat.drasgow.it1.it2 <- matrix(NA, nrow = max(C) + 1, ncol = max(C) + 1)
      observed.mat.drasgow.it1.it2[1:(C[item1] + 1), 1:(C[item2] + 1)] <- table(factor(data[, item1], levels = 0:C[item1]), factor(data[, item2], levels = 0:C[item2]))
    }
    # Merge cells with expected frequencies < 5:
    expected.order.it1.it2             <- order(c(expected.mat.drasgow.it1.it2))
    expected.mat.drasgow.it1.it2       <- c(expected.mat.drasgow.it1.it2)[expected.order.it1.it2]
    observed.mat.drasgow.it1.it2       <- c(observed.mat.drasgow.it1.it2)[expected.order.it1.it2]
    expected.mat.drasgow.it1.it2.less5 <- expected.mat.drasgow.it1.it2 < 5
    if (length(C) == 1) {df.it1.it2 <- (C + 1)^2 - 1} else {df.it1.it2 <- (C[item1] + 1) * (C[item2] + 1) - 1}
    
    if (sum(expected.mat.drasgow.it1.it2.less5, na.rm = TRUE) > 0)
    {
      pos.sum <- max(which(expected.mat.drasgow.it1.it2.less5 == 1))
      if (sum(expected.mat.drasgow.it1.it2[1:pos.sum]) < 5) {pos.sum <- pos.sum + 1}
      expected.mat.drasgow.it1.it2[pos.sum] <- sum(expected.mat.drasgow.it1.it2[1:pos.sum])
      expected.mat.drasgow.it1.it2          <- expected.mat.drasgow.it1.it2[pos.sum : ((max(C) + 1)^2)]
      observed.mat.drasgow.it1.it2[pos.sum] <- sum(observed.mat.drasgow.it1.it2[1:pos.sum])
      observed.mat.drasgow.it1.it2          <- observed.mat.drasgow.it1.it2[pos.sum : ((max(C) + 1)^2)]
      if (length(C) == 1) {df.it1.it2 <- (C + 1)^2 - pos.sum} else {df.it1.it2 <- (C[item1] + 1) * (C[item2] + 1) - pos.sum}
    }
    # Compute (adjusted) chi squares (/df):
    chisq        <- sum(((observed.mat.drasgow.it1.it2 - expected.mat.drasgow.it1.it2)^2) / expected.mat.drasgow.it1.it2, na.rm = TRUE)
    chisq.df     <- chisq / df.it1.it2
    chisq.adj    <- max(0, 3000 * (chisq - df.it1.it2) / N.NAs.d + df.it1.it2)
    chisq.adj.df <- chisq.adj / df.it1.it2
    c(Item1 = item1, Item2 = item2, N = N.NAs.d, df = df.it1.it2, chisq = chisq, chisq.df = chisq.df, chisq.adj = chisq.adj, chisq.adj.df = chisq.adj.df)
  }))
  doublets.res <- cbind(Doublet = 1:nrow(doublets), doublets.res)
  
  # Triplets:
  triplets.NAs <- cbind(triplets, N.NAs.triplets)
  triplets.res <-  t(apply(triplets.NAs, 1, function(vec)
  {
    item1   <- vec[1]
    item2   <- vec[2]
    item3   <- vec[3]
    N.NAs.t <- vec[4]
    probs.array.aber.drasgow.item1     <- probs.array.aber.drasgow[, item1, ]
    probs.array.aber.drasgow.item2     <- probs.array.aber.drasgow[, item2, ]
    probs.array.aber.drasgow.item3     <- probs.array.aber.drasgow[, item3, ]
    probs.array.aber.drasgow.item1.arr <- array(rep(probs.array.aber.drasgow.item1, (max(C) + 1) * (max(C) + 1)), c(N.nodes.chi, max(C) + 1, max(C) + 1, max(C) + 1))
    probs.array.aber.drasgow.item2.arr <- array(rep(probs.array.aber.drasgow.item2, (max(C) + 1) * (max(C) + 1)), c(N.nodes.chi, max(C) + 1, max(C) + 1, max(C) + 1))
    probs.array.aber.drasgow.item2.arr <- aperm(probs.array.aber.drasgow.item2.arr, c(1, 3, 2, 4))
    probs.array.aber.drasgow.item3.arr <- array(rep(probs.array.aber.drasgow.item3, (max(C) + 1) * (max(C) + 1)), c(N.nodes.chi, max(C) + 1, max(C) + 1, max(C) + 1))
    probs.array.aber.drasgow.item3.arr <- aperm(probs.array.aber.drasgow.item3.arr, c(1, 4, 3, 2))
    weights.arr3                       <- array(rep(weights.arr, (max(C) + 1)^3), c(N.nodes.chi, max(C) + 1, max(C) + 1, max(C) + 1))
    expected.mat.drasgow.it1.it2.it3   <- N.NAs.t * apply(probs.array.aber.drasgow.item1.arr * probs.array.aber.drasgow.item2.arr * 
                                                      probs.array.aber.drasgow.item3.arr * weights.arr3, 2:4, sum)
    if (length(C) > 1)
    {
      if (C[item1] < max(C)) expected.mat.drasgow.it1.it2.it3[(C[item1] + 2):(max(C) + 1), , ] <- NA
      if (C[item2] < max(C)) expected.mat.drasgow.it1.it2.it3[, (C[item2] + 2):(max(C) + 1), ] <- NA
      if (C[item3] < max(C)) expected.mat.drasgow.it1.it2.it3[, , (C[item3] + 2):(max(C) + 1)] <- NA
    }
    if (length(C) == 1) 
    {
      observed.mat.drasgow.it1.it2.it3   <- table(factor(data[, item1], levels = 0:C), factor(data[, item2], levels = 0:C), factor(data[, item3], levels = 0:C)) # table(data[, item1], data[, item2], data[, item3])
    } else 
    {
      observed.mat.drasgow.it1.it2.it3 <- array(NA, dim = c(max(C) + 1, max(C) + 1, max(C) + 1))
      observed.mat.drasgow.it1.it2.it3[1:(C[item1] + 1), 1:(C[item2] + 1), 1:(C[item3] + 1)] <- 
        table(factor(data[, item1], levels = 0:C[item1]), factor(data[, item2], levels = 0:C[item2]), factor(data[, item3], levels = 0:C[item3]))
    }
    # Merge cells with expected frequencies < 5:
    expected.order.it1.it2.it3             <- order(c(expected.mat.drasgow.it1.it2.it3))
    expected.mat.drasgow.it1.it2.it3       <- c(expected.mat.drasgow.it1.it2.it3)[expected.order.it1.it2.it3]
    observed.mat.drasgow.it1.it2.it3       <- c(observed.mat.drasgow.it1.it2.it3)[expected.order.it1.it2.it3]
    expected.mat.drasgow.it1.it2.it3.less5 <- expected.mat.drasgow.it1.it2.it3 < 5
    if (length(C) == 1) {df.it1.it2.it3 <- (C + 1)^3 - 1} else {df.it1.it2.it3 <- (C[item1] + 1) * (C[item2] + 1) * (C[item3] + 1) - 1}
    if (sum(expected.mat.drasgow.it1.it2.it3.less5, na.rm = TRUE) > 0)
    {
      pos.sum <- max(which(expected.mat.drasgow.it1.it2.it3.less5 == 1))
      if (sum(expected.mat.drasgow.it1.it2.it3[1:pos.sum]) < 5) {pos.sum <- pos.sum + 1}
      expected.mat.drasgow.it1.it2.it3[pos.sum] <- sum(expected.mat.drasgow.it1.it2.it3[1:pos.sum])
      expected.mat.drasgow.it1.it2.it3          <- expected.mat.drasgow.it1.it2.it3[pos.sum : ((max(C) + 1)^3)]
      observed.mat.drasgow.it1.it2.it3[pos.sum] <- sum(observed.mat.drasgow.it1.it2.it3[1:pos.sum])
      observed.mat.drasgow.it1.it2.it3          <- observed.mat.drasgow.it1.it2.it3[pos.sum : ((max(C) + 1)^3)]
      if (length(C) == 1) {df.it1.it2.it3 <- (C + 1)^3 - pos.sum} else {df.it1.it2.it3 <- (C[item1] + 1) * (C[item2] + 1) * (C[item3] + 1) - pos.sum}
    }
    # Compute (adjusted) chi squares (/df):
    chisq        <- sum(((observed.mat.drasgow.it1.it2.it3 - expected.mat.drasgow.it1.it2.it3)^2) / expected.mat.drasgow.it1.it2.it3, na.rm = TRUE)
    chisq.df     <- chisq / df.it1.it2.it3
    chisq.adj    <- max(0, 3000 * (chisq - df.it1.it2.it3) / N.NAs.t + df.it1.it2.it3)
    chisq.adj.df <- chisq.adj / df.it1.it2.it3
    c(Item1 = item1, Item2 = item2, Item3 = item3, N = N.NAs.t, df = df.it1.it2.it3, 
      chisq = chisq, chisq.df = chisq.df, chisq.adj = chisq.adj, chisq.adj.df = chisq.adj.df)
  }))
  triplets.res <- cbind(Triplet = 1:nrow(triplets), triplets.res)
  
  # Summarize results:
  f.int          <- function(x) {if (x < 1) 1 else (if (x < 2) 2 else (if (x < 3) 3 else (if (x < 4) 4 else (if (x < 5) 5 else (if (x < 7) 6 else 7)))))}
  singlets.table <- c(table(factor(sapply(singlets.res[, 7], f.int), levels=1:7)), round(mean(singlets.res[, 7]), 4), round(sd(singlets.res[, 7]), 4))
  doublets.table <- c(table(factor(sapply(doublets.res[, 9], f.int), levels=1:7)), round(mean(doublets.res[, 9]), 4), round(sd(doublets.res[, 9]), 4))
  triplets.table <- c(table(factor(sapply(triplets.res[, 10], f.int), levels=1:7)), round(mean(triplets.res[, 10]), 4), round(sd(triplets.res[, 10]), 4))
  all.table      <- rbind(singlets.table, doublets.table, triplets.table)
  rownames(all.table) <- c("Singlets", "Doublets", "Triplets")
  colnames(all.table) <- c("Less_1", "1_to_2", "2_to_3", "3_to_4", "4_to_5","5_to_7","Larger_7", "Mean", "SD")
  
  res <- list(Singlets      = round(singlets.res, precision), 
              Doublets      = round(doublets.res, precision), 
              Triplets      = round(triplets.res, precision), 
              Summary.table = round(all.table, precision))
  class(res) <- "MODFIT"
  return(res)
}

# Export data in MODFIT friendly format ----
Export.MODFIT <- function(data, C, IP, file.name = "MyData") {
  # Missing values: NA -> 9
  data[is.na(data)] <- 9
  write.xlsx2(data, paste0(file.name, "SCORES.xlsx"), col.names = FALSE, row.names = FALSE)
  write.xlsx2(cbind(IP$alpha, IP$delta, IP$taus[, 1:C]), paste0(file.name, "IPs.xlsx"), col.names = FALSE, row.names = FALSE)
}
