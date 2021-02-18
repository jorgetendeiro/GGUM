# There are eight possible unfolding models in GGUM2004.
# Here we only implement the two more popular options:
# - GUM  (Model 3)
# - GGUM (Model 8)
# Model 4 is implemented internally to assist with the estimation of suitable 
# initial values for the GGUM.
# 

#' @import stats utils graphics psych abind viridis xlsx
#' @importFrom Rdpack reprompt

# GenData.GGUM ----
#' @title Generate data from the GUM/GGUM
#'   
#' @description \code{GenData.GGUM} generates all model parameters (items and 
#'   persons) and item scores.
#'   
#' @param N Number of persons (rows).
#' @param I Number of items (columns).
#' @param C \eqn{C} is the number of observable response categories minus 1 
#'   (i.e., the item scores will be in the set \eqn{\{0, 1, ..., C\}}{{0, 1, 
#'   ..., C}}). It should either be a vector of \eqn{I} elements or a scalar. In
#'   the latter, case it is assumed that \eqn{C} applies to all items.
#' @param model A string identifying the model. Possible values are "GUM" or 
#'   "GGUM" (default).
#' @param seed An integer, allowing the user to control the generation process 
#'   (for replication purposes).
#'   
#' @return The function returns a list with five elements: \item{alpha.gen}{The 
#'   discrimination parameters.} \item{delta.gen}{The difficulty parameters.} 
#'   \item{taus.gen}{The threshold parameters.} \item{theta.gen}{The person 
#'   parameters.} \item{data}{The (NxI) data matrix. The item scores are coded 
#'   0, 1, ..., C for an item with (C+1) observable response categories.}
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
#'   Parameters \eqn{\alpha_i}{alpha_i} are randomly uniformly drawn from the 
#'   (.5, 2) interval. Parameters \eqn{\delta_i}{delta_i} are randomly drawn 
#'   from the standard normal distribution bounded between \eqn{-2} and 2. The 
#'   threshold parameters are generated following the same procedure of Roberts,
#'   Donoghue, and Laughlin (2002). Finally, the person parameters are randomly 
#'   drawn from the standard normal distribution.
#'   
#'   If \code{model = "GUM"} the data based on the GUM (Roberts and Laughlin, 
#'   1996) model are generated. The GUM is a constrained version of the GGUM, 
#'   where all discrimination parameters are equal to 1 and the item thresholds
#'   are shared by all items.
#'   
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' gen1 <- GenData.GGUM(500, 10, 5, seed = 456)
#' gen1$data      # Retrieve the data.
#' gen1$alpha.gen # The discrimination parameters.
#' 
#' # Generate data based on items varying in the number of observable response categories:
#' gen2 <- GenData.GGUM(500, 5, c(5, 5, 5, 4, 4), seed = 789)
#' 
#' @export
GenData.GGUM <- function(N, I, C, model = "GGUM", seed = 123)
{
    set.seed(seed)
    
    # Sanity check - C:
    Sanity.C(C, I)
    # Sanity check - model:
    Sanity.model(model)
    
    # Alphas:
    if (model == "GGUM") alpha <- round(runif(I, .5, 2), 4) 
    if (model == "GUM")  alpha <- rep(1, I)
    
    # Deltas:
    delta             <- sort(round(rnorm(I, 0, 1), 4))
    delta[delta < -2] <- -2
    delta[delta > 2]  <- 2
    
    if (length(C) == 1) C <- rep(C, I)
    C.max <- max(C)
    
    # Taus (GGUM):
    if (model == "GGUM")
    {
        tau.half      <- matrix(NA, nrow = I, ncol = C.max)
        tau.half[, 1] <- round(runif(I, .4, 1), 4)
        if (C.max >= 2)
        {
            for (i in 2:C.max) 
            {
                tau.half[, i] <- (i <= C) * (tau.half[, i - 1] + .25 + round(rnorm(I, 0, .04), 4))
            }
        }
        taus <- cbind(-tau.half[, C.max:1], 0, tau.half)
    }
    
    # Taus (GUM):
    if (model == "GUM")
    {
        tau.half    <- rep(NA, C.max)
        tau.half[1] <- round(runif(1, .4, 1), 4)
        if (C.max >= 2)
        {
            for (i in 2:C.max) {
                tau.half[i] <- tau.half[i - 1] + .25 + round(rnorm(1, 0, .04), 4)
            }
        }
        taus <- c(0, tau.half)
        taus <- matrix(rep(taus, I), nrow = I, byrow = TRUE)
        for (i in 1:I) {
            if (C[i] < C.max) taus[i, (C[i] + 2):(C.max + 1)] <- 0
        }
        taus <- cbind(-taus[, (C.max + 1):2], taus)
    }
    
    # Thetas:
    theta <- round(rnorm(N, 0, 1), 4)
    
    # Generate data:
    M           <- 2 * C + 1
    probs.array <- array(NA, dim = c(N, I, C.max + 1))
    for (z in 0:C.max)
    {
        probs.array[, , z + 1] <- P.GGUM(z, alpha, delta, taus, theta, C)
    }
    res <- apply(probs.array, 1:2, function(vec) which( rmultinom(1, 1, vec) == 1) - 1)
    
    return(list(
        alpha.gen = alpha, 
        delta.gen = delta, 
        taus.gen  = taus, 
        theta.gen = theta, 
        data      = res))
}

# GGUM.data.condense ----
GGUM.data.condense <- function(data)
{
    # Sanity check - data:
    Sanity.data(data)
        
    I       <- ncol(data)
    idFac   <- as.integer(factor(apply(data, 1, paste, collapse = "")))
    mat.ord <- cbind(data, idFac)[order(idFac), ]
    mat.sub <- mat.ord[c(1, which(diff(sort(idFac)) == 1) + 1), 1:I]
    return(list(data.condensed = cbind(mat.sub, table(idFac)), ind = idFac))
}

# GUM.initprm ----
GUM.initprm <- function(data.condensed, C, threshold = 2)
{
    N              <- sum(data.condensed[, ncol(data.condensed)])
    n.row          <- nrow(data.condensed)
    I              <- ncol(data.condensed) - 1
    data.full      <- matrix(unlist(sapply(1:n.row, function(x)
    {
        rep(data.condensed[x, 1:I], data.condensed[x, I + 1])
    }
    )), ncol = I, byrow = TRUE)
    data.dich      <- matrix(NA, nrow = N, ncol = I)
    threshold.int  <- rep(threshold, I)
    threshold.int[which(C <= 2)] <- 1
    threshold.mat  <- matrix(rep(threshold.int, N), nrow = N, byrow = TRUE)
    data.dich[data.full >= threshold.mat] <- 1
    data.dich[data.full <  threshold.mat] <- 0
    delta.ini      <- rep(NA, I)
    s.vec          <- colSums(data.dich, na.rm = TRUE)
    sh             <- max(s.vec)
    h              <- which(s.vec == sh)[1]
    delta.ini[h]   <- 0
    tau.B          <- -log(sh / (N - sh))
    if (tau.B > 0) {tau.B <- -tau.B}
    if ((tau.B < -3))
    {
        rescale.f <- (N * exp(2) / (1 + exp(2))) / sh
        s.vec     <- rescale.f * s.vec
        sh        <- max(s.vec)
        tau.B     <- -log(sh / (N - sh))
    }
    for (i in (1:I)[-h])
    {
        func <- function(dlt)
        {
            (exp(-dlt - tau.B) + exp(-2 * dlt - tau.B)) / 
                (1 + exp(-dlt - tau.B) + exp(-2 * dlt - tau.B) + 
                     exp(-3 * dlt)) - s.vec[i] / N
        }
        my.root      <- try(uniroot(f = func, interval = c(0, 6), tol = 1e-20)$
                                root, silent = TRUE)
        delta.ini[i] <- if("try-error" %in% class(my.root)) 2 else my.root
    }
    # Delta signs:
    delta.sgn <- sign(principal(data.full, nfactors = 1, rotate = "none", 
                                scores = FALSE)$loadings[, 1])
    delta.ini <- delta.ini * delta.sgn
    # Taus:
    taus.ini <- matrix(0, nrow = I, ncol = 2 * max(C) + 1)
    for (i in 1:I)
    {
        taus.ini[i, (max(C) - C[i] + 1): max(C)] <- 
            (C[i] + 1 - 1:C[i]) * tau.B / (C[i] - threshold.int[i] + 1)
    }
    taus.ini[, (max(C) + 2):(2 * max(C) + 1)] <- -taus.ini[, max(C):1]
    #
    return(list(delta.ini = delta.ini, taus.ini = taus.ini))
}

# P(Xs|theta) ----
Ls   <- function(data.condensed, alpha, delta, taus, nodes, C)
{
    I           <- length(alpha)
    C.max       <- max(C)
    data.s      <- data.condensed[, -(I + 1), drop = FALSE]
    S           <- nrow(data.condensed)
    N.nodes     <- length(nodes)
    probs.array <- array(NA, dim = c(N.nodes, I, C.max + 1))                 # N.nodes x I x (C + 1)
    for (z in 0:C.max)
    {
        probs.array[, , z + 1] <- P.GGUM(z, alpha, delta, taus, nodes, C)
    }
    #
    mat.ind <- t(data.s) + 1
    ind.arr <- array(0, dim = c(I, S, C.max + 1))
    for (i in 1:I) ind.arr[i, , ][cbind(1:S, mat.ind[i, ])] <- 1
    #
    res     <- matrix(NA, nrow = S, ncol = N.nodes)
    for (f in 1:N.nodes)
    {
        probs.array.ext <- aperm(array(rep(probs.array[f, , ], S), 
                                       dim = c(I, C.max + 1, S)), c(1, 3, 2))
        res[, f] <- apply(probs.array.ext * ind.arr + (1 - ind.arr), 2, prod)
    }
    return(res)
}

# P(Xs) ----
P.tilde.s.vec <- function(Ls.mat, weights)
{
    return(as.vector(Ls.mat %*% weights))
}

# P(Zi = z|theta) ----
P.izf   <- function(alpha, delta, taus, nodes, C)
{
    I           <- length(alpha)
    C.max       <- max(C)
    N.nodes     <- length(nodes)
    probs.array <- array(NA, dim = c(N.nodes, I, C.max + 1))
    for (z in 0:C.max)
    {
        probs.array[, , z + 1] <- P.GGUM(z, alpha, delta, taus, nodes, C)
    }
    return(probs.array)
}

# Auxiliary functions ----
# a.arr
a.arr <- function(alpha, delta, taus, nodes, C)
{
    I        <- length(alpha)
    N.nodes  <- length(nodes)
    C.max    <- max(C)
    taus.ext <- cbind(0, taus)[, 1:(C.max + 1), drop = FALSE]
    res      <- array(NA, dim = c(I, C.max + 1, N.nodes))
    for (z in 0:C.max)
    {
        sum.taus <- sapply(1:I, function(i)
        {
            if (((C.max - C[i]) + z+1) <= (C.max + 1)) {
                sum(taus.ext[i, ((C.max - C[i]) + 1):((C.max - C[i]) + z +1), 
                             drop = FALSE])
            } else 0
        })
        vec.ind <- (z <= C)
        for (f in 1:N.nodes)
        {
            res[, z+1, f] <- vec.ind * 
                exp( alpha * ( z * (nodes[f] - delta) - sum.taus ) )
        }
    }
    return(res)
}
# b.arr
b.arr <- function(alpha, delta, taus, nodes, C)
{
    I        <- length(alpha)
    N.nodes  <- length(nodes)
    C.max    <- max(C)
    M        <- 2 * C + 1
    taus.ext <- cbind(0, taus)[, 1:(C.max + 1), drop = FALSE]
    res      <- array(NA, dim = c(I, C.max + 1, N.nodes))
    for (z in 0:C.max)
    {
        sum.taus <- sapply(1:I, function(i)
        {
            if (((C.max - C[i]) + z+1) <= (C.max + 1)) {
                sum(taus.ext[i, ((C.max - C[i]) + 1):((C.max - C[i]) + z +1), 
                             drop = FALSE])
            } else 0
        })
        vec.ind <- (z <= C)
        for (f in 1:N.nodes)
        {
            res[, z+1, f] <- vec.ind * 
                exp( alpha * ( (M - z) * (nodes[f] - delta) - sum.taus ) )
        }
    }
    return(res)
}
# q.mat
q.mat <- function(taus, C)
{
    I        <- length(C)
    C.max    <- max(C)
    taus.ext <- cbind(0, taus)[, 1:(C.max + 1), drop = FALSE]
    res       <- matrix(NA, nrow = I, ncol = C.max + 1)
    for (z in 0:C.max)
    {
        sum.taus <- sapply(1:I, function(i)
        {
            if (((C.max - C[i]) + z +1) <= (C.max + 1)) {
                sum(taus.ext[i, ((C.max - C[i]) + 1):((C.max - C[i]) + z+1), 
                             drop = FALSE])
            } else 0
        })
        res[, z+1] <- sum.taus
    }
    return(res)
}
# g.mat
g.mat <- function(alpha, delta, taus, nodes, C)
{
    res <- apply(a.arr(alpha, delta, taus, nodes, C) + 
                     b.arr(alpha, delta, taus, nodes, C), 
                 c(1, 3), sum, na.rm = TRUE)
    return(res)
}
# t.mat
t.mat <- function(delta, nodes)
{
    I   <- length(delta)
    res <- t(sapply(1:I, function(i) {nodes - delta[i]}))
    return(res)
}

# dP.alpha.arr ----
dP.alpha.arr <- function(alpha, delta, taus, nodes, C)
{
    I       <- length(alpha)
    N.nodes <- length(nodes)
    if (length(C) == 1) {C <- rep(C, I)}
    C.max   <- max(C)
    M       <- 2 * C + 1
    #
    arr.ind <- matrix(0, nrow = I, ncol = C.max + 1)
    for (i in 1:I) arr.ind[i, 1:(C[i] + 1)] <- 1
    arr.ind <- array(rep(arr.ind, N.nodes), dim = c(I, C.max + 1, N.nodes))
    #
    z.arr   <- arr.ind * aperm(array(rep(0:C.max, I * N.nodes), 
                                     dim = c(C.max + 1, I, N.nodes)), c(2, 1, 3))
    t.arr   <- aperm(array(rep(t.mat(delta, nodes), C.max + 1), 
                           dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2))
    g.arr   <- aperm(array(rep(g.mat(alpha, delta, taus, nodes, C), C.max + 1), 
                           dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2))
    q.arr   <- array(rep(q.mat(taus, C), N.nodes), 
                     dim = c(I, C.max + 1, N.nodes))
    #
    a.res   <- a.arr(alpha, delta, taus, nodes, C)
    b.res   <- b.arr(alpha, delta, taus, nodes, C)
    #
    num1    <- a.res * (z.arr * t.arr - q.arr) + b.res * 
        ((M - z.arr) * t.arr - q.arr)
    den1    <- g.arr
    sum.w   <- apply(num1, c(1, 3), sum, na.rm = TRUE)
    sum.w   <- aperm(array(rep(sum.w, C.max + 1), 
                           dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2))
    num2    <- (a.res + b.res) * sum.w
    den2    <- g.arr^2
    return(aperm(num1 / den1 - num2 / den2, c(3, 1, 2)))
}

# dP.delta.arr ----
dP.delta.arr <- function(alpha, delta, taus, nodes, C)
{
    I         <- length(alpha)
    N.nodes   <- length(nodes)
    if (length(C) == 1) {C <- rep(C, I)}
    C.max     <- max(C)
    M         <- 2 * C + 1
    #
    arr.ind <- matrix(0, nrow = I, ncol = C.max + 1)
    for (i in 1:I)
    {
        arr.ind[i, 1:(C[i] + 1)] <- 1
    }
    arr.ind <- array(rep(arr.ind, N.nodes), dim = c(I, C.max + 1, N.nodes))
    #
    alpha.arr <- array(rep(alpha, (C.max + 1) * N.nodes), 
                       dim = c(I, C.max + 1, N.nodes))
    z.arr     <- arr.ind * aperm(array(rep(0:C.max, I * N.nodes), 
                                       dim = c(C.max + 1, I, N.nodes)), c(2, 1, 3))
    g.arr     <- aperm(array(
        rep(g.mat(alpha, delta, taus, nodes, C), C.max + 1), 
        dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2))
    #
    a.res     <- a.arr(alpha, delta, taus, nodes, C)
    b.res     <- b.arr(alpha, delta, taus, nodes, C)
    #
    num1      <- a.res * (-alpha.arr * z.arr) + b.res * 
        (-alpha.arr * (M - z.arr))
    den1      <- g.arr
    sum.w     <- apply(num1, c(1, 3), sum, na.rm = TRUE)
    sum.w     <- aperm(array(rep(sum.w, C.max + 1), 
                             dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2))
    num2      <- (a.res + b.res) * sum.w
    den2      <- g.arr^2
    return(aperm(num1 / den1 - num2 / den2, c(3, 1, 2)))
}

# dP.tau.arr ----
dP.tau.arr <- function(alpha, delta, taus, nodes, C)
{
    I         <- length(alpha)
    N.nodes   <- length(nodes)
    if (length(C) == 1) {C <- rep(C, I)}
    C.max     <- max(C)
    #
    arr.ind <- matrix(0, nrow = I, ncol = C.max + 1)
    for (i in 1:I)
    {
        arr.ind[i, 1:(C[i] + 1)] <- 1
    }
    arr.ind <- array(rep(arr.ind, N.nodes), dim = c(I, C.max + 1, N.nodes))
    #
    a.res     <- a.arr(alpha, delta, taus, nodes, C)
    b.res     <- b.arr(alpha, delta, taus, nodes, C)
    #
    alpha.arr <- array(rep(alpha, (C.max + 1) * N.nodes), 
                       dim = c(I, C.max + 1, N.nodes))
    z.arr     <- arr.ind * aperm(array(rep(0:C.max, I * N.nodes), 
                                       dim = c(C.max + 1, I, N.nodes)), c(2, 1, 3))
    g.arr     <- aperm(array(
        rep(g.mat(alpha, delta, taus, nodes, C), C.max + 1), 
        dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2))
    #
    res <- array(NA, dim = c(I, C.max + 1, N.nodes, C.max))
    for (k in 1:C.max)
    {
        Uzk.arr <- array(0, dim = c(I, C.max + 1, N.nodes))
        Uzk.arr[ , (k + 1):(C.max + 1), ] <- 1
        sum.w   <- apply(a.res * (-alpha.arr * Uzk.arr) + b.res * 
                             (-alpha.arr * Uzk.arr), c(1, 3), sum, na.rm = TRUE)
        sum.w   <- aperm(array(rep(sum.w, C.max + 1), 
                               dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2))
        num     <- (a.res + b.res) * (g.arr * (-alpha.arr * Uzk.arr) - sum.w)
        den     <- g.arr^2
        res[ , , , k] <- num / den
    }
    return(aperm(res, c(3, 1, 2, 4)))
}

# dP.phi ----
dP.phi <- function(alpha, delta, taus, nodes, C, param = "alphadelta")
{
    if (param == "alphadelta")
    {
        res <- list(alpha = dP.alpha.arr(alpha, delta, taus, nodes, C),
                    delta = dP.delta.arr(alpha, delta, taus, nodes, C),
                    taus  = NULL)
    }
    if (param == "delta")
    {
        res <- list(alpha = NULL,
                    delta = dP.delta.arr(alpha, delta, taus, nodes, C),
                    taus  = NULL)
    }
    if (param == "taus")
    {
        res <- list(alpha = NULL,
                    delta = NULL,
                    taus  = dP.tau.arr(alpha, delta, taus, nodes, C))
    }
    return(res)
}

# DlogL.dphi ----
DlogL.dphi <- function(param = "alphadelta", dP, r.bar.izf, P.izf.arr)
{
    if (param == "alphadelta")
    {
        res.alpha <- apply(r.bar.izf * dP$alpha /  P.izf.arr, 2, sum, 
                           simplify = FALSE, na.rm = TRUE)
        res.delta <- apply(r.bar.izf * dP$delta /  P.izf.arr, 2, sum, 
                           simplify = FALSE, na.rm = TRUE)
        res.taus  <- NULL
    }
    if (param == "delta")
    {
        res.alpha <- NULL
        res.delta <- apply(r.bar.izf * dP$delta /  P.izf.arr, 2, sum, 
                           simplify = FALSE, na.rm = TRUE)
        res.taus  <- NULL
    }
    if (param == "taus")
    {
        res.alpha <- NULL
        res.delta <- NULL
        res.taus  <- apply(r.bar.izf * dP$taus / P.izf.arr, c(2, 4), sum, 
                           simplify = FALSE, na.rm = TRUE)
    }
    res <- list(alpha = res.alpha, delta = res.delta, taus = res.taus)
    return(res)
}

# Estimate thetas and their SEs ----
#' @title Estimate thetas and their SEs (GUM, GGUM)
#'   
#' @description \code{Theta.EAP} estimates the person theta parameters via EAP.
#'   
#' @param IP Object of class \code{GGUM}. The GUM/ GGUM estimated item parameters
#'   via functions \code{GUM()}/ \code{GGUM()}, respectively.
#' @param SE Logical value: Estimate the standard errors of the theta estimates?
#'   Default is \code{TRUE}.
#' @param precision Number of decimal places of the results (default = 4).
#' @param N.nodes Number of nodes for numerical integration (default = 30).
#'   
#' @return If \code{SE = TRUE}, the function returns an \eqn{N\times 2}{Nx2} 
#'   matrix with two columns (thetas, SEs), where \eqn{N}{N} is the number of 
#'   rows in the data matrix (i.e., persons). If \code{SE = FALSE}, the function
#'   returns the theta estimates as a vector of length \eqn{N}{N}.
#'   
#' @section Details: The EAP procedure used here is based on Roberts, Donoghue,
#'   and Laughlin (2000), namely Equation 25 for the \eqn{\theta}{theta}
#'   estimates and Equation 26 for corresponding standard errors. The EAP
#'   estimate is the posterior mean of the \eqn{\theta}{theta} distribution for
#'   the corresponding response pattern. The standard error is computed as an
#'   approximation to the standard deviation of the posterior distribution. See
#'   Roberts et al. (2000) for more details.
#'   
#' @author Jorge N. Tendeiro, \email{tendeiro@hiroshima-u.ac.jp}
#'   
#' @examples
#' # For GUM:
#' # Generate data
#' #   (toy example: Too few items (due to computation time constraints) for 
#' #   accurate estimation of person parameters; larger number of items is 
#' #   required in practice):
#' gen1 <- GenData.GGUM(400, 5, 3, "GUM", seed = 139)
#' # Fit the GUM:
#' fit1 <- GUM(gen1$data, 3)
#' # Estimate the theta parameters:
#' Theta.EAP(fit1)
#' \dontrun{
#' # For GGUM:
#' # Generate data:
#' set.seed(1); C <- sample(3:5, 10, replace = TRUE)
#' gen2 <- GenData.GGUM(2000, 10, C, "GGUM", seed = 156)
#' # Fit the GGUM:
#' fit2 <- GGUM(gen2$data, C)
#' # Estimate the theta parameters:
#' Theta.EAP(fit2)
#' }
#' 
#' @export
Theta.EAP <- function(IP, SE = TRUE, precision = 4, N.nodes = 30)
{
    # Sanity check - class of IP:
    Sanity.class(IP)
    
    data <- IP$data
    I    <- ncol(data)
    N    <- nrow(data)
    C    <- IP$C
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
    ind            <- tmp$ind
    # Nodes and weights:
    nodes   <- seq(-4, 4, length.out = N.nodes)
    weights <- dnorm(nodes) / sum(dnorm(nodes))
    #
    S           <- nrow(data.condensed)
    Ls.mat      <- Ls(data.condensed, IP$alpha, IP$delta, IP$taus, 
                      nodes, C)
    nodes.mat   <- matrix(rep(nodes, S),   ncol = N.nodes, byrow = TRUE) 
    weights.mat <- matrix(rep(weights, S), ncol = N.nodes, byrow = TRUE) 
    num         <- rowSums(nodes.mat * Ls.mat * weights.mat)
    den         <- rowSums(            Ls.mat * weights.mat)
    res         <- num / den
    #
    Th.condensed <- res
    Th.full      <- res[ind]
    Th.full.all  <- rep(NA, N)
    
    if (length(rows.rm) > 0) {
      Th.full.all[-rows.rm] <- Th.full
    } else {
      Th.full.all <- Th.full
    }
    
    if (SE)
    {
      thetas.mat  <- matrix(rep(Th.condensed, N.nodes), ncol = N.nodes,
                            byrow = FALSE)
      num.SE      <- rowSums(((nodes.mat - thetas.mat)^2) * Ls.mat * weights.mat)
      Th.SE.condensed <- sqrt(num.SE / den)
      Th.SE.full  <- sqrt(num.SE / den)[ind]
      # 
      Th.SE.full.all           <- rep(NA, N)
      
      if (length(rows.rm) > 0) {
        Th.SE.full.all[-rows.rm] <- Th.SE.full
      } else {
        Th.SE.full.all <- Th.SE.full
      }
      
      return(cbind(
        Person   = (1:N), 
        Theta    = round(Th.full.all   , precision), 
        Theta.SE = round(Th.SE.full.all, precision)))
    }
    
    return(cbind(
      Person   = (1:N), 
      Theta    = round(Th.full.all , precision)))
}

# d2logP.dtheta2.arr ----
d2logP.dtheta2.arr <- function(data, alpha, delta, taus, theta, C)
{
    I <- length(alpha)
    if (length(C) == 1) {C <- rep(C, I)}
    C.max   <- max(C)
    M       <- 2 * C + 1
    # Nodes (weights seem not needed):
    N.nodes <- 100
    nodes   <- seq(min(-4, min(theta, na.rm = TRUE)), 
                   max(4, max(theta, na.rm = TRUE)), length.out = N.nodes)
    #
    arr.ind <- matrix(0, nrow = I, ncol = C.max + 1)
    for (i in 1:I)
    {
        arr.ind[i, 1:(C[i] + 1)] <- 1
    }
    arr.ind   <- array(rep(arr.ind, N.nodes), dim = c(I, C.max + 1, N.nodes))
    #
    alpha.arr <- array(rep(alpha, (C.max + 1) * N.nodes), dim = c(I, C.max + 1, N.nodes))
    z.arr     <- arr.ind * aperm(array(rep(0:C.max, I * N.nodes), dim = c(C.max + 1, I, N.nodes)), c(2, 1, 3))
    g.arr     <- aperm(array(rep(g.mat(alpha, delta, taus, nodes, C), C.max + 1), dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2))
    #
    a.res     <- a.arr(alpha, delta, taus, nodes, C)
    b.res     <- b.arr(alpha, delta, taus, nodes, C)
    #
    part1     <- (alpha.arr * (z.arr^2) * a.res + alpha * ((M - z.arr)^2) * b.res) / (a.res + b.res)
    part2     <- (z.arr * a.res + (M - z.arr) * b.res) / (a.res + b.res)
    part3     <- (alpha.arr * z.arr * a.res + alpha * (M - z.arr) * b.res) / (a.res + b.res)
    part4     <- apply(part1 * (a.res + b.res), c(1, 3), sum, na.rm = TRUE)
    part4     <- aperm(array(rep(part4, C.max + 1), dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2)) / g.arr
    part5     <- apply(part3 * (a.res + b.res), c(1, 3), sum, na.rm = TRUE)
    part5     <- aperm(array(rep(part5, C.max + 1), dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2)) / g.arr
    part6     <- apply(part2 * (a.res + b.res), c(1, 3), sum, na.rm = TRUE)
    part6     <- aperm(array(rep(part6, C.max + 1), dim = c(I, N.nodes, C.max + 1)), c(1, 3, 2)) / g.arr
    #
    res  <- alpha.arr * (part1 - (part2 * part3) - part4 + (part5 * part6))
    res  <- aperm(res, c(3, 1, 2))
    return(list(d2logP.dtheta2 = res, N.nodes = N.nodes, nodes = nodes))
}

