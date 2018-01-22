# Sanity.model(): Sanity check - model ----
Sanity.model <- function(model)
{
  if (!(model %in% c("GUM", "GGUM")))
  {
    stop('Parameter "model" can only be "GUM" or "GGUM". Aborted.')
  }
}

# Sanity.C(): Sanity check - C ----
Sanity.C <- function(C, I)
{
  if (!(is.numeric(C) & all(floor(C) == C)& min(C) >= 1 & (length(C) == 1 | length(C) == I)))
  {
    stop('Parameter "C" must be either a positive integer (thus, the same 
             number of observable response categories is assumed for all items) 
             or a vector of positive integers of length I. Aborted.')
  }
}

# Sanity.Cfixed(): Sanity check - C fixed ----
Sanity.Cfixed <- function(C) 
{
  if (max(C) - min(C) > 0)
  {
    stop('C needs to be constant across items under the GUM. Aborted.')
  }
}

# Sanity.data(): Sanity check - data ----
Sanity.data <- function(data)
{
  if (!is.numeric(data) | min(data, na.rm = TRUE) < 0 | any(floor(data) != data, na.rm = TRUE))
  {
    stop('The data matrix must be numeric with non-negative integer entries only. Aborted.')
  }
}

# Sanity.class(): Sanity check - class ----
Sanity.class <- function(obj)
{
  if (class(obj) != "GGUM")
  {
    warning('Parameter "IP" is not of class "GGUM". Please make sure the item parameters 
                object is correctly created.')
  }
}

# Sanity.params(): Sanity check - Iparameters ----
Sanity.params <- function(alpha, delta, taus, theta, C)
{
  alpha.vec <- is.vector(alpha)
  delta.vec <- is.vector(delta)
  if (alpha.vec & delta.vec) 
  {
    length.vecs <- length(alpha) == length(delta)
  } else length.vecs <- FALSE
  taus.mat <- is.matrix(taus)
  if (taus.mat & length.vecs)
  {
    I          <- length(alpha)
    taus.nrows <- (nrow(taus) == I)
    
  } else taus.nrows <- FALSE
  theta.vec <- is.vector(theta)
  if (length.vecs)
  {
    I     <- length(alpha)
    C.vec <- is.vector(C) & all(floor(C) == C) & min(C) >= 1 & (length(C) == 1 | length(C) == I)
  } else C.vec <- FALSE
  if (C.vec)
  {
    taus.ncols <- ncol(taus) == (2*max(C) + 1)
  } else taus.ncols <- FALSE
  
  if (!alpha.vec) stop('alpha should be a vector. Aborted.')
  if (!delta.vec) stop('delta should be a vector. Aborted.')
  if (!length.vecs) stop('alpha and delta need to be vectors of equal length. Aborted.')
  if (!C.vec) 
  {
    stop('"C" must be either a positive integer (thus, the same 
    number of observable response categories is assumed for all items) 
    or a vector of positive integers of length I (= number items). 
    Aborted.')
  }
  if (!(taus.mat & taus.nrows & taus.ncols)) 
  {
    stop('taus should be a matrix with number rows = number items and 
    number columns = 2*max(C)+1. Also, note the strict structure of the 
    taus parameters (see Roberts et al., 2000). Aborted.')
  }
  if (!theta.vec) stop('theta should be a vector. Aborted.')
}

# Sanity.OS(): Sanity check - OS ----
Sanity.OS <- function()
{
  if (Sys.info()["sysname"] != "Windows")
  {
    warning('
This function was written for Windows users because GGUM2004 is a Windows 
executable. It may be of no use in other platforms (e.g., Linux, Mac) unless 
some sort of emulator is used.')
    
  }
}
