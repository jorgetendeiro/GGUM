# Define summary() method for class "GGUM":
#' @export
summary.GGUM <- function(object, ...)
{
  x <- object
  I <- ncol(x$data)
  M <- 2*max(x$C) + 1
  
  summ <- cbind(x$alpha, x$delta, x$taus)
  colnames(summ) <- c("Alpha", "Delta", paste0("Tau", 1:M))
  print(summ)
}

# Define print() method for class "GGUM":
#' @export
print.GGUM <- function(x, ...)
{
  I <- ncol(x$data)
  M <- 2*max(x$C) + 1
  
  summ <- cbind(x$alpha, x$delta, x$taus)
  colnames(summ) <- c("Alpha", "Delta", paste0("Tau", 1:M))
  print(summ)
}

# Define plot() method for class "GGUM":
#' @export
plot.GGUM <- function(x, ...)
{
  plotCRC(x)
}

# Define summary() method for class "MODFIT":
#' @export
summary.MODFIT <- function(object, ...)
{
  x <- object
  print(x$Summary.table)
}

# Define print() method for class "MODFIT":
#' @export
print.MODFIT <- function(x, ...)
{
  print(x$Summary.table)
}