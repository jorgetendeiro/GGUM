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
    stop('C needs to be constant across items under the the GUM. Aborted.')
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