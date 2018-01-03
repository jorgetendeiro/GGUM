#' @title Writes a command file for GGUM2004
#'   
#' @description \code{write.GGUM2004} creates a GGUM2004 command file according 
#' to the test characteristics. The file is saved in the GGUM2004 predefined 
#' installation folder.
#'  
#' @param I The number of items.
#' @param C \eqn{C} is the number of observable response 
#' categories minus 1 (i.e., the item scores will be in the set 
#' \eqn{\{0, 1, ..., C\}}{{0, 1, ..., C}}). It should either be a vector of 
#' \eqn{I} elements or a scalar. In the latter case it is assumed that \eqn{C} 
#' applies to all items.
#' @param cutoff Either a number or a vector of \eqn{I} elements which defines 
#' the cutoff value. Default is 2.
#' @param model A string identifying the model. Possible values are "GUM" or 
#'   "GGUM" (default). 
#' @param cmd.file A character string defining the name to give to the command 
#' file. No file extension is required.
#' @param data.file A character string defining the name of the data file. No 
#' file extension is required.
#' @param dir A character string defining the directory where GGUM2004 is 
#' installed (default: "C:/GGUM2004"). The data file identified by the 
#' \code{data.file} parameter (exported by function 
#' \code{\link[GGUM]{export.GGUM2004}}) should exist in this 
#' directory.
#'   
#' @return A script file is saved in the directory where GGUM2004 is installed.
#' 
#' @section Details:
#' This function prepares a GGUM2004 friendly command script, which 
#' may be used to run the GGUM2004 program (Roberts, Donoghue, & Laughlin, 
#' 2000; Roberts et al., 2006). GGUM2004 may be executed directly or may be 
#' called from R, see \code{\link[GGUM]{run.GGUM2004}}.
#' 
#' By default, both the data file (exported via 
#' \code{\link[GGUM]{export.GGUM2004}}) and the GGUM2004 command 
#' script exported by this function are saved in the directory 
#' \code{C:/GGUM2004}, where GGUM2004 is expected to be installed (unzipped).
#' 
#' @references
#' \insertRef{Robertsetal2000}{GGUM}
#' 
#' \insertRef{Robertsetal2006}{GGUM}
#' 
#' @author Sebastian Castro Alvarez, \email{s.castro.alvarez@student.rug.nl}
#' 
#' @examples
#' 1+1
#' @export
write.GGUM2004 <- function(I, C, cutoff = 2, model = "GGUM", 
                           cmd.file = "cmd", 
                           data.file = "data", 
                           dir = "C:/GGUM2004")
{
  # Sanity check - C:
  Sanity.C(C, I)
  
  # Sanity check - C fixed:
  if (model == "GUM") Sanity.Cfixed(C)
  
  # Sanity check - model:
  Sanity.model(model)
  
  if (model == "GGUM") model <- 8 else model <- 3
  
  if (max(C) == min(C)) {
    ConstantC <- "Y"
    response  <- paste(C[1]+1, "NUMBER OF RESPONSE CATEGORIES")
  } else {
    ConstantC <- "N"
    response  <- C + 1
  }
  
  if (max(cutoff) == min(cutoff)) {
    ConstantCO <- "Y"
    responseCO <- paste(cutoff[1], "RESPONSE CUTOFF")
  } else {
    ConstantCO <- "N"
    responseCO <- cutoff
  }
  
  out<-list(paste(model, "ESTIMATE PARAMETERS OF MODEL", model),
            "N CONSTRAINTS ARE NOT USED",
            "N DO NOT CHANGE THE SIGN OF INITIAL PARAMETER ESTIMATES",
            "30 NUMBER OF QUADRATURE POINTS",
            paste0(dir, "/", data.file, ".txt"),
            paste0("(i4,1x,",I,"i2)"),
            paste(I,"NUMBER OF ITEMS"),
            paste(ConstantC,"IS NUMBER OF CATEGORIES CONSTANT?"),
            response,
            "N DO YOU WANT TO RECODE THE DATA?",
            paste(ConstantCO,"IS RESPONSE CUTOFF CONSTANT?"),
            responseCO,
            "N DISCARD ANY ITEMS",
            "N DISCARD ANY PEOPLE",
            "N SIGNS OF INITIAL LOCATION ESTIMATES NOT MANUALLY ASSIGNED",
            "100 NUMBER OF OUTER CYCLES",
            "50 NUMBER OF INNER CYCLES",
            "50 NUMBER OF FISHER SCORING ITERATIONS FOR THRESHOLDS",
            "50 NUMBER OF FISHER SCORING ITERATIONS FOR DELTAS & ALPHAS",
            "0,001 CRITERION",
            "N WANT TO PLOT",
            "N WANT FIT STATISTICS")
  
  writeLines(unlist(out), con = paste0(dir, "/", cmd.file, ".txt"))
}

#' @title Read GGUM2004 person estimates into R
#'
#' @description \code{read.person.GGUM2004} reads the output file from GGUM2004 
#' with the person parameters. Both the person parameter estimates and their 
#' standard errors are imported into R.
#'
#' @param N Number of persons (rows).
#' @param temp.dir The directory where GGUM2004 saved the output. By default it 
#' is "C:/GGUM2004/TEMPFILE".
#' @param precision Number of decimal places of the results (default = 4).
#' 
#' @return An \eqn{N\times 3}{Nx3} matrix is returned. The first column is the 
#' person ID, the second column has the person parameter estimates, and the 
#' last column has the standard errors.
#' 
#' @section Details:
#' The computations are based on the formulas from Roberts, Donoghue, & 
#' Laughlin (2000).
#' 
#' @references
#' \insertRef{Robertsetal2000}{GGUM}
#' 
#' @author Sebastian Castro Alvarez, \email{s.castro.alvarez@student.rug.nl}
#' 
#' @examples
#' 1+1
#' @export
read.person.GGUM2004 <- function(N, temp.dir = "C:/GGUM2004/TEMPFILE", 
                                 precision = 4)
{
  GGUM.file <- readLines(paste0(temp.dir, "/", "FT17F001"))
  GGUM.file <- gsub("=", "= ", GGUM.file)
  GGUM.file <- gsub("#", "", GGUM.file)
  GGUM.file <- gsub("************", "NaN", GGUM.file, fixed = TRUE)
  th        <- grep(pattern = "THETA", GGUM.file)
  theta     <- read.table(textConnection(GGUM.file[th]))
  
  out           <- matrix(NA, N, 3)
  out[, 1]      <- 1:N
  colnames(out) <- c("Person", "Theta","Theta.SE")
  for(i in 1:N){
    if (i %in% theta$V2)
    {
      pos         <- which(theta$V2 == i)
      out[i, 2:3] <- unlist(theta[pos, c(4, 6)])
    }
  }
  
  return(round(out, precision))
}

#' @title Read GGUM2004 Item Estimates
#'
#' @description \code{read.item.GGUM2004} reads the item parameter file from
#'   GGUM2004. It reads the delta estimates, the alpha estimates, the taus
#'   estimates, and their standard errors.
#'
#' @param I the number of items.
#' @param C either a number or a vector. C is the number of observable response
#'   categories minus 1.
#' @param model A string identifying the model. Possible values are "GUM" or 
#' "GGUM" (default).
#' @param temp.dir the path where GGUM2004 save its output.By default it is
#'   "C:/GGUM2004/TEMPFILE".
#' @return \code{read.item.GGUM2004} returns a list cointaning the following
#'   components: \itemize{ \item \code{delta} a vector with delta estimates
#'   \item \code{alpha} a vector with alpha estimates \item \code{taus} a matrix
#'   with taus estimates \item \code{deltaSE} a vector with delta standard
#'   errors \item \code{alphaSE} a vector with alpha standard errors \item
#'   \code{tausSE} a matrix with taus standard errors }
#' @export
read.item.GGUM2004 <- function(I, C, model = "GGUM", 
                             temp.dir="C:/GGUM2004/TEMPFILE")
{
  Sanity.model(model)
  Sanity.C(C, I)
  if(model == "GUM"){ Sanity.Cfixed(C)}
  
  
  
  olddir <- getwd()
  setwd(temp.dir)
  # Clean up of GGUM2004 output item's parameters file
  GGUM.file<-readLines("FT16F001")
  GGUM.file<-gsub("=", "= ", GGUM.file)
  GGUM.file<-gsub("************", "NaN", GGUM.file,fixed=TRUE) # Replace NA when GGUM2004 don't converge
  writeLines(GGUM.file,"FT16F001_Copy")
  it.est<-read.table("FT16F001_Copy",skip=4,header=FALSE,comment.char ="",fill=TRUE)
  setwd(olddir)
  it.est<-split(it.est,it.est$V1)
  delta<-it.est$`ITEM#`[,6]
  deltaSE<-it.est$`ITEM#`[,8]
  if(model == "GGUM"){
    alpha<-it.est$`ITEM#`[,10]
    alphaSE<-it.est$`ITEM#`[,12]
  }
  if(length(C) == 1){C <- rep(C,I)}
  Cplus1 <- C+1
  threshold <- it.est$`THRESHOLD#`[,6]
  thresholdSE <- it.est$`THRESHOLD#`[,8]
  if(model == "GUM"){
    taushalf <- matrix(rep(threshold[2:max(Cplus1)], I), nrow = I, ncol = max(C),
                       byrow = TRUE)
    taushalfSE <- matrix(rep(thresholdSE[2:max(Cplus1)], I), nrow = I, ncol = max(C),
                         byrow = TRUE)
    taus <- cbind(taushalf, rep(0, I), taushalf[, max(C):1]*(-1))
    tausSE <- taushalfSE[, max(C):1]
    
    return(list("delta"=delta,"taus"=taus,
                "deltaSE"=deltaSE,"tausSE"=tausSE))
  }
  else
  {
    taushalf<- matrix(NA, nrow = I, ncol = max(C))
    taushalfSE<-taushalf
    taushalf[1,] <- c(rep(0, max(C)- C[1]), threshold[2:Cplus1[1]])
    taushalfSE[1,] <- c(rep(0, max(C)- C[1]), thresholdSE[2:Cplus1[1]])
    for(i in 2:I){
      taushalf[i,] <- c(rep(0, max(C)- C[i]), threshold[(sum(Cplus1[1:(i-1)]) + 2):sum(Cplus1[1:i])])
      taushalfSE[i,] <- c(rep(0, max(C)- C[i]), thresholdSE[(sum(Cplus1[1:(i-1)]) + 2):sum(Cplus1[1:i])])
    }
    taus <- cbind(taushalf, rep(0, I), taushalf[, max(C):1]*(-1))
    tausSE <- taushalfSE[, max(C):1]
    
    return(list("delta"=delta,"alpha"=alpha,"taus"=taus,
                "deltaSE"=deltaSE,"alphaSE"=alphaSE,"tausSE"=tausSE))
  }
}


#' @title Runs GGUM2004
#'
#' @description \code{run.GGUM2004} sends the command file to GGUM2004 and runs
#'   it. It returns the time and the parameter estimates.
#'
#' @param cmd.file the name of the GGUM2004 command file.
#' @param data.file A character string defining the name of the data file. No 
#' file extension is required.
#' @param dir the directory of GGUM2004 program. It is predefined as
#'   "C:/GGUM2004".
#' @param model A string identifying the model. Possible values are "GUM" or 
#' "GGUM" (default).
#' @return \code{run.GGUM2004} returns a list cointaning the following
#'   components: \itemize{ \item \code{time} the spent time in the analysis
#'   \item \code{delta} a vector with delta estimates \item \code{alpha} a
#'   vector with alpha estimates \item \code{taus} a matrix with taus estimates
#'   \item \code{theta} a vector with theta estimates}
#' @export
run.GGUM2004<-function(cmd.file = "cmd", 
                       data.file = "data", 
                       dir="C:/GGUM2004", 
                       model = "GGUM")
{
  # Prepare dir/TEMPFILE to store the outputs from GGUM2004:
  tempfolder <- paste(dir ,"TEMPFILE", sep="/")
  if (!dir.exists(tempfolder)) dir.create(tempfolder)
  file.remove(dir(tempfolder, full.names = TRUE))
  
  cmd <- paste(paste(dir ,"ggumnsf", sep="/"), 
               paste0(dir, "/", cmd.file, ".txt"), 
               tempfolder)
  
  # Run GGUM2004:
  cat("\n")
  t0 <- proc.time()
  system(cmd)
  t1 <- proc.time()
  
  # Retrieve I, N, C from files:
  tmp <- readLines(paste0(dir, "/", data.file, ".txt"))
  I <- (nchar(tmp[1]) - 5) / 2
  N <- length(tmp)
  rm(tmp); 
  tmp <- readLines(paste0(dir, "/", cmd.file, ".txt"))
  pos <- which(grepl("Y IS NUMBER OF CATEGORIES CONSTANT?", tmp))
  if (length(pos) > 0) 
  {
    C <- as.numeric(substr(tmp[pos+1], 1, 1)) - 1
  } else {
    pos1 <- which(grepl("N IS NUMBER OF CATEGORIES CONSTANT?", tmp))
    pos2 <- which(grepl("N DO YOU WANT TO RECODE THE DATA?", tmp))
    C    <- as.numeric(tmp[(pos1+1):(pos2-1)]) - 1
  }
  C.max <- max(C)
  
  # Retrieve item parameters and SEs:
  items <- read.item.GGUM2004(I, C, model, tempfolder)
  theta <- read.person.GGUM2004(N, tempfolder)
  
  if (model == "GUM")
  {
    SE.out <- cbind(items$deltaSE, items$tausSE)
    colnames(SE.out) <- c("SE.delta", paste0("SE.tau", 1:C.max))
    return(list("time" = t1 - t0, "delta" = items$delta,
                "taus" = items$taus, "theta" = theta, "SE" = SE.out))
  }
  
  if (model == "GGUM")
  {
    SE.out <- cbind(items$alphaSE, items$deltaSE, items$tausSE)
    colnames(SE.out) <- c("SE.alpha", "SE.delta", paste0("SE.tau", 1:C.max))
    return(list("time" = t1 - t0, "alpha" = items$alpha, "delta" = items$delta,
                "taus" = items$taus, "theta" = theta, "SE" = SE.out))
    
  }
}

#' @title Exports data in GGUM2004 friendly format
#'
#' @description \code{export.GGUM2004} exports the data from R to a text file 
#' according to the format required by GGUM2004.
#'
#' @param data The R data matrix to be exported. 
#' @param data.file A character string defining the name of the data file. No 
#' file extension is required.
#' @param dir A character string defining the directory where GGUM2004 is 
#' installed (default: "C:/GGUM2004"). 
#' 
#' @section Details:
#' This function exports the R matrix \code{data} in GGUM2004 (Roberts, 
#' Donoghue, & Laughlin, 2000; Roberts et al., 2006) friendly format. This data 
#' file is to be used together with a GGUM2004 command script (or using the GUI 
#' itself, of course). GGUM2004 may be executed directly or may be called from 
#' R, see \code{\link[GGUM]{run.GGUM2004}}.
#' 
#' By default, both the data file exported by this function and the GGUM2004 
#' command script (exported via \code{\link[GGUM]{write.GGUM2004}}) 
#' are saved in the directory \code{C:/GGUM2004}, where GGUM2004 is expected to 
#' be installed (unzipped).
#' 
#' @references
#' \insertRef{Robertsetal2000}{GGUM}
#' 
#' \insertRef{Robertsetal2006}{GGUM}
#' 
#' @author Sebastian Castro Alvarez, \email{s.castro.alvarez@student.rug.nl}
#' 
#' @examples
#' 1+1
#' @export
export.GGUM2004 <- function(data, data.file = "data", dir = "C:/GGUM2004") 
{
  # Sanity check - data:
  Sanity.data(data)
  
  N <- nrow(data)
  if (N > 2000) stop("GGUM2004 is limited to 2000 subjects or less. Aborted")
  
  I <- ncol(data)
  if (I > 100) stop("GGUM2004 is limited to 100 items or less. Aborted")
  
  if (max(data, na.rm = TRUE) > 9) {
    stop("GGUM2004 is limited to 10 response categories per item (coded 0,...,9) or less. Aborted")
  }
  
  
  data.export <- data
  data.export[is.na(data.export)] <- -9
  address <- paste0(dir, "/", data.file, ".txt")
  
  if (file.exists(address)) file.remove(address)
  
  for (n in 1:N) {
    data.str <- gsub(" -", "-", paste(data.export[n, ], collapse = " "))
    if (substr(data.str, 1, 1) == "-") {space <- " "} else {space <- "  "}
    if (n < 10)                  {write.table(paste(c("000", n, space, data.str), collapse = ""), address,
                                              col.names = FALSE, row.names = FALSE, sep = "", append = TRUE, quote = FALSE)}
    if ((n >= 10) & (n < 100))   {write.table(paste(c("00" , n, space, data.str), collapse = ""), address,
                                              col.names = FALSE, row.names = FALSE, sep = "", append = TRUE, quote = FALSE)}
    if ((n >= 100) & (n < 1000)) {write.table(paste(c("0"  , n, space, data.str), collapse = ""), address,
                                              col.names = FALSE, row.names = FALSE, sep = "", append = TRUE, quote = FALSE)}
    if (n >= 1000)               {write.table(paste(c(       n, space, data.str), collapse = ""), address,
                                              col.names = FALSE, row.names = FALSE, sep = "", append = TRUE, quote = FALSE)}
  }
}
