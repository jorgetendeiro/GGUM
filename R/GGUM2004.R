#' @title Writes a Generic Command File for GGUM2004
#'
#' @description \code{write.GGUM2004} allows to partially customise the GGUM2004
#'   command file according to the test characteristics. The file is
#'   automatically saved in the GGUM2004 predefined installation folder.
#'
#' @param cmd.file a character string which give a name to the command file.
#' @param data.file a character string with the name of the data file. It is
#'   mandatory to write the extension.
#' @param I the number of items.
#' @param C either a number or a vector. C is the number of observable response
#'   categories minus 1.
#' @param ConstantC indicates whether C is constant or not. Use "Y" if C is
#'   constant, otherwise use "N".
#' @param cutoff either a number or a vector which defines the cutoff value.
#' @param model indicates what unfolding model will be used. It could be "GUM"
#'   or "GGUM".
#' @param cmd.dir the directory for GGUM2004 program. It is predefined as
#'   "C:/GGUM2004".
#' @export

write.GGUM2004<-function(cmd.file,data.file, I, C, ConstantC = c("Y","N"), cutoff = 2,
                         model = c("GUM", "GGUM"), cmd.dir="C:/GGUM2004"){
  if(model == "GUM" && ConstantC == "N"){
    stop("GUM does not allow different number of categories across items. If you choose GUM, ConstantC must be 'N'.")
  }
  if(model == "GUM"){model = 3} else { model = 8}
  if(ConstantC=="Y"){response<-paste(C[1]+1,"NUMBER OF RESPONSE CATEGORIES",sep=" ")}
  if(ConstantC=="N"){response<-C+1}
  ifelse(length(cutoff)>=2,ConstantCO <-"N",ConstantCO<-"Y")
  if(ConstantCO=="Y"){responseCO<-paste(cutoff,"RESPONSE CUTOFF",sep=" ")}
  if(ConstantCO=="N"){responseCO<-cutoff}
  olddir<-getwd()
  out<-list(paste(model, "ESTIMATE PARAMETERS OF MODEL", model, sep = " "),
            "N CONSTRAINTS ARE NOT USED",
            "N DO NOT CHANGE THE SIGN OF INITIAL PARAMETER ESTIMATES",
            "30 NUMBER OF QUADRATURE POINTS",
            paste0("C:\\GGUM2004\\",data.file),
            paste0("(i4,1x,",I,"i2)"),
            paste(I,"NUMBER OF ITEMS",sep=" "),
            paste(ConstantC,"IS NUMBER OF CATEGORIES CONSTANT?",sep=" "),
            response,
            "N DO YOU WANT TO RECODE THE DATA?",
            paste(ConstantCO,"IS RESPONSE CUTOFF CONSTANT?",sep=" "),
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
  setwd(cmd.dir)
  writeLines(unlist(out), con=cmd.file)
  setwd(olddir)}

#' @title Read GGUM2004 Person Estimates
#'
#' @description \code{read.person.GGUM2004} reads the person parameter file from
#'   GGUM2004. It reads the estimates and the standard errors.
#'
#' @param N the number of persons to be read.
#' @param tempfolder the path where GGUM2004 save its output.By default it is
#'   "C:/GGUM2004/TEMPFILE".
#' @return A matrix (dim N*2) is returned with the estimates in the first column
#'   and their standard errors in the second column.
#' @export
read.person.GGUM2004<-function(N,tempfolder="C:/GGUM2004/TEMPFILE"){
  out<-matrix(NA,N,2)
  colnames(out)<-c("theta","thetaSE")
  olddir <- getwd()
  setwd(tempfolder)
  GGUM.file<-readLines("FT17F001")
  GGUM.file<-gsub("=", "= ", GGUM.file)
  GGUM.file<-gsub("#","",GGUM.file)
  GGUM.file<-gsub("************", "NaN", GGUM.file,fixed=TRUE) #To prevent
  th<-grep(pattern= "THETA",GGUM.file)
  theta<-read.table(textConnection(GGUM.file[th]))
  for(i in 1:N){
    if(length(theta$V2[theta$V2==i])!= 0){out[i,]<-c(theta$V4[theta$V2==i],theta$V6[theta$V2==i])}
  }
  setwd(olddir)
  return(out)
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
#' @param model indicates what unfolding model was used. It could be "GUM"
#'   or "GGUM".
#' @param tempfolder the path where GGUM2004 save its output.By default it is
#'   "C:/GGUM2004/TEMPFILE".
#' @return \code{read.item.GGUM2004} returns a list cointaning the following
#'   components: \itemize{ \item \code{delta} a vector with delta estimates
#'   \item \code{alpha} a vector with alpha estimates \item \code{taus} a matrix
#'   with taus estimates \item \code{deltaSE} a vector with delta standard
#'   errors \item \code{alphaSE} a vector with alpha standard errors \item
#'   \code{tausSE} a matrix with taus standard errors }
#' @export
read.item.GGUM2004<-function(I, C, model = c("GUM", "GGUM"), tempfolder="C:/GGUM2004/TEMPFILE"){
  olddir <- getwd()
  setwd(tempfolder)
  # Clean up of GGUM2004 output item's parameters file
  GGUM.file<-readLines("FT16F001")
  GGUM.file<-gsub("=", "= ", GGUM.file)
  GGUM.file<-gsub("************", "NaN", GGUM.file,fixed=TRUE) # Replace NA when GGUM2004 don't converge
  writeLines(GGUM.file,"FT16F001_Copy")
  it.est<-read.table("FT16F001_Copy",skip=4,h=FALSE,comment.char ="",fill=TRUE)
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
    tausSE <- cbind(taushalfSE, rep(0, I), taushalfSE[, max(C):1])

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
    tausSE <- cbind(taushalfSE, rep(0, I), taushalfSE[, max(C):1])

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
#' @param I the number of items.
#' @param C either a number or a vector. C is the number of observable response
#'   categories minus 1.
#' @param N the number of persons to be read.
#' @param model indicates what unfolding model will be used. It could be "GUM"
#'   or "GGUM".
#' @param cmd.dir the directory of GGUM2004 program. It is predefined as
#'   "C:/GGUM2004"
#' @return \code{run.GGUM2004} returns a list cointaning the following
#'   components: \itemize{ \item \code{time} the spent time in the analysis
#'   \item \code{delta} a vector with delta estimates \item \code{alpha} a
#'   vector with alpha estimates \item \code{taus} a matrix with taus estimates
#'   \item \code{theta} a vector with theta estimates}
#' @export
run.GGUM2004<-function(cmd.file, I, C, N, model = c("GUM", "GGUM"), cmd.dir="C:/GGUM2004"){
  tempfolder<-paste(cmd.dir ,"TEMPFILE", sep="/")
  cmd<-paste(paste(cmd.dir ,"ggumnsf", sep="/"), paste(cmd.dir,cmd.file,sep="/"), tempfolder,sep=" ")
  cat("\n")
  t0<-proc.time()
  system(cmd)
  t1<-proc.time()
  items <- read.item.GGUM2004(I, C, model = model, tempfolder = tempfolder)
  theta <- read.person.GGUM2004(N, tempfolder = tempfolder)[,1]
  if (model == GUM){
    return(list("time" = t1 - t0, "delta" = items$delta,
                "taus" = items$taus, "theta" = theta))
  }
  else
  {
    return(list("time" = t1 - t0, "delta" = items$delta, "alpha" = items$alpha,
              "taus" = items$taus, "theta" = theta))
    }
}

#' @title Exports Data to GGUM2004 Format
#'
#' @description \code{export.GGUM2004} writes the data to a text file according
#'   to the format required by GGUM2004.
#'
#' @param data a matrix with the data.
#' @param file.name a character string which give the name to the data file.
#' @export
export.GGUM2004 <- function(data, file.name = "MyData") {
  N <- nrow(data)
  if (N > 2000) {
    stop("GGUM2004 is limited to 2000 subjects or less. Aborted")
  }
  #
  I <- ncol(data)
  if (I > 100) {
    stop("GGUM2004 is limited to 100 items or less. Aborted")
  }
  #
  max.RespCat <- max(data, na.rm = TRUE)
  if (max.RespCat > 9) {
    stop("GGUM2004 is limited to 10 response categories per item (coded 0,...,9) or less. Aborted")
  }
  #
  data.export <- data
  data.export[is.na(data.export)] <- -9
  for (n in 1:N) {
    data.str <- gsub(" -", "-", paste(data.export[n, ], collapse = " "))
    if (substr(data.str, 1, 1) == "-") {space <- " "} else {space <- "  "}
    if (n < 10)                  {write.table(paste(c("000", n, space, data.str), collapse = ""), paste0(file.name, ".txt"),
                                              col.names = FALSE, row.names = FALSE, sep = "", append = TRUE, quote = FALSE)}
    if ((n >= 10) & (n < 100))   {write.table(paste(c("00" , n, space, data.str), collapse = ""), paste0(file.name, ".txt"),
                                              col.names = FALSE, row.names = FALSE, sep = "", append = TRUE, quote = FALSE)}
    if ((n >= 100) & (n < 1000)) {write.table(paste(c("0"  , n, space, data.str), collapse = ""), paste0(file.name, ".txt"),
                                              col.names = FALSE, row.names = FALSE, sep = "", append = TRUE, quote = FALSE)}
    if (n >= 1000)               {write.table(paste(c(       n, space, data.str), collapse = ""), paste0(file.name, ".txt"),
                                              col.names = FALSE, row.names = FALSE, sep = "", append = TRUE, quote = FALSE)}
  }
}
