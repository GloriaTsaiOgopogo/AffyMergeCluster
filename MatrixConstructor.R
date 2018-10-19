### install packages needed

required_pkgs <- c("dplyr","tibble")
new <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if(length(new) > 0){
  install.packages(new, dependency = T)
}


pkgs_bio <- c("oligo","hgu133plus2cdf")
new <- pkgs_bio[!(pkgs_bio %in% installed.packages()[, "Package"])]
if(length(new) > 0){
  source("http://bioconductor.org/biocLite.R")
  biocLite(new)
}

### Select the path of AffyMergeCluster folder.

if (interactive())
    setwd(choose.dir(getwd(), "Please choose AffyMergeCluster folder."))

### input: A directory contains CEL.gz files(downloaded from NCBI GEO)

require(oligo)

read.cel.files <- function(celPath){

  CEL.list <- list.files(celPath, pattern = "\\.[c|C][e|E][l|L].[g|G][z|Z]$", full.name=TRUE)

  if(length(CEL.list) <= 0){
    cat("Make sure that all the CEL files are of the same platform/array type.\n\n")
  }else{
    cat("If this is the first time processing this platform, it takes some time to install the required package. \n\n")
    CEL.raw <- read.celfiles(filenames = file.path(CEL.list), checkType = F)  ## takes some time

    gpl <- CEL.raw@annotation   ### platfrom eg. "pd.hg.u133.plus.2"
    if(!(gpl %in% installed.packages())){
      source("http://bioconductor.org/biocLite.R")
      biocLite(gpl)
    }
    CEL.raw
  }
}

##supported_gpl <- read.table("./gpl_tab/gpl_map", sep = "\t", header = T)

## Convert expression raw data into present/absent(value 1/0) binary matrix

binary.expr <- function(CEL.raw){

  typeofarray <- class(CEL.raw)[1]

  if(typeofarray == "GeneFeatureSet"){  ## slow
    cat(paste("This is a ",
              typeofarray,
              ". \nIt takes more time than an ExpressionFeatureSet to process.
	\nprocessing...\n\n",
              sep = ""))
    PAcalls <- paCalls(CEL.raw, method = "PSDABG", verbose = FALSE)  ## faster
    ##PAcalls <- ifelse(PAcalls<=0.06,ifelse(PAcalls < 0.04,"P","M"),"A")
    PAcalls <- ifelse(PAcalls < 0.04, 1, 0)  ## 参数调整
    storage.mode(PAcalls) <- "integer"
  }else if(typeofarray == "ExpressionFeatureSet"){
    cat(paste("This is an ",
              typeofarray,
              ". \nprocessing...\n\n",
              sep = ""))
    PAcalls <- paCalls(CEL.raw, method = "MAS5", verbose = FALSE)  ## faster   ## 参数调整
    PAcalls <-PAcalls$calls
    PAcalls <- ifelse(PAcalls == "P", 1, 0)
    storage.mode(PAcalls) <- "integer"	
    PAcalls <- as.data.frame(PAcalls)
	
	cat('0-1 probe matrix has been built.\n\n') 
	
  }else{
    cat(paste("This is a(an) ",
              typeofarray,
              ". \nOnly Gene/Expression-FeatureSet can be processed.\n\n",
              sep = ""))
  }  
  PAcalls
}

## Summary probesets to PDO
require(dplyr)
require(tibble)

probesetID2PDO <- function(PAcalls, gpl){
  load("./AffyData.rda")
  
  if(gpl %in% supported_gpl$Platform_design)
  {
    tmp_data <- get(ls(pattern = gpl))

    probe_matrix <- as_tibble(PAcalls[tmp_data$ProbeID,])

    probe_matrix$pdo <- tmp_data$PDO

    combinecalls <- function(x) ifelse(sum(x)==0,0,1)

    pdo_matrix <- probe_matrix %>% group_by(pdo) %>% summarise_all(funs(combinecalls)) 

    pdo_matrix	<- as.data.frame(pdo_matrix)
	
	pdo_matrix

  }else(
  
    cat("Please check your data. \nType \"supported_gpl\" to see if your data is supported by AffyMergeCluster. \n\n")
  )
    
}

matrix.shrinker <- function(pdo_matrix){

  tmp_matrix <- pdo_matrix[, !(colnames(pdo_matrix) == "pdo")]
  
  tmp_matrix <- tmp_matrix[which(rowMeans(tmp_matrix)!=0),]

  tmp_matrix <- tmp_matrix[which(rowMeans(tmp_matrix)!=1),]
  
  tmp_matrix
}

rename.phy <- function(taxnames){

  ## ensure that the length of sample names are 10 characters.
  taxnames <- substr(taxnames,1,10)
  
  for(i in length(taxnames)){
    while(nchar(taxnames[i]) < 10){
        taxnames[i]<- paste("_",taxnames[i] ,sep = "")
    }
  }
  
  ## deal with duplicated sample names.  
  if(length(unique(taxnames)) < length(taxnames)){
	cat("Duplicated sample names found! \n\n")		
  }  
  sub <-c(LETTERS,seq(0,9,1))
  
  while(!identical(unique(taxnames),taxnames)){
    for(i in which(duplicated(taxnames))){
	  taxnames[i]<- substr(taxnames[i], 1, 7)
	  while(nchar(taxnames[i])<10){
	    taxnames[i] <- paste(taxnames[i], sample(sub, 1),sep = "")
	  }
    }	
  }
  taxnames
}

write.phy <- function(bimatrix, taxnames, outFile){
  if(is.vector(taxnames)){
    cat("Input sample names are:\n  ")
	cat(taxnames)
	cat("\n\n")
	
    if(length(taxnames)!=ncol(bimatrix)){
      cat("Wrong number of sample names. \n\nSet sample names to CEL file names.\n\n")
	  taxnames <- colnames(bimatrix)	  
	}
  }else{
    cat("Wrong input sample names.\n\nSet sample names to CEL file names.\n\n");
	taxnames <- colnames(bimatrix)
  }
  
  taxnames <- rename.phy(taxnames)  
  
  cat("Actual sample names are:\n  ")
  cat(taxnames)
  cat("\n\n")  
  
  cat(paste(dim(bimatrix)[2], dim(bimatrix)[1], sep = "\t"), append = TRUE, file = outFile, sep = "\n")
  for(i in 1:dim(bimatrix)[2]){
    tmp <- character()
    for(j in 1:dim(bimatrix)[1]){	
	  tmp <- paste(tmp, bimatrix[j,i], sep = "\t")
	}	
    cat(paste(taxnames[i], tmp, sep = ""), append = TRUE, file = outFile, sep = "\n")	
  }  
}
