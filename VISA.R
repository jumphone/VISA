

visa.read10Xatac <- function(PATH){
    PATH=PATH
    library(Matrix)
    indata <- Matrix::readMM(paste0(PATH,"/matrix.mtx")) 
    # binarize the matrix
    indata@x[indata@x > 0] <- 1

    # format cell info
    cellinfo <- read.table(paste0(PATH,"/barcodes.tsv"))
    row.names(cellinfo) <- cellinfo$V1
    names(cellinfo) <- "cells"

    # format peak info
    peakinfo <- read.table(paste0(PATH,"/peaks.bed"))
    names(peakinfo) <- c("chr", "bp1", "bp2")
    peakinfo$site_name <- paste0(peakinfo$chr,':', peakinfo$bp1,'-',peakinfo$bp2)
    row.names(peakinfo) <- peakinfo$site_name

    row.names(indata) <- row.names(peakinfo)
    colnames(indata) <- row.names(cellinfo)
    
    return(indata)
    }


.getLast <-function(x){
    y=unlist(strsplit(x, "_"))
    y=y[length(y)]
    return(y)
    }


visa.getBarcode <- function(CN){
    y=apply(matrix(CN,ncol=1),1,.getLast )
    return(y)    
    }

