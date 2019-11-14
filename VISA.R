

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

visa.norm.100k <- function(x){
    if(sum(x)!=0){
        y=x/sum(x)*100000
        }else{
        y=x
        }
    return(y)
    }


.getChr <-function(x){
    y=unlist(strsplit(x, ":"))
    y=y[1]
    return(y)
    }

.getStart <-function(x){
    y=unlist(strsplit(x, ":"))
    y=y[2]
    y=unlist(strsplit(y, "-"))
    y=y[1]
    return(y)
    }

.getEnd <-function(x){
    y=unlist(strsplit(x, ":"))
    y=y[2]
    y=unlist(strsplit(y, "-"))
    y=y[2]
    return(y)
    }

visa.signal2bdg <- function(x){
    CN=names(x)
    CHR=apply(matrix(CN,ncol=1),1,.getChr )
    START=apply(matrix(CN,ncol=1),1,.getStart )
    END=apply(matrix(CN,ncol=1),1,.getEnd )
    SIGNAL=round(x,2)
    OUT=cbind(CHR,START,END,SIGNAL)
    return(OUT)
    }



visa.col <- function(TAG){
    TAG=as.factor(TAG)
    require(scales)
    my_color_palette <- hue_pal()(length(unique(TAG)))
    COL=my_color_palette[TAG]
    return(COL)
    }

visa.3d <- function(VEC, COL){
    VEC=VEC
    COL=COL
    library("rgl")
    library("car")
    scatter3d(VEC[,1], VEC[,2], VEC[,3], point.col = COL, surface=FALSE)  
    }


