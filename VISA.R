

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


###############
#11.14,2019

visa.lcol <- function(TAG){
    TAG=as.factor(TAG)
    require(scales)
    my_color_palette <- hue_pal()(length(unique(TAG)))
    COL=my_color_palette[TAG]
    return(COL)
    }

visa.col=visa.lcol


visa.vcol <- function(VALUE, CV, CN){
    VALUE=VALUE 
    library('circlize')
    CRF=colorRamp2(CV, CN)
    COL=CRF(VALUE)
    return(COL)
    }



visa.plot3d <- function(VEC, COL){
    VEC=VEC
    COL=COL
    library("rgl")
    library("car")
    scatter3d(VEC[,1], VEC[,2], VEC[,3], 
              xlab=colnames(VEC)[1], ylab=colnames(VEC)[2], zlab=colnames(VEC)[3],
              point.col = COL, surface=FALSE)  
    }


visa.id3d <- function(VEC, COL='blue'){
    VEC=VEC
    LABEL=rownames(VEC)
    COL=COL
    library("rgl")
    library("car")
    x=VEC[,1]
    y=VEC[,2]
    z=VEC[,3]
    OUT=Identify3d(x, y, z, axis.scales=TRUE, groups = NULL, labels = LABEL,
           col = 'blue',
           offset = ((100/length(x))^(1/3)) * 0.02)
    return(OUT)
    }



############################################

#2020.09.07

visa.generate_mean <- function(exp_sc_mat, TAG, print_step=100){
    print_step=print_step
    exp_sc_mat=exp_sc_mat
    TAG=TAG
    
    NewRef=matrix(0,ncol=length(unique(TAG)),nrow=nrow(exp_sc_mat))
    
    TAG=as.character(TAG)
    refnames=unique(TAG)
    total_num=length(refnames)
    outnames=c()
    i=1
    while(i<=length(refnames)){
        one=refnames[i]
        this_col=which(TAG==one)
        outnames=c(outnames,one)
        if(length(this_col) >1){   
            #this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
            this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
            }else{
            this_new_ref = exp_sc_mat[,this_col]
            }
        NewRef[,i]=this_new_ref
        if(i%%print_step==1){print(paste0(i,' / ' ,total_num ))}
        i=i+1       
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }


visa.densityClustPlot <- function(VEC, CLUST, SIZE=5){
    library(MASS)
    library(ggplot2)
    ################################
    VEC=VEC
    CLUST=CLUST
    SIZE=SIZE
    ###############################
    CLUST.VEC=t(visa.generate_mean( t(VEC),as.character(CLUST) ))
    ###############################
    n <- 1000
    x <- VEC
    df = data.frame(x); colnames(df) = c("x","y")

    commonTheme = list(labs(color="Density",fill="Density",
                        x="UMAP_1",
                        y="UMAP_2"),
                   theme_bw()
                   )
    #############################################
    rrr.x <-CLUST.VEC
    rrr.df = data.frame(rrr.x); colnames(rrr.df) = c("x","y")
    ##############################################
    OUT=ggplot(data=df,aes(x,y)) + 
        stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
        scale_fill_continuous(low="green",high="red") + 
        ############################
        geom_text(data=rrr.df, size=SIZE, label=rownames(CLUST.VEC), check_overlap = TRUE) +
        #############################
        xlim(c(min(VEC[,1])-1,max(VEC[,1])+1))+ylim(c(min(VEC[,2])-1,max(VEC[,2])+1))+
        guides(alpha="none") +commonTheme
    
    return(OUT)  
    }


visa.densityPlot <- function(VEC){
    library(MASS)
    library(ggplot2)
    ################################
    VEC=VEC
    
    ###############################
    
    ###############################
    n <- 1000
    x <- VEC
    df = data.frame(x); colnames(df) = c("x","y")

    commonTheme = list(labs(color="Density",fill="Density",
                        x="UMAP_1",
                        y="UMAP_2"),
                   theme_bw()
                   )
    #############################################
    #COL='grey40'
    ##############################################
    OUT=ggplot(data=df,aes(x,y)) + 
        stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
        scale_fill_continuous(low="green",high="red") + 
        ############################
        #geom_point(colour=COL,size=0.8)+
        #############################
        xlim(c(min(VEC[,1])-1,max(VEC[,1])+1))+ylim(c(min(VEC[,2])-1,max(VEC[,2])+1))+
        guides(alpha="none") +commonTheme
    
    return(OUT)  
    }





