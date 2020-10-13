

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
                        x=colnames(VEC)[1],
                        y=colnames(VEC)[2]),
                   theme_bw()
                   )
    #############################################
    rrr.x <-CLUST.VEC
    rrr.df = data.frame(rrr.x); colnames(rrr.df) = c("x","y")
    ##############################################
    OUT=ggplot(data=df,aes(x,y)) + 
        geom_point(colour='grey60',size=0.8)+
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
                        x=colnames(VEC)[1],
                        y=colnames(VEC)[2]),
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




#################################
#2020.10.13
#fastDM

#########################
fastDM.logNorm <- function(x){
    y=x
    y[which(x<0)]=0
    y=x/(sum(x)+1)
    y=y*1000000
    y=log(y+1,10)
    return(y)
    }
#########################

####################################
fastDM.simple_combine <- function(exp_sc_mat1, exp_sc_mat2){    
        exp_sc_mat=exp_sc_mat1
        exp_ref_mat=exp_sc_mat2
        exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
        exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
        gene_sc=rownames(exp_sc_mat)
        gene_ref=rownames(exp_ref_mat)
        gene_over= gene_sc[which(gene_sc %in% gene_ref)]
        exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
        exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
        colname_sc=colnames(exp_sc_mat)
        colname_ref=colnames(exp_ref_mat)
        OUT=list()
        OUT$exp_sc_mat1=exp_sc_mat
        OUT$exp_sc_mat2=exp_ref_mat
        OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
        return(OUT)
        }  
        

########################################################################
fastDM.generate_agg <- function(exp_sc_mat, TAG, print_step=100){
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
            this_new_ref=apply(exp_sc_mat[,this_col],1,sum)
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

########################################################################
fastDM.generate_mean <- function(exp_sc_mat, TAG, print_step=100){
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
    
#########################

fastDM.adjustEXP<- function(EXP, GS){
    library(gmodels)
    VAR=apply(EXP,1,var)
    EXP=EXP[which(VAR>0),]
    N.EXP=apply(EXP,2,fastDM.logNorm)
    S.EXP=t(apply(N.EXP,1,scale))
    rownames(S.EXP)=rownames(EXP)
    colnames(S.EXP)=colnames(EXP)
    EXP.sum=apply(EXP,2,sum)
    RATIO.MAT=c()
    i=1
    while(i<=length(GS)){
        this_gs=GS[[i]]
        this_ratio=apply(EXP[which(rownames(EXP)%in% this_gs),],2,sum)/EXP.sum
        RATIO.MAT=cbind(RATIO.MAT, this_ratio)
        i=i+1}
    
    #################
    R.EXP=S.EXP
    i=1
    while(i<=nrow(EXP)){
        this_exp=S.EXP[i,]
        this_data=cbind(this_exp, EXP.sum, RATIO.MAT)
        colnames(this_data)=c('SE','SUM',paste0('R',1:ncol(RATIO.MAT)))
        this_data=as.data.frame(this_data)
        fit=lm(SE~.,data=this_data)
        out_exp=this_exp-predict(fit)
        R.EXP[i,]=out_exp
        if(i %%1000==1){print(i)}
        i=i+1}
    ###################
    rownames(R.EXP)=rownames(EXP)
    colnames(R.EXP)=colnames(EXP)
    return(R.EXP)
    }
###############

fastDM.smoothEXP <- function(EXP, VEC, CUT=0.1, COR.N=100, SEED=123){
    EXP=EXP
    VEC=VEC
    SEED=SEED
    CUT=CUT
    COR.N=COR.N
    ###############
    set.seed(SEED)
    library(umap)
    custom.config = umap.defaults
    custom.config$n_neighbors=15
    custom.config$n_components=1
    custom.config$spread=1
    custom.config$min_dist=0.1
    UMAP1.OUT=umap(VEC,config =custom.config,method='umap-learn')
    UMAP1=UMAP1.OUT$layout
    #plot(UMAP1,pch=16,col='grey70',cex=0.5)
    U.EXP=EXP
    O.UMAP1=order(UMAP1[,1])
    i=1
    while(i<=nrow(EXP)){
        this_smooth=smooth.spline(EXP[i,O.UMAP1])
        U.EXP[i,O.UMAP1[this_smooth$x]]=this_smooth$y
        if(i %%1000==1){print(i)}
        i=i+1}

    #plot(UMAP1[,1],U.EXP[1,],pch=16,col='grey70',cex=0.5)
    SU.EXP=t(apply(U.EXP,1,scale))
    rownames(SU.EXP)=rownames(U.EXP)
    colnames(SU.EXP)=colnames(U.EXP) 
    #######################
    set.seed(SEED)
    USED.INDEX=sample(1:ncol(EXP),COR.N)
    SCOR=c()
    i=1
    while(i<=nrow(SU.EXP)){
        this_scor=cor(EXP[i,USED.INDEX],SU.EXP[i,USED.INDEX],method='spearman')
        SCOR=c(SCOR,this_scor)
        if(i%%1000==1){print(i)}
        i=i+1}
    #######################
    F.SU.EXP=SU.EXP[which(SCOR>CUT),]
    ######################
    RESULT=list()
    RESULT$smoothedEXP=SU.EXP
    RESULT$filteredSmoothedEXP=F.SU.EXP
    RESULT$spearmanCor=SCOR
    #plot(UMAP1[,1],SU.EXP[1,],pch=16,col='grey70',cex=0.5)
    return(RESULT)
    }



fastDM.smoothPosEXP <- function(EXP, VEC, CUT=0.1, COR.N=100, SEED=123){
    EXP=EXP
    VEC=VEC
    SEED=SEED
    CUT=CUT
    COR.N=COR.N
    ###############
    set.seed(SEED)
    library(umap)
    custom.config = umap.defaults
    custom.config$n_neighbors=15
    custom.config$n_components=1
    custom.config$spread=1
    custom.config$min_dist=0.1
    UMAP1.OUT=umap(VEC,config =custom.config,method='umap-learn')
    UMAP1=UMAP1.OUT$layout
    #plot(UMAP1,pch=16,col='grey70',cex=0.5)
    U.EXP=EXP
    O.UMAP1=order(UMAP1[,1])
    i=1
    while(i<=nrow(EXP)){
        this_ordered_exp=EXP[i,O.UMAP1]
        this_smooth=smooth.spline(this_ordered_exp)
        this_y=this_smooth$y
        this_y[which(this_y<=0)]=0
        this_y[which(this_ordered_exp<=0)]=0
        U.EXP[i,O.UMAP1[this_smooth$x]]=this_y
        if(i %%1000==1){print(i)}
        i=i+1}

    #plot(UMAP1[,1],U.EXP[1,],pch=16,col='grey70',cex=0.5)
    #SU.EXP=t(apply(U.EXP,1,scale))
    SU.EXP=U.EXP
    rownames(SU.EXP)=rownames(U.EXP)
    colnames(SU.EXP)=colnames(U.EXP) 
    #######################
    set.seed(SEED)
    USED.INDEX=sample(1:ncol(EXP),COR.N)
    SCOR=c()
    i=1
    while(i<=nrow(SU.EXP)){
        this_scor=cor(EXP[i,USED.INDEX],SU.EXP[i,USED.INDEX],method='spearman')
        SCOR=c(SCOR,this_scor)
        if(i%%1000==1){print(i)}
        i=i+1}
    #######################
    F.SU.EXP=SU.EXP[which(SCOR>CUT),]
    ######################
    RESULT=list()
    RESULT$smoothedEXP=SU.EXP
    RESULT$filteredSmoothedEXP=F.SU.EXP
    RESULT$spearmanCor=SCOR
    #plot(UMAP1[,1],SU.EXP[1,],pch=16,col='grey70',cex=0.5)
    return(RESULT)
    }



###########################












