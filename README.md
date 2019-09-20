# VISA

### VISuAlize single cell data 

# Load ATAC Data

    

    library(Seurat)
    source('https://raw.githubusercontent.com/jumphone/VISA/master/VISA.R')

    #load ATAC data from matrix file:
    peaks=visa.read10Xatac(PATH='filtered_peak_bc_matrix')
    
    activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "Mus_musculus.GRCm38.97.chr.gtf",
        seq.levels = c(1:19, "X", "Y"), upstream = 2000, verbose = TRUE)
    
# Visualize Peaks

    #setwd('F:/BEER_REVISE/ATAC')
    
    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
    source('https://raw.githubusercontent.com/jumphone/VISA/master/VISA.R')
    
    library(Seurat)
    peaks <- Read10X_h5("../data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
    
    mybeer = readRDS('mybeer.final.RDS')
    
    pbmc <- mybeer$seurat
    PCUSE=mybeer$select
    pbmc=BEER.combat(pbmc) 
    umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
    pbmc@reductions$umap@cell.embeddings=umap
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F) 
    
    
<img src="https://github.com/jumphone/VISA/raw/master/img/VISA1.png" width="500">



    set.seed(123)
    KC=kmeans(pbmc@reductions$umap@cell.embeddings,centers=30)
    Idents(pbmc)=as.character(KC$cluster)
    DimPlot(pbmc, reduction.use='umap', pt.size=0.1,label=T) 
    

<img src="https://github.com/jumphone/VISA/raw/master/img/VISA2.png" width="500">

    
    KCREF=.generate_mean(as.matrix(pbmc@assays$RNA@data), as.character(KC$cluster))  
    library('gplots')
    C=cor(KCREF,method='spearman')
    H=heatmap.2(C,scale=c("none"),dendrogram='both',Colv=TRUE,Rowv=TRUE,trace='none',
        col=colorRampPalette(c('blue','white','red')),margins=c(5,5))
    HC=as.hclust(H$colDendrogram)
    CLUST=cutree(HC,k=3)
    
    RC=rep('black',ncol(C))
    RC[which(CLUST==1)]='red'
    RC[which(CLUST==2)]='blue'
    RC[which(CLUST==3)]='green'
    
    heatmap.2(C,scale=c("none"),dendrogram='both',Colv=TRUE,Rowv=TRUE,trace='none',
        col=colorRampPalette(c('blue','white','red')),RowSideColors=RC,margins=c(5,5))
  
  
  
<img src="https://github.com/jumphone/VISA/raw/master/img/VISA2.2.png" width="400">
    
    
    CELL.CLUST=rep(NA,ncol(pbmc))
    i=1
    while(i<=max(CLUST)){
        CELL.CLUST[which(KC$cluster %in% names(CLUST[which(CLUST==i)]))]=i
        i=i+1
        }
        
    Idents(pbmc)=as.character(CELL.CLUST)
    DimPlot(pbmc, reduction.use='umap', pt.size=0.1,label=T) 

     
    
<img src="https://github.com/jumphone/VISA/raw/master/img/VISA3.png" width="500">

    
    ATAC.C1.CN=colnames(pbmc)[which(pbmc@meta.data$batch=='ATAC' & Idents(pbmc)=='1')]
    ATAC.C1.BC=visa.getBarcode(ATAC.C1.CN)    
    peaks.BC=visa.getBarcode(colnames(peaks))
    peaks.C1=as.matrix(peaks[,which(peaks.BC %in% ATAC.C1.BC)])
    
    
    peaks.C1.percent=peaks.C1
    peaks.C1.percent[which(peaks.C1>0)]=1
    
    peaks.C1.signal=round(apply(peaks.C1,1,mean)*100)
       
    BDG=visa.signal2bdg(peaks.C1.signal)
    
    
    
    write.table(BDG,file='ATAC.C1.bedgraph',sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
    
    #Then, load data to IGV
 
# IGV

<img src="https://github.com/jumphone/VISA/raw/master/img/VISA4.png" width="1000">

 
 
 



    CHRs=unique(BDG[,1])
    i=1
    while(i<=length(CHRs)){
        this_chr=CHRs[i]
        this_index=which(BDG[,1] == this_chr)
        this_n=min(1000,length(this_index))
        this_start=as.numeric(BDG[this_index,2])
        this_end=as.numeric(BDG[this_index,3])
        this_signal=as.numeric(BDG[this_index,4])
        this_data=cbind(this_start,this_end,this_signal)
        this_data_scale=apply(this_data,2,scale)
        KM=kmeans(this_data_scale[,c(1:2)],centers = this_n )      
        
        i=i+1
        }
    
