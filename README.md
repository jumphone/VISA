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

    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

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


    COMKC=rep(NA,length(KC$cluster))
    COMKC[which(KC$cluster %in% c(21,2,10,30,28,15,12,29,7,6,3,26))]='C1'
    COMKC[which(KC$cluster %in% c(25,1,11,5))]='C2'
    COMKC[which(KC$cluster %in% c(13,16,9,27,4,20,24,22,23,18,19,14,8,17))]='C3'
    
    Idents(pbmc)=as.character(COMKC)
    DimPlot(pbmc, reduction.use='umap', pt.size=0.1,label=T) 
    
<img src="https://github.com/jumphone/VISA/raw/master/img/VISA3.png" width="500">

    
    ATAC.C1.CN=colnames(pbmc)[which(pbmc@meta.data$batch=='ATAC' & Idents(pbmc)=='C1')]
    ATAC.C1.BC=visa.getBarcode(ATAC.C1.CN)    
    peaks.BC=visa.getBarcode(colnames(peaks))
    peaks.C1=peaks[,which(peaks.BC %in% ATAC.C1.BC)]
    peaks.C1.signal=apply(peaks.C1,1,sum)
    peaks.C1.signal.normal=peaks.C1.signal/length(ATAC.C1.BC)
    
    BDG=visa.signal2bdg(peaks.C1.signal.norm)
    write.table(BDG,file='ATAC.C1.bedgraph',sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
    
    #Then, load data to IGV
  
<img src="https://github.com/jumphone/VISA/raw/master/img/VISA4.png" width="500">

    
    
