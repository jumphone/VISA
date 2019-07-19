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


    library(Seurat)
    peaks <- Read10X_h5("../data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
    
    mybeer = readRDS('mybeer.final.RDS')
    
    pbmc = readRDS('') 
    
