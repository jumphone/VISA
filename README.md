# VISA


### VIsualize Single-cell Atac data 

    library(Seurat)
    source('https://raw.githubusercontent.com/jumphone/VISA/master/VISA.R')

    #load ATAC data from matrix file:
    peaks=visa.read10Xatac(PATH='filtered_peak_bc_matrix')
    
    activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "Mus_musculus.GRCm38.97.chr.gtf",
        seq.levels = c(1:20, "X", "Y"), upstream = 2000, verbose = TRUE)
    
    
    
    
    
