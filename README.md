# VISA


### VIsualize Single-cell Atac data 

    library(Seurat)
    source('https://raw.githubusercontent.com/jumphone/VISA/master/VISA.R')

    #load ATAC data from matrix file:
    peaks=visa.read10Xatac('filtered_peak_bc_matrix')
    
    
    
