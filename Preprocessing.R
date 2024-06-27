################################################################################.
##### --------                Melanoma gdT                     --------- #######
################################################################################.

###### Clean workspace ######---------------------------------------------------
Sys.setenv(LANG = "en")
rm(list=ls())
gc()

###### Libraries ######---------------------------------------------------------
library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)
library(DESeq2)
library(scater)
library(table1)
library(readxl)
library(viridis)


###### variables used through script ######---------------  --------------------
#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-",""), 3, -1)

# working directory, careful, qmd changes wd!!! write it manually
wd <- "D:/Sequencing/Melanoma-project"

# sample
lane <- "L"

xpt <- "10xRNA_H-melanoma-PBMC-gdT" 

sample <- paste(xpt, 
                lane,
                sep = "-") ### sample name (name of folder!!!)


data_path <- str_c(wd,
                   "/_raw/",
                   sample)

results_path <- str_c(wd,
                      "/data/_preprocessing/QC")


final_path <- str_c(wd,
                    "/data/_preprocessing/")



##### #######################   read data    ####################################
################################################################################.

filtered_data <- Read10X(data.dir = str_c(data_path,"/filtered_feature_bc_matrix"))

filtered_data_ADT <- filtered_data$`Antibody Capture`


data <- CreateSeuratObject(counts = filtered_data$`Gene Expression`, 
                           project = sample, 
                           min.cells = 3, 
                           min.features = 100)

data




##### #######################      QC       ####################################
################################################################################.

### Add metadata ---------------------------------------------------------------
# Mitochondrial content
data <- PercentageFeatureSet(data,
                             pattern = "^MT-",
                             col.name = "percent_mt")

# Ribosomal content
data <- PercentageFeatureSet(data,
                             pattern = "^RP[SL]",
                             col.name = "percent_ribo")

# Hemoglobin content
data <- PercentageFeatureSet(data,
                             pattern = "^HB[^(P)]",
                             col.name = "percent_hb")

data("cc.genes")
sgenes <- cc.genes$s.genes
g2mgenes <- cc.genes$g2m.genes

data <- NormalizeData(data)
data <- CellCycleScoring(data,
                         assay = "RNA", 
                         s.features = sgenes, 
                         g2m.features = g2mgenes, 
                         set.ident = FALSE)



##### #######################    GEX QC     ####################################
################################################################################.

## setting thresholds 
{r}
nFeature_max <- 4000  # genes
nFeature_min <- 500
percent_mito_max <- 7 ## 10 is recommended; sometimes ppl go to 5 but be careful
percent_hb_max <- 0.005 ## filter out cells that express hemoglobin genes --> red blood cell contamination
nCount_max <- NULL ## read, you can decide if using it or not
nCount_min <- NULL

#Adjust filtering depending on the thresholds you want to use
dataQC_gex <- subset(data,
                     subset = nFeature_RNA > nFeature_min &
                         nFeature_RNA < nFeature_max & 
                         percent_mt < percent_mito_max # &
                     # percent_hb < percent_hb_max
)


##### #######################    HTO QC     ####################################
################################################################################.

## select only HTO data for UMI present in pat_data ("clean" dataset)
hto <- filtered_data_ADT[, colnames(filtered_data_ADT) %in% colnames(dataQC_gex)]

# add HTO data to seurat object
dataQC_gex[["HTO"]] <- CreateAssayObject(counts = hto)





## Normalize the HTO expression levels.
dataQC_gex <- NormalizeData(dataQC_gex,
                            assay = "HTO", 
                            normalization.method = "CLR")



## Demultiplex cells based on HTO enrichment
# here we use default methods but for very large datasets other parameters might be better
dataQC_gex <- HTODemux(dataQC_gex, 
                       assay = "HTO",
                       positive.quantile = 0.99)


hto_threshold <- "0.99"

## filtering 
Idents(dataQC_gex) <- "HTO_classification.global"

dataQC_gex_hto <- subset(dataQC_gex, idents = "Singlet")

## filter hemoglobin genes
dataQC_gex_hto_hb  <- subset(dataQC_gex_hto,
                             subset = percent_hb < percent_hb_max)




##### #######################     Save      ####################################
################################################################################.

## Make sure that all cells have sample name in metadata
dataQC_gex_hto_hb  <- RenameCells(dataQC_gex_hto, add.cell.id = sample) 

saveRDS(object = dataQC_gex_hto_hb, 
        file = paste(final_path, sample,"_preprocessed_hb.rds",sep=""))
