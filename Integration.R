################################################################################.
#### ----------------      Melanoma gdT - integration     -------------- #######
################################################################################.

###### Clean workspace ######---------------------------------------------------
Sys.setenv(LANG = "en")
rm(list=ls())
gc()

###### Libraries ######---------------------------------------------------------
library(gplots)
library(dplyr)
library(stringr)
library(scales)
library(ggplot2)
library(patchwork)
library(viridis)
library(ccRemover)
library(rgl)
#library(Seurat, lib.loc = "/Library/Frameworks/R.framework/Versions/3.5/Resources/other_libs")
library(Seurat)
library(tidyr)
library(tidyverse)
library(dsb)



###### variables used through script ######-------------------------------------
rm(list=ls())
#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-",""), 3, -1)

# HD
path <- "D:/Sequencing/Melanoma-project"

setwd(path)

data_path <- paste(path, "/data/_preprocessing",sep = "")

results_path <- paste(path, "/data",sep = "")

plots_path <- paste(path, "/data/plots/02.integration",sep = "")




############################      Read data     ################################
################################################################################.

# number of melanoma samples
mel_smpls <- 12

# generate character vector with all the samples except for the 1st one
samples <- c(str_c("10xRNA_H-melanoma-PBMC-gdT",
                   LETTERS[2:mel_smpls],
                   sep = "-"),
             "10xRNA_H-healthy-PBMC-gdT-HD1", # exclude HD1 for now cause seq was not deep enough
             "10xRNA_H-healthy-PBMC-gdT-HD2")

files <- c(str_c(data_path,
                 "/",
                 samples,
                 "_preprocessed_hb.rds"))

# Read in 1st sample manually
gdT.melanoma <- readRDS(paste(data_path,"/10xRNA_H-melanoma-PBMC-gdT-A_preprocessed_hb.rds",
                              sep = ""))

# Read in and merge the rest of the samples
for (i in files) {
    
    x <- readRDS(i)
    
    gdT.melanoma <- merge(gdT.melanoma, y = x)
    
    rm(x)
}

VlnPlot(gdT.melanoma,
        features = "percent_hb",
        group.by = "orig.ident")

gdT.melanoma <- subset(gdT.melanoma,
                       subset = percent_hb < 0.005)

VlnPlot(gdT.melanoma,
        features = "percent_hb",
        group.by = "orig.ident")




#####################      Clusters pre-integration     ########################
################################################################################.

gdT.melanoma <-  NormalizeData(gdT.melanoma,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)

gdT.melanoma <- FindVariableFeatures(gdT.melanoma,
                                     selection.method = "vst",
                                     nfeatures = 2000)

gdT.melanoma <- ScaleData(gdT.melanoma,
                          verbose = TRUE,
                          vars.to.regress = c("nCount_RNA",
                                              "percent_mt",
                                              "G2M.Score",
                                              "S.Score"),
                          do.scale = TRUE)


gdT.melanoma <- RunPCA(gdT.melanoma,
                       npcs = 60,
                       verbose = TRUE)

ElbowPlot(gdT.melanoma, ndims = 60)



gdT.melanoma <- FindNeighbors(gdT.melanoma, dims = 1:25)

gdT.melanoma <- FindClusters(gdT.melanoma, resolution = 0.5, 
                             cluster.name = "unintegrated_clusters")

gdT.melanoma <- RunUMAP(gdT.melanoma,
                        dims = 1:25,
                        reduction.name = "umap.unintegrated")


DimPlot(gdT.melanoma, group.by = c("orig.ident", "seurat_clusters"), 
        label = F, reduction = "umap.unintegrated")

ggsave(paste(dato,"unintegrated.png",sep = "_"),
       path = plots_path,
       width = 21,
       height = 8,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)





##########################      Integration     ################################
################################################################################.
gc()
gdT.melanoma <- IntegrateLayers(object = gdT.melanoma, 
                                method = HarmonyIntegration,
                                orig.reduction = "pca", 
                                new.reduction = "harmony",
                                verbose = T)




#######################      check integration     #############################
################################################################################.
gc()

gdT.melanoma <- FindNeighbors(gdT.melanoma, dims = 1:25,
                              reduction = "harmony")

gdT.melanoma <- FindClusters(gdT.melanoma, resolution = 0.5, 
                             cluster.name = "harmony.clusters")

gdT.melanoma <- RunUMAP(gdT.melanoma,
                        dims = 1:25,
                        reduction = "harmony",
                        reduction.name = "umap.harmony")

DimPlot(gdT.melanoma, group.by = c("orig.ident", "harmony.clusters"), 
        label = F, reduction = "umap.harmony")

ggsave(paste(dato,"harmony.png",sep = "_"),
       path = plots_path,
       width = 21,
       height = 8,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)






##########################      Join Layers     ################################
################################################################################.

# seurat v5 is organized in layers
gdT.melanoma <- JoinLayers(gdT.melanoma)
gc()




#############################     Metadata     #################################
################################################################################.

##### add metadata -------------------------------------------------------------

gdT.melanoma@meta.data <- gdT.melanoma@meta.data %>% 
    mutate(outcome = case_when(str_detect(hash.ID, "CR") ~ "R",                 # basic keywords
                               str_detect(hash.ID, "PD") ~ "NR",
                               str_detect(hash.ID, "BC") ~ "HD"),
           
           timepoint = case_when(str_detect(hash.ID, "Baseline") ~ "Baseline",
                                 str_detect(hash.ID, "Eval") ~ "Eval",
                                 str_detect(hash.ID, "BC") ~ "Healthy")) 





gdT.melanoma@meta.data <- gdT.melanoma@meta.data %>% 
    mutate(outcome.long = case_when(outcome == "R" ~ "Responder",               # combinations
                                    outcome == "NR" ~ "Non-responder",               # combinations
                                    outcome == "HD" ~ "Healthy donor"),
           
           outcome.lane = paste(gdT.melanoma$outcome,
                                str_sub(gdT.melanoma$orig.ident, start = -1),
                                sep = "_"),
           
           treatment = case_when(str_detect(orig.ident, "-HD") ~ "healthy",
                                 str_detect(orig.ident, "-A|-B|-C|-D|-E|-F|-H|-J") ~ "mono",      
                                 str_detect(orig.ident, "-G|-I|-K|-L") ~ "combi"),
           outcome.timepoint = str_c(outcome, "_", timepoint))


gdT.melanoma@meta.data <- gdT.melanoma@meta.data %>% 
    mutate(outcome.treatment = paste(gdT.melanoma$outcome,
                                     treatment,
                                     sep = "_"),
           
           
           group = paste(gdT.melanoma$treatment,
                         gdT.melanoma$outcome,
                         gdT.melanoma$timepoint,
                         sep = "_"))


gdT.melanoma@meta.data <- gdT.melanoma@meta.data %>% 
    mutate(sampleID = case_when(!str_detect(hash.ID, "BC") ~ paste(gdT.melanoma$outcome,                               # sampleID + patient ID
                                                                   gdT.melanoma$timepoint,
                                                                   str_sub(gdT.melanoma$orig.ident, start = -1),
                                                                   sep = "_"),
                                str_detect(hash.ID, "BC01") ~ "Healthy_1",
                                str_detect(hash.ID, "BC09") ~ "Healthy_2",
                                str_detect(hash.ID, "BC15") ~ "Healthy_3",
                                str_detect(hash.ID, "BC16") ~ "Healthy_4",
                                str_detect(hash.ID, "BC18") ~ "Healthy_5",
                                str_detect(hash.ID, "BC23") ~ "Healthy_6",
                                str_detect(hash.ID, "BC25") ~ "Healthy_7",
                                str_detect(hash.ID, "BC27") ~ "Healthy_8"))


################################################################################.
##########################      Save data     ##################################
################################################################################.

saveRDS(object = gdT.melanoma, 
        file = paste(results_path, "/", "10xRNA_H-melanoma-PBMC-gdT_int-hb.rds",sep=""))

