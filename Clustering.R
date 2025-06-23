################################################################################.
#### --------            Melanoma gdT   -  Clustering          --------- #######
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

# colour string for imputation and overlays
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))


# HD
# path <- "E:/Sequencing/Melanoma-project"


setwd(path)

data_path <- paste(path, "/data",sep = "")
plots_path <- paste(path, "/data/plots/03.cluster-res_tcr-filt",sep = "")
#plots_path <- paste(path, "/data/plots/03.cluster-res",sep = "")




##### #########################   read data    ####################################
################################################################################.

gdT.melanoma <- readRDS(paste(data_path,"/10xRNA_H-melanoma-PBMC-gdT_int-hb.rds",
                              sep = ""))




# ###########################   initial clustering    ##########################
################################################################################.
gc()
# Basic clustering

gdT.melanoma <- FindVariableFeatures(gdT.melanoma,
                                     selection.method = "vst",
                                     nfeatures = 2000) #increased it if necessary


top10 <- head(VariableFeatures(gdT.melanoma))
p1 <- VariableFeaturePlot(gdT.melanoma)
p2 <- LabelPoints(plot = p1,
                  points = top10,
                  repel = T)

p1+p2


gdT.melanoma <- ScaleData(gdT.melanoma,
                          verbose = TRUE,
                          vars.to.regress = c("nCount_RNA",
                                              "percent_mt",
                                              "G2M.Score",
                                              "S.Score"),
                          do.scale = TRUE,
                          features = rownames(gdT.melanoma) # not enough memory to scale everything
                          #features = x #10000 variable features to reduce size of object
) #


gdT.melanoma <- RunPCA(gdT.melanoma,
                       npcs = 60,
                       verbose = TRUE)



ElbowPlot(gdT.melanoma, ndims = 60)
# ggsave(paste(dato,"ElbowPlot_60PCs.png",sep = "_"), 
#        path = plots_path,
#        width = 20, 
#        height = 15, 
#        units = "cm",
#        dpi = 600, 
#        scale=1.5, 
#        bg = "white",
#        limitsize = FALSE)


gdT.melanoma <- FindNeighbors(gdT.melanoma, dims = 1:30,
                              reduction = "harmony")

gdT.melanoma <- FindClusters(gdT.melanoma,
                             resolution = 0.1,
                             cluster.name = "integrated_res.0.1")


# Calculate UMAP again and name it as default
gdT.melanoma <- RunUMAP(gdT.melanoma,
                        dims = 1:30,
                        reduction = "harmony",
                        reduction.name = "umap")

DimPlot(subset(gdT.melanoma, downsample=10000), group.by = "integrated_res.0.1", 
        label = TRUE, reduction = "umap")



## ##################      Filter contaminating cells    ########################
################################################################################.

 
##### ###############        TCRab expressing cells   ###########################

# filter only TRAV and TRAJ genes because gdT can express rearanged beta chain
aTCR_genes <-  c(rownames(gdT.melanoma[["RNA"]]$data)[startsWith(rownames(gdT.melanoma[["RNA"]]$data),"TRAV")],
                 rownames(gdT.melanoma[["RNA"]]$data)[startsWith(rownames(gdT.melanoma[["RNA"]]$data),"TRAJ")])


TCRd<- c("TRDV1",    
         "TRDV2",    
         "TRDV3")

DoHeatmap(gdT.melanoma,
          features = c(TCRd, aTCR_genes))

 


###### unwanted cells------------------------------------------------------------------------
aTCR_genes
aTCR_genes <- aTCR_genes[! aTCR_genes %in% c("TRAV14DV4", # Vd4 cells
                                             "TRAV29DV5", # Vd5 cells
                                             "TRAV23DV6", # Vd6 cells
                                             
                                             "TRAV36DV7",
                                             "TRAV38-2DV8")]
genelist <- aTCR_genes
### init empty list
genelist_subset_query <- ""

### for loop to construct query
for(i in 1:length(genelist)){
    if(i != length(genelist)){
        genelist_subset_query <- paste(genelist_subset_query, "`", genelist[i], "`", "> 0 |", sep = "") #### can change & to | (logical "or") if needed
    }
    if(i == length(genelist)){
        genelist_subset_query <- paste(genelist_subset_query, "`", genelist[i], "`", "> 0", sep = "")
    }
}

genelist_subset_query

abCells <- subset(gdT.melanoma,
                  subset = `TRAV18`> 0 |`TRAV21`> 0 |`TRAV22`> 0 |`TRAV24`> 0 |
                  `TRAV41`> 0 |`TRAV5`> 0 |`TRAV8-4`> 0 |`TRAV8-5`> 0 |
                  `TRAV13-2`> 0 |`TRAV19`> 0 |`TRAV26-1`> 0 |`TRAV27`> 0 |
                  `TRAV38-1`> 0 |`TRAV26-2`> 0 |`TRAV35`> 0 |`TRAV39`> 0 |
                  `TRAV40`> 0 |`TRAV4`> 0 |`TRAV12-2`> 0 |`TRAV13-1`> 0 |
                  `TRAV17`> 0 |`TRAV30`> 0 |`TRAV12-3`> 0 |`TRAV31`> 0 |
                  `TRAV34`> 0 |`TRAV12-1`> 0 |`TRAV1-2`> 0 |`TRAV9-2`> 0 |
                  `TRAJ49`> 0 |`TRAJ46`> 0 |`TRAJ45`> 0 |`TRAJ44`> 0 |
                  `TRAJ42`> 0 |`TRAJ38`> 0 |`TRAJ35`> 0 |`TRAJ58`> 0 |
                  `TRAJ52`> 0 |`TRAJ41`> 0 |`TRAJ39`> 0 |`TRAJ54`> 0 |
                  `TRAJ37`> 0 |`TRAJ43`> 0 |`TRAJ36`> 0 |`TRAJ34`> 0)



##### ###############            TCRVD expression      #########################
delta_genes <-  c(rownames(gdT.melanoma[["RNA"]]$data)[startsWith(rownames(gdT.melanoma[["RNA"]]$data),"TRDV")])

DoHeatmap(gdT.melanoma,
          features = c(TCRd))

VlnPlot(gdT.melanoma,
        features =  delta_genes,
        pt.size = 0.1,
        alpha = 0.1,
        ncol = 2)

TCRd_genes <- c(delta_genes,
                "TRAV14DV4", # Vd4 cells
                "TRAV29DV5", # Vd5 cells
                "TRAV23DV6", # Vd6 cells
                
                "TRAV36DV7",
                "TRAV38-2DV8")


###### no delta Tcr -------------------------------------------------------------

noDelta <- subset(gdT.melanoma,
                  subset = `TRDV1`== 0 &`TRDV2`== 0 &`TRDV3`== 0 &  
                      `TRAV14DV4`== 0 & `TRAV29DV5`== 0 & `TRAV23DV6`== 0 & 
                      `TRAV36DV7`== 0 & `TRAV38-2DV8`== 0)


table(noDelta$integrated_res.0.1)


##### ##############        Identify cycling cells   ###########################


Idents(gdT.melanoma) <- gdT.melanoma@meta.data$integrated_res.0.1
VlnPlot(gdT.melanoma, features =  c("nCount_RNA", "percent_mt","G2M.Score",
                                    "S.Score", "percent_ribo", "percent_hb"),
        pt.size = 0.1,
        alpha = 0.1,
        ncol = 2)

VlnPlot(gdT.melanoma, features =  c("nCount_RNA", "percent_mt","G2M.Score",
                                    "S.Score", "percent_ribo", "percent_hb"),
        pt.size = 0,
        ncol = 2)
ggsave(paste(dato,"gen-vln.png",sep = "_"),
       path = paste0(plots_path, "/cycle-macrophages"),
       width = 9,
       height = 10,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)



VlnPlot(gdT.melanoma, features =  c("S.Score"),
        pt.size = 0.1,
        alpha = 0.1,
        ncol = 2)+
    geom_hline(yintercept=0.2,
               color = "red", 
               size=1)

ggsave(paste(dato,"S-score_vln.png",sep = "_"),
       path = paste0(plots_path, "/cycle-macrophages"),
       width = 9,
       height = 4,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)



FeaturePlot(gdT.melanoma,
            features = c("G2M.Score","S.Score") ,
            reduction = "umap")


###### unwanted cells-----------------------------------------------------------
cyclingCells <- subset(gdT.melanoma,
                       subset = S.Score > 0.2)


 



##### ##############       Monocytes / Macrophages      ########################

# Make sure we are calculating DEGs for the right identities
Idents(gdT.melanoma) <- gdT.melanoma@meta.data$integrated_res.0.1

table(gdT.melanoma@meta.data$integrated_res.0.1)

DEG.pos <- FindAllMarkers(gdT.melanoma,
                          only.pos = T,           # only including UPregulated DEGs per cluster?
                          test.use = "wilcox",    # test to use
                          min.pct = 0.25,         # genes included in results should be expressed in atleast 25% of the cluster
                          logfc.threshold = 0.15 # genes included should be increased more than a logFC 0.15 compared to the other cells
                          #max.cells.per.ident = 5000 # downsample if necessary
) 

top10.pos <- DEG.pos %>% 
    filter(p_val_adj < 0.05) %>% 
    group_by(cluster) %>% 
    top_n(n = 10, wt = avg_log2FC)


b <- DoHeatmap(subset(gdT.melanoma, downsample=400), 
               features = top10.pos$gene) + 
    scale_fill_gradientn(colours = mycols) 

ggsave(paste(dato,"DEG-Heatmap_top10.png",sep = "_"),
       plot = b,
       path = paste0(plots_path, "/cycle-macrophages"),
       width = 15,
       height = 25,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)



VlnPlot(gdT.melanoma, features =  c("FCN1", "VCAN", "MS4A6A"),
        pt.size = 0.3,
        alpha = 0.1,
        ncol = 2)

ggsave(paste(dato,"macrophage_vln.png",sep = "_"),
       path = paste0(plots_path, "/cycle-macrophages"),
       width = 9,
       height = 10,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)



###### unwanted cells-----------------------------------------------------------
macro <- subset(gdT.melanoma,
                subset = FCN1 > 0 & VCAN > 0 & MS4A6A > 0)



##### ##################        Filtering      #################################
ncol(gdT.melanoma)
# before filtering there is 97091 cells in dataset

gdT.melanoma <- subset(gdT.melanoma, # filter only TRA genes because some gd can have re-arranged beta chain
                       subset = `TRAV18`== 0 &`TRAV21`== 0 &`TRAV22`== 0 &`TRAV24`== 0 &
                           `TRAV41`== 0 &`TRAV5`== 0 &`TRAV8-4`== 0 &`TRAV8-5`== 0 &
                           `TRAV13-2`== 0 &`TRAV19`== 0 &`TRAV26-1`== 0 &`TRAV27`== 0 &
                           `TRAV38-1`== 0 &`TRAV26-2`== 0 &`TRAV35`== 0 &`TRAV39`== 0 &
                           `TRAV40`== 0 &`TRAV4`== 0 &`TRAV12-2`== 0 &`TRAV13-1`== 0 &
                           `TRAV17`== 0 &`TRAV30`== 0 &`TRAV12-3`== 0 &`TRAV31`== 0 &
                           `TRAV34`== 0 &`TRAV12-1`== 0 &`TRAV1-2`== 0 &`TRAV9-2`== 0 &
                           `TRAJ49`== 0 &`TRAJ46`== 0 &`TRAJ45`== 0 &`TRAJ44`== 0 &
                           `TRAJ42`== 0 &`TRAJ38`== 0 &`TRAJ35`== 0 &`TRAJ58`== 0 &
                           `TRAJ52`== 0 &`TRAJ41`== 0 &`TRAJ39`== 0 &`TRAJ54`== 0 &
                           `TRAJ37`== 0 &`TRAJ43`== 0 &`TRAJ36`== 0 &`TRAJ34`== 0)

ncol(gdT.melanoma)
# afer TCRa filtering there is 94373 cells left; we filtered out 2718 cells


gdT.melanoma <- subset(gdT.melanoma,
                       subset = S.Score < 0.2)

ncol(gdT.melanoma)
# after S.score filtering we have 94116 cells left; we filtered out 257 cells



gdT.melanoma <- subset(gdT.melanoma,
                       subset = FCN1 == 0 | VCAN == 0 | MS4A6A == 0)

ncol(gdT.melanoma)
# 94070 cells left, 46 filtered

VlnPlot(gdT.melanoma, features =  c("FCN1", "VCAN", "MS4A6A"),
        pt.size = 0.3,
        alpha = 0.1,
        ncol = 2)

ggsave(paste(dato,"macrophage-after-filter_vln.png",sep = "_"),
       path = paste0(plots_path, "/cycle-macrophages"),
       width = 9,
       height = 10,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)



gdT.melanoma <- subset(gdT.melanoma,
                       subset = `TRDV1`> 0 | `TRDV2`> 0 |`TRDV3`> 0 |  
                           `TRAV14DV4`> 0 | `TRAV29DV5`> 0 | `TRAV23DV6`> 0 | 
                           `TRAV36DV7`> 0 | `TRAV38-2DV8`> 0)


ncol(gdT.melanoma)
# 84537 cells left, 9533 filtered

table(gdT.melanoma$integrated_res.0.1)




# ###########################   Clustering    ###################################
################################################################################.
gc()
# Re-cluster at different resolutions to pick best one

gdT.melanoma <- FindVariableFeatures(gdT.melanoma,
                                     selection.method = "vst",
                                     nfeatures = 2000) #increased it if necessary



top10 <- head(VariableFeatures(gdT.melanoma))
p1 <- VariableFeaturePlot(gdT.melanoma)
p2 <- LabelPoints(plot = p1,
                  points = top10,
                  repel = T)

p1+p2



gdT.melanoma <- RunPCA(gdT.melanoma,
                       npcs = 60,
                       verbose = TRUE)



ElbowPlot(gdT.melanoma, ndims = 60)
ggsave(paste(dato,"ElbowPlot_60PCs.png",sep = "_"), 
       path = plots_path,
       width = 20, 
       height = 15, 
       units = "cm",
       dpi = 600, 
       scale=1.5, 
       bg = "white",
       limitsize = FALSE)


gdT.melanoma <- FindNeighbors(gdT.melanoma, dims = 1:40,
                              reduction = "harmony")

gdT.melanoma <- FindClusters(gdT.melanoma,
                             resolution = c(0.1, 0.2, 0.3, 0.4, 0.5,
                                            0.6, 0.7, 0.8, 0.9, 1),
                             cluster.name = c("integrated_res.0.1",
                                              "integrated_res.0.2",
                                              "integrated_res.0.3",
                                              "integrated_res.0.4",
                                              "integrated_res.0.5",
                                              "integrated_res.0.6",
                                              "integrated_res.0.7",
                                              "integrated_res.0.8",
                                              "integrated_res.0.9",
                                              "integrated_res.1"))


# Calculate UMAP again and name it as default
gdT.melanoma <- RunUMAP(gdT.melanoma,
                        dims = 1:40,
                        reduction = "harmony",
                        reduction.name = "umap")







#### ############################      QC     ###################################
################################################################################.

Idents(gdT.melanoma) <- gdT.melanoma@meta.data$integrated_res.0.1

DimPlot(gdT.melanoma, reduction = "umap", label = TRUE, repel = TRUE)

DimPlot(gdT.melanoma,
        group.by = "orig.ident",
        reduction = "umap",
        label = F)

ggsave(paste(dato,"post integration.png",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 15,
       height = 10,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)



DimPlot(gdT.melanoma,
        group.by = "outcome",
        reduction = "umap",
        label = F)

ggsave(paste(dato,"outcome.png",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 12,
       height = 10,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)



DimPlot(gdT.melanoma,
        group.by = "timepoint",
        reduction = "umap",
        label = F)

ggsave(paste(dato,"timepoint.png",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 12,
       height = 10,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)



DimPlot(gdT.melanoma, 
        reduction = "umap", 
        split.by = "sampleID",
        raster = F, # display more than 100000 cells
        ncol = 12) 

ggsave(paste(dato,"sampleID.png",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 30,
       height = 15,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)



DimPlot(gdT.melanoma, 
        reduction = "umap", 
        split.by ="orig.ident",
        ncol = 5,
        raster = F) 
ggsave(paste(dato,"lane.png",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 25,
       height = 15,
       units = "cm",
       dpi = 600,
       scale=2,
       limitsize = FALSE)

ggsave(paste(dato,"lane.svg",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 25,
       height = 15,
       units = "cm",
       dpi = 600,
       scale=2,
       limitsize = FALSE)


DimPlot(gdT.melanoma, 
        reduction = "umap", 
        split.by ="outcome.treatment",
        ncol = 2,
        raster = F) 

ggsave(paste(dato,"outcome.treatment.png",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 15,
       height = 15,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)




DimPlot(gdT.melanoma, 
        reduction = "umap", 
        split.by ="group") 

ggsave(paste(dato,"hashtag.png",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 30,
       height = 6,
       units = "cm",
       dpi = 600,
       scale=1.5,
       limitsize = FALSE)






##########################     TRDV subsets     ################################
################################################################################.
gc()
vd1 <- subset(gdT.melanoma, 
              TRDV1 > 0 & TRDV2 == 0 & TRDV3 == 0 &
                  `TRAV14DV4`== 0 & `TRAV29DV5`== 0 & 
                  `TRAV23DV6`== 0 & `TRAV36DV7`== 0 & 
                  `TRAV38-2DV8`== 0) 

vd2 <- subset(gdT.melanoma, 
              TRDV2 > 0 & TRDV1 == 0 & TRDV3 == 0 & 
                  `TRAV14DV4`== 0 & `TRAV29DV5`== 0 & 
                  `TRAV23DV6`== 0 & `TRAV36DV7`== 0 & 
                  `TRAV38-2DV8`== 0)

vd3 <- subset(gdT.melanoma, 
              TRDV3 > 0 & TRDV1 == 0 & TRDV2 == 0 &
                  `TRAV14DV4`== 0 & `TRAV29DV5`== 0 & 
                  `TRAV23DV6`== 0 & `TRAV36DV7`== 0 & 
                  `TRAV38-2DV8`== 0)



##### #######################       Vd1      ###################################
gc()
vd1 <- FindVariableFeatures(vd1,
                            selection.method = "vst",
                            nfeatures = 2000) #increased it if necessary



vd1 <- RunPCA(vd1,
              npcs = 60,
              verbose = TRUE)



ElbowPlot(vd1, ndims = 60)
ggsave(paste(dato,"vd1_ElbowPlot_60PCs.png",sep = "_"), 
       path = plots_path_vd1,
       width = 20, 
       height = 15, 
       units = "cm",
       dpi = 600, 
       scale=1.5, 
       bg = "white",
       limitsize = FALSE)


vd1 <- FindNeighbors(vd1, dims = 1:30,
                     reduction = "harmony")

vd1 <- FindClusters(vd1,
                    resolution = c(0.1, 0.2, 0.3, 0.4, 0.5,
                                   0.6, 0.7, 0.8, 0.9, 1),
                    cluster.name = c("integrated_res.0.1",
                                     "integrated_res.0.2",
                                     "integrated_res.0.3",
                                     "integrated_res.0.4",
                                     "integrated_res.0.5",
                                     "integrated_res.0.6",
                                     "integrated_res.0.7",
                                     "integrated_res.0.8",
                                     "integrated_res.0.9",
                                     "integrated_res.1"))


vd1 <- RunUMAP(vd1,
               dims = 1:30,
               reduction = "harmony",
               reduction.name = "umap")


Idents(vd1) <- vd1@meta.data$integrated_res.0.2
DimPlot(vd1, 
        reduction = "umap", 
        split.by ="orig.ident",
        ncol = 5,
        raster = F) 

ggsave(paste(dato,"vd1_lane.svg",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 25,
       height = 15,
       units = "cm",
       dpi = 600,
       scale=2,
       limitsize = FALSE)


##### #######################       Vd2      ###################################
gc()
vd2 <- FindVariableFeatures(vd2,
                            selection.method = "vst",
                            nfeatures = 2000) #increased it if necessary


vd2 <- RunPCA(vd2,
              npcs = 60,
              verbose = TRUE)



ElbowPlot(vd2, ndims = 60)
ggsave(paste(dato,"vd2_ElbowPlot_60PCs.png",sep = "_"), 
       path = plots_path_vd2,
       width = 20, 
       height = 15, 
       units = "cm",
       dpi = 600, 
       scale=1.5, 
       bg = "white",
       limitsize = FALSE)


vd2 <- FindNeighbors(vd2, dims = 1:30,
                     reduction = "harmony")

vd2 <- FindClusters(vd2,
                    resolution = c(0.1, 0.2, 0.3, 0.4, 0.5,
                                   0.6, 0.7, 0.8, 0.9, 1),
                    cluster.name = c("integrated_res.0.1",
                                     "integrated_res.0.2",
                                     "integrated_res.0.3",
                                     "integrated_res.0.4",
                                     "integrated_res.0.5",
                                     "integrated_res.0.6",
                                     "integrated_res.0.7",
                                     "integrated_res.0.8",
                                     "integrated_res.0.9",
                                     "integrated_res.1"))


vd2 <- RunUMAP(vd2,
               dims = 1:30,
               reduction = "harmony",
               reduction.name = "umap")


Idents(vd2) <- vd2@meta.data$integrated_res.0.2
DimPlot(vd2, 
        reduction = "umap", 
        split.by ="orig.ident",
        ncol = 5,
        raster = F) 

ggsave(paste(dato,"vd2_lane.svg",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 25,
       height = 15,
       units = "cm",
       dpi = 600,
       scale=2,
       limitsize = FALSE)




##### #######################       vd3      ###################################
gc()
vd3 <- FindVariableFeatures(vd3,
                            selection.method = "vst",
                            nfeatures = 2000) #increased it if necessary

vd3 <- RunPCA(vd3,
              npcs = 60,
              verbose = TRUE)



ElbowPlot(vd3, ndims = 60)
ggsave(paste(dato,"vd3_ElbowPlot_60PCs.png",sep = "_"), 
       path = plots_path_vd3,
       width = 20, 
       height = 15, 
       units = "cm",
       dpi = 600, 
       scale=1.5, 
       bg = "white",
       limitsize = FALSE)


vd3 <- FindNeighbors(vd3, dims = 1:30,
                     reduction = "harmony")

vd3 <- FindClusters(vd3,
                    resolution = c(0.1, 0.2, 0.3, 0.4, 0.5,
                                   0.6, 0.7, 0.8, 0.9, 1),
                    cluster.name = c("integrated_res.0.1",
                                     "integrated_res.0.2",
                                     "integrated_res.0.3",
                                     "integrated_res.0.4",
                                     "integrated_res.0.5",
                                     "integrated_res.0.6",
                                     "integrated_res.0.7",
                                     "integrated_res.0.8",
                                     "integrated_res.0.9",
                                     "integrated_res.1"))


# Calculate UMAP again and name it as default
vd3 <- RunUMAP(vd3,
               dims = 1:30,
               reduction = "harmony",
               reduction.name = "umap")




Idents(vd3) <- vd3@meta.data$integrated_res.0.2
DimPlot(vd3, 
        reduction = "umap", 
        split.by ="orig.ident",
        ncol = 5,
        raster = F) 

ggsave(paste(dato,"vd3_lane.svg",sep = "_"),
       path = paste(plots_path, "integration_check", sep = "/"),
       width = 25,
       height = 15,
       units = "cm",
       dpi = 600,
       scale=2,
       limitsize = FALSE)







################################################################################.
##########################      Save data     ##################################
################################################################################.
gc()

saveRDS(object = gdT.melanoma, 
        file = paste(data_path, "/", "10xRNA_H-melanoma-PBMC-gdT.rds",sep=""))

saveRDS(object = vd1, 
        file = paste(data_path, "/", "10xRNA_Vd1.rds",sep=""))

saveRDS(object = vd2, 
        file = paste(data_path, "/", "10xRNA_Vd2.rds",sep=""))

saveRDS(object = vd3, 
        file = paste(data_path, "/", "10xRNA_Vd3.rds",sep=""))