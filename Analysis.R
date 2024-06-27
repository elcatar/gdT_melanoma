################################################################################.
##### --------                Melanoma gdT analysis               --------- ####
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
library(Seurat)
library(tidyr)
library(tidyverse)
library(RColorBrewer)




###### variables used through script ######---------------  --------------------

#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-",""), 3, -1)

# colour string for imputation and overlays
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))


# VB lab folder
# path <- paste("U:/Sund/Public/T-Cell-Signalling-and-Development/VB Lab/10x_data/10xRNA_H-melanoma-PBMC-gdT",
#               sep = "/")

# External HD
# External HD
path <- "E:/Sequencing/Melanoma-project"

setwd(path)

data_path <- paste(path, "/data",sep = "")



plots_path_gen <- paste(path, "/data/plots/99.paper",sep = "")

plots_path1 <- paste(plots_path_gen, "/10x-data-1",sep = "")

plots_path_vd1 <- paste(plots_path_gen, "10x-data-2/vd1", sep = "/")
plots_path_vd2 <- paste(plots_path_gen, "10x-data-2/Vd2", sep = "/")
plots_path_vd3 <- paste(plots_path_gen, "10x-data-2/Vd3", sep = "/")



##### #######################   read data    ###################################
################################################################################.


gdT.melanoma <- readRDS(paste(data_path,"/10xRNA_H-melanoma-PBMC-gdT.rds",
                              sep = ""))


Idents(gdT.melanoma) <- gdT.melanoma$integrated_res.0.1

vd1 <- readRDS(paste(data_path,
                     "/10xRNA_vd1.rds",
                     sep = ""))


vd2 <- readRDS(paste(data_path,
                     "/10xRNA_Vd2.rds",
                     sep = ""))

vd3 <- readRDS(paste(data_path,
                     "/10xRNA_Vd3.rds",
                     sep = ""))


##### #######################   functions    ###################################
################################################################################.
### calculate frequencies e.g. clusters -----------------------------------------
SmplFreq <- function(object, group.by){
    
    group.by <- enquo(group.by)
    
    # get sample names
    grps <- object$sampleID %>% 
        unique %>%  
        as.character()
    
    
    for (i in 1:length(grps)) {
        
        if (i ==1){
            
            ## initialize tibble with all the groups we want (so we can correct for NA)
            smpls <- object@meta.data %>% 
                select(sampleID, 
                       hash.ID, 
                       outcome,
                       timepoint,
                       outcome.long,
                       outcome.lane,
                       treatment,   
                       outcome.timepoint,
                       outcome.treatment,
                       group,
                       sampleID,
                       patientID) %>% 
                unique() %>% 
                remove_rownames() 
            
            frm <- object@meta.data %>% 
                select(!!group.by) %>% 
                unique() %>% 
                remove_rownames() %>% 
                cross_join(smpls)
            
            
            ## calculate the first sample
            tmp <- object@meta.data %>%  
                filter(sampleID == grps[i]) %>% 
                group_by(!!group.by) %>% 
                summarize(count = n()) %>% 
                mutate(sampleID = grps[i])
            
            freq <- frm %>% 
                filter(sampleID == grps[i]) %>% 
                full_join(tmp) %>% 
                mutate(count = if_else(is.na(count), 0, count)) %>% 
                mutate(relative_freq = count/sum(count)) %>% 
                mutate(freq = round(relative_freq*100,
                                    digits = 2))
            
            
        } else {
            ## calculate the first sample
            tmp <- object@meta.data %>%  
                filter(sampleID == grps[i]) %>% 
                group_by(!!group.by) %>% 
                summarize(count = n()) %>% 
                mutate(sampleID = grps[i])
            
            tmp2 <- frm %>% 
                filter(sampleID == grps[i]) %>% 
                full_join(tmp) %>% 
                mutate(count = if_else(is.na(count), 0, count)) %>% 
                mutate(relative_freq = count/sum(count)) %>% 
                mutate(freq = round(relative_freq*100,
                                    digits = 2))
            
            freq <- freq %>% 
                full_join(tmp2)
            
        }
    }
    
    return(freq)
    
}



##### ######################       DEGs      ###################################
################################################################################.
################################################################################.


###### ########################     Genearl      ###############################

##### DEGs between groups ......................................................

gc()
Idents(gdT.melanoma) <- gdT.melanoma$group
DEG.grp.gen <- FindAllMarkers(gdT.melanoma,
                              only.pos = F,          
                              test.use = "wilcox",   
                              min.pct = 0.25,         
                              logfc.threshold = 0.5 
) %>%
    filter(p_val_adj < 0.05)

write.csv(DEG.grp.gen, 
          file = str_c(table_path, "/1.DEG-between-groups_gen.csv"))




##### DEGs between clusters ....................................................

gc()
Idents(gdT.melanoma) <- gdT.melanoma$integrated_res.0.1
DEG.clusters.gen <- FindAllMarkers(gdT.melanoma,
                                   only.pos = F,          
                                   test.use = "wilcox",   
                                   min.pct = 0.25,         
                                   logfc.threshold = 0.5) %>%
    filter(p_val_adj < 0.05)

write.csv(DEG.clusters.gen, 
          file = str_c(table_path, "/5.DEG-between-Clusters_gen.csv"))






###### ########################     Vd1      ###################################

##### DEGs between groups ......................................................

gc()
Idents(vd1) <- vd1$group
DEG.grp.vd1 <- FindAllMarkers(vd1,
                              only.pos = F,          
                              test.use = "wilcox",   
                              min.pct = 0.25,         
                              logfc.threshold = 0.5) %>%
    filter(p_val_adj < 0.05)

write.csv(DEG.grp.vd1, 
          file = str_c(table_path, "/2.DEG-between-groups_Vd1.csv"))




##### DEGs between clusters ......................................................

gc()
Idents(vd1) <- vd1$integrated_res.0.2
DEG.clusters.vd1 <- FindAllMarkers(vd1,
                                   only.pos = F,          
                                   test.use = "wilcox",   
                                   min.pct = 0.25,         
                                   logfc.threshold = 0.5) %>%
    filter(p_val_adj < 0.05)

write.csv(DEG.clusters.vd1, 
          file = str_c(table_path, "/5.DEG-between-Clusters_gen.csv"))




###### ########################     Vd2      ###################################

##### DEGs between groups ......................................................

gc()
Idents(vd2) <- vd2$group
DEG.grp.vd2 <- FindAllMarkers(vd2,
                              only.pos = F,          
                              test.use = "wilcox",   
                              min.pct = 0.25,         
                              logfc.threshold = 0.5) %>%
    filter(p_val_adj < 0.05)

write.csv(DEG.grp.vd2, 
          file = str_c(table_path, "/2.DEG-between-groups_vd2.csv"))




##### DEGs between clusters ......................................................

gc()
Idents(vd2) <- vd2$integrated_res.0.2
DEG.clusters.vd2 <- FindAllMarkers(vd2,
                                   only.pos = F,          
                                   test.use = "wilcox",   
                                   min.pct = 0.25,         
                                   logfc.threshold = 0.5) %>%
    filter(p_val_adj < 0.05)

write.csv(DEG.clusters.vd2, 
          file = str_c(table_path, "/5.DEG-between-Clusters_gen.csv"))






###### ########################     Vd3      ###################################

##### DEGs between groups ......................................................

gc()
Idents(vd3) <- vd3$group
DEG.grp.vd3 <- FindAllMarkers(vd3,
                              only.pos = F,          
                              test.use = "wilcox",   
                              min.pct = 0.25,         
                              logfc.threshold = 0.5) %>%
    filter(p_val_adj < 0.05)

write.csv(DEG.grp.vd3, 
          file = str_c(table_path, "/2.DEG-between-groups_vd3.csv"))




##### DEGs between clusters ......................................................

gc()
Idents(vd3) <- vd3$integrated_res.0.2
DEG.clusters.vd3 <- FindAllMarkers(vd3,
                                   only.pos = F,          
                                   test.use = "wilcox",   
                                   min.pct = 0.25,         
                                   logfc.threshold = 0.5) %>%
    filter(p_val_adj < 0.05)

write.csv(DEG.clusters.vd3, 
          file = str_c(table_path, "/5.DEG-between-Clusters_gen.csv"))







##### ###################     General plots     ################################
################################################################################.
################################################################################.
###### General UMAP ------------------------------------------------------------
Idents(gdT.melanoma) <- gdT.melanoma@meta.data$integrated_res.0.1


DimPlot(subset(gdT.melanoma, downsample=10000), group.by = "integrated_res.0.1", 
        label = TRUE, reduction = "umap")
ggsave(paste(dato,"UMAP.res.0.1.svg",sep = "_"),
       path = plots_path1,
       width = 8.5,
       height = 8,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)


###### TCR dotplots ------------------------------------------------------------
TCR_genes <- c(rownames(gdT.melanoma[["RNA"]]$data)[startsWith(rownames(gdT.melanoma[["RNA"]]$data),"TRDV")],
               rownames(gdT.melanoma[["RNA"]]$data)[startsWith(rownames(gdT.melanoma[["RNA"]]$data),"TRGV")])
DotPlot(gdT.melanoma, features = TCR_genes) + RotatedAxis() + coord_flip()
ggsave(paste(dato,"TCR-dotplot_all.svg",sep = "_"), 
       path = plots_path1,
       width = 10, 
       height = 8, 
       units = "cm",
       dpi = 300, 
       scale=1.5, 
       limitsize = FALSE)

TCR_genes2 <- TCR_genes[! TCR_genes %in% c("TRGV6", "TRGV5P", "TRGV1")]
DotPlot(gdT.melanoma, features = TCR_genes2) + RotatedAxis() + coord_flip()
ggsave(paste(dato,"TCR-dotplot.svg",sep = "_"), 
       path = plots_path1,
       width = 10, 
       height = 8, 
       units = "cm",
       dpi = 300, 
       scale=1.5, 
       limitsize = FALSE)




###### Feature plots -----------------------------------------------------------
Idents(gdT.melanoma) <- gdT.melanoma@meta.data$integrated_res.0.1

feature_genes <- c("TRDV1", "TRDV2", "TRDV3",
                   "TBX21", "SLAMF7", "CCL5",
                   "GATA3",
                   "RORC","IL23R", "CCR6",
                   "GZMA", "GZMB", "GZMK", "GZMH",
                   "PDCD1","CTLA4","TIGIT", 
                   "LAG3", "HAVCR2", "TOX", "ZBTB16")

ps <- FeaturePlot(subset(gdT.melanoma, downsample=10000),
                  features = feature_genes,
                  reduction = "umap")


for (i in 1:length(feature_genes)) {
    
    p <- ps[[i]]
    gene <- feature_genes[i]
    
    ggsave(paste(dato, "feature", gene, ".svg",sep = "_"),
           plot = p,
           path = str_c(plots_path1, "/FeaturePlots"),
           width = 8.5,
           height = 8,
           units = "cm",
           dpi = 300,
           scale=1.5,
           limitsize = FALSE)
    
}


###### Cluster abundancies -----------------------------------------------------
meta.data.gdT.melanoma <- gdT.melanoma@meta.data

## calculate freq using custum function
c.freq <- SmplFreq(gdT.melanoma,integrated_res.0.1)

## save data
# write.csv(c.freq, 
#           file = str_c(plots_path1, "/cluster-freq.csv"))

# filtering out R_Baseline_A because of low counts (less than 50 cells)
c.freq %>% 
    filter(!str_detect(sampleID, "R_Baseline_A")) %>% 
    ggplot(aes(x=integrated_res.0.1, 
               y= freq,
               fill=group)) + 
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    scale_fill_brewer(palette = "Paired") +
    labs(y = "% of cells in each individual patient",
         x = "Cluster")

ggsave(paste(dato,"Clusters-frequencies.svg",sep = "_"),
       path = plots_path1,
       width = 25,
       height = 10,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)




## save individual plots for each cluster
pts <- c.freq %>% select(integrated_res.0.1) %>% unique %>% pull 

for (i in pts) {
    
    c.freq %>% 
        filter(!str_detect(sampleID, "R_Baseline_A")) %>% 
        filter(integrated_res.0.1 == i) %>% 
        ggplot(aes(x = group, 
                   y = freq,
                   fill = group)) + 
        geom_boxplot(outlier.shape = NA) +
        geom_jitter()+
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, 
                                         hjust = 1),
              plot.title = element_text(hjust = 0.5),
              text = element_text(size = 10))+
        scale_fill_brewer(palette = "Paired") +
        labs(y = "% of cells in each individual patient",
             x = "",
             title = i)
    
    
    
    ggsave(paste(dato,  "cluster", i, "freq.svg",sep = "_"),
           path = plots_path1,
           width = 15,
           height = 8,
           units = "cm",
           dpi = 300,
           scale=1.5,
           limitsize = FALSE)
    
    
    
}




###### Heatmaps         --------------------------------------------------------
cluster_genes <- c("KLRC1",
                   "IL7R",
                   "GZMK",
                   "THEMIS",
                   "KLRB1",
                   "TNFAIP3",
                   "NFKBIA",
                   "RPS26",
                   "RORA",
                   "CD69",
                   
                   "TIGIT",    
                   "KIR2DL3",
                   "KLRC3",
                   "NCR1",
                   "KIR3DL1",
                   "KIR3DL2",
                   "TOX",
                   "CD244",
                   "EOMES",
                   "CD8A",
                   
                   "ARL15",
                   "ANK3",
                   "MED13L",
                   "RAD51B",
                   "KLF12",
                   "IMMP2L",
                   "RUNX1",
                   "NCOA2",
                   "HIVEP2",
                   "NCOA1",
                   
                   "CCR7",  
                   "LEF1",
                   "NELL2",
                   "LTB",
                   "TCF7",
                   "TNFSF8",
                   "BACH2",
                   "CD27",
                   "SELL",
                   "CD27",
                   "S1PR1")


Idents(gdT.melanoma) <- gdT.melanoma@meta.data$integrated_res.0.1
DoHeatmap(subset(gdT.melanoma, downsample=400),
          features = cluster_genes)+ 
    scale_fill_gradientn(colours = mycols)

ggsave(paste(dato,"DEG_clusters_Heatmap.pdf",sep = "_"),
       path = plots_path1,
       width = 15,
       height = 15,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)





ap1.genes <- c("JUN", "JUNB", "JUND",
               "FOS", "FOSB", "FOSL1", "FOSL2",
               "ATF1", "ATF2", "ATF3", "ATF4", "ATF5", "ATF6", "ATF6B", "ATF7", 
               "BATF", "BATF2", "BATF3", "JDP2",
               "MAF", "MAFA", "MAFB", "MAFF", "MAFG" , "MAFK",
               "NR4A1", "NR4A2", "NR4A3")


Idents(gdT.melanoma) <- gdT.melanoma$group
DoHeatmap(subset(gdT.melanoma, downsample=400),
          features = ap1.genes)+ 
    scale_fill_gradientn(colours = mycols)

ggsave(paste(dato,"AP-1_Heatmap.pdf",sep = "_"),
       path = plots_path1,
       width = 15,
       height = 10,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)










##### ####################       Vd1 plots      ################################
################################################################################.
################################################################################.

######  UMAP ------------------------------------------------------------
Idents(vd1) <- vd1@meta.data$integrated_res.0.2


DimPlot(vd1, group.by = "integrated_res.0.2", 
        label = TRUE, reduction = "umap")
ggsave(paste(dato,"UMAP.res.0.2.svg",sep = "_"),
       path = plots_path1,
       width = 8.5,
       height = 8,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)




######  Heatmaps     -------------------------------------------------------

vd1.genes <- read.csv(paste(data_path, "vd1_cluster_deg_selected.csv", sep = "/")) %>% 
    pull(gene)


Idents(vd1) <- vd1$integrated_res.0.2
DoHeatmap(subset(vd1, downsample=400),
          features = vd1.genes)+ 
    scale_fill_gradientn(colours = mycols)

ggsave(paste(dato,"Vd1_DEG_clusters_Heatmap.pdf",sep = "_"),
       path = plots_path_vd1,
       width = 15,
       height = 10,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)


#########   Vln Vd1 -------------------------------------------------------
Idents(vd1) <- vd1$integrated_res.0.2
vln_genes <- c( "PDCD1","CTLA4","TIGIT", 
                "LAG3", "HAVCR2", "TOX",
                "KIR3DX1",
                "KIR2DL3",
                "KIR2DL1",
                "KIR2DL4",
                "KIR3DL1",
                "KIR3DL2",
                "KIRREL3",
                "KIR3DL3",
                "KIRREL2"
)

clrs <-brewer.pal(n = 9, name = "Paired")

ps <- VlnPlot(vd1,
              features = vln_genes,
              pt.size = 0,
              cols = clrs)

for (i in 1:length(vln_genes)) {
    
    p <- ps[[i]]
    gene <- vln_genes[i]
    
    ggsave(paste(dato,"vd1",  "vln", gene, ".svg",sep = "_"),
           plot = p,
           path = str_c(plots_path_vd1, "/ICR_Plots"),
           width = 8.5,
           height = 8,
           units = "cm",
           dpi = 300,
           scale=1.5,
           limitsize = FALSE)
    
}







###### Log2FC plots ---------------------------------------------------------------
interesting_genes_vd1 <- read.csv(paste(data_path, "vd1_group_deg_selected.csv", sep = "/")) %>% 
    pull(gene)

selected_vd1 <- DEG.grp.vd1 %>% 
    filter(gene %in% interesting_genes_vd1)

selected_all_vd1 <- DEG.grp.vd1 %>% 
    filter(gene %in% interesting_genes_all)




#Create a copy of the data and a jittered version of the x variable
datJit_vd1  <- selected_vd1 
datJit_vd1$group <- jitter(as.numeric(factor(selected_vd1$cluster)), 1.8)


#Create a boxplot, overlay the jittered points and
# label the bottom 1% points
ggplot(selected_vd1,
       aes(x = cluster,
           y = avg_log2FC)) +
    geom_blank()+
    geom_point(data=datJit_vd1,
               aes(x=group)) +
    theme_bw() +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5))+
    geom_hline(yintercept = 0) +
    geom_text_repel(data = datJit_vd1,
                    aes(x = group,
                        label = gene),
                    size = 3,
                    max.overlaps = 50,
                    min.segment.length = 0) 


ggsave(paste(dato,"vd1_interesting_genes.svg",sep = "_"),
       path = plots_path_vd1,
       width = 20,
       height = 8,
       units = "cm",
       dpi = 600,
       # useDingbats=T,
       scale=1.5,
       limitsize = FALSE)





##### ####################       Vd2 plots      ################################
################################################################################.
################################################################################.

######  UMAP ------------------------------------------------------------
Idents(vd2) <- vd2@meta.data$integrated_res.0.2


DimPlot(vd2, group.by = "integrated_res.0.2", 
        label = TRUE, reduction = "umap")
ggsave(paste(dato,"UMAP.res.0.2.svg",sep = "_"),
       path = plots_path1,
       width = 8.5,
       height = 8,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)




######  Heatmaps     -------------------------------------------------------

vd2.genes <- read.csv(paste(data_path, "vd2_cluster_deg_selected.csv", sep = "/")) %>% 
    pull(gene)


Idents(vd2) <- vd2$integrated_res.0.2
DoHeatmap(subset(vd2, downsample=400),
          features = vd2.genes)+ 
    scale_fill_gradientn(colours = mycols)

ggsave(paste(dato,"vd2_DEG_clusters_Heatmap.pdf",sep = "_"),
       path = plots_path_vd2,
       width = 15,
       height = 10,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)



###### Log2FC plots ---------------------------------------------------------------
interesting_genes_vd2 <- read.csv(paste(data_path, "vd2_group_deg_selected.csv", sep = "/")) %>% 
    pull(gene)

selected_vd2 <- DEG.grp.vd2 %>% 
    filter(gene %in% interesting_genes_vd2)

selected_all_vd2 <- DEG.grp.vd2 %>% 
    filter(gene %in% interesting_genes_all)




#Create a copy of the data and a jittered version of the x variable
datJit_vd2  <- selected_vd2 
datJit_vd2$group <- jitter(as.numeric(factor(selected_vd2$cluster)), 1.8)


#Create a boxplot, overlay the jittered points and
# label the bottom 1% points
ggplot(selected_vd2,
       aes(x = cluster,
           y = avg_log2FC)) +
    geom_blank()+
    geom_point(data=datJit_vd2,
               aes(x=group)) +
    theme_bw() +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5))+
    geom_hline(yintercept = 0) +
    geom_text_repel(data = datJit_vd2,
                    aes(x = group,
                        label = gene),
                    size = 3,
                    max.overlaps = 50,
                    min.segment.length = 0) 


ggsave(paste(dato,"vd2_interesting_genes.svg",sep = "_"),
       path = plots_path_vd2,
       width = 20,
       height = 8,
       units = "cm",
       dpi = 600,
       # useDingbats=T,
       scale=1.5,
       limitsize = FALSE)





##### ####################       Vd3 plots      ################################
################################################################################.
################################################################################.

######  UMAP ------------------------------------------------------------
Idents(vd3) <- vd3@meta.data$integrated_res.0.2


DimPlot(vd3, group.by = "integrated_res.0.2", 
        label = TRUE, reduction = "umap")
ggsave(paste(dato,"UMAP.res.0.2.svg",sep = "_"),
       path = plots_path1,
       width = 8.5,
       height = 8,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)




######  Heatmaps     -------------------------------------------------------

vd3.genes <- read.csv(paste(data_path, "vd3_cluster_deg_selected.csv", sep = "/")) %>% 
    pull(gene)


Idents(vd3) <- vd3$integrated_res.0.2
DoHeatmap(subset(vd3, downsample=400),
          features = vd3.genes)+ 
    scale_fill_gradientn(colours = mycols)

ggsave(paste(dato,"vd3_DEG_clusters_Heatmap.pdf",sep = "_"),
       path = plots_path_vd3,
       width = 15,
       height = 10,
       units = "cm",
       dpi = 300,
       scale=1.5,
       limitsize = FALSE)



###### Log2FC plots ---------------------------------------------------------------
interesting_genes_vd3 <- read.csv(paste(data_path, "vd3_group_deg_selected.csv", sep = "/")) %>% 
    pull(gene)

selected_vd3 <- DEG.grp.vd3 %>% 
    filter(gene %in% interesting_genes_vd3)

selected_all_vd3 <- DEG.grp.vd3 %>% 
    filter(gene %in% interesting_genes_all)




#Create a copy of the data and a jittered version of the x variable
datJit_vd3  <- selected_vd3 
datJit_vd3$group <- jitter(as.numeric(factor(selected_vd3$cluster)), 1.8)


#Create a boxplot, overlay the jittered points and
# label the bottom 1% points
ggplot(selected_vd3,
       aes(x = cluster,
           y = avg_log2FC)) +
    geom_blank()+
    geom_point(data=datJit_vd3,
               aes(x=group)) +
    theme_bw() +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5))+
    geom_hline(yintercept = 0) +
    geom_text_repel(data = datJit_vd3,
                    aes(x = group,
                        label = gene),
                    size = 3,
                    max.overlaps = 50,
                    min.segment.length = 0) 


ggsave(paste(dato,"vd3_interesting_genes.svg",sep = "_"),
       path = plots_path_vd3,
       width = 20,
       height = 8,
       units = "cm",
       dpi = 600,
       # useDingbats=T,
       scale=1.5,
       limitsize = FALSE)

