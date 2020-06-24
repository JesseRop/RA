#!/usr/bin/env Rscript
library(dplyr)
library(httr)
library(MAST)
library(unixtools)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(ggmin)
library(pheatmap)
library(ggrepel)
library(shiny)
library(shinythemes)
library(plotly)
library(magrittr)
library(rlang)
library(tidyr)
library(tibble)
library(tidyverse)
library(grid)
library(data.table)
library(shinycssloaders)
library(DT)


# setwd("/home/jr345y/ra_devt/")

##Reading in the list of precomputed Seurat objects (Resolution 0.15 0.25&0.35, 0.45&0.55 computed separately read in and combined then written to disk and )
alveri = readRDS("alivernini.8k.rds")
alveri = list(alveri)
# ra_macrophage_combined_umap_list_res<-readRDS("ra_macrophage_combined_umap_list_res.rds")

##Reading in the list of precomputed table of cluster markers
ra_macrophage_combined_clusters_tables_res<-readRDS("ra_macrophage_combined_clusters_tables_res.rds")

##Reading in the list of precomputed table of differential expressed (DE) genes
ra_macrophage_combined_de_tables = readRDS("ra_macrophage_combined_de_tables.rds")

##Reading in the list of precomputed table of cluster markers
ra_macrophage_combined_de_ggplots_table = readRDS("ra_macrophage_combined_de_ggplots_table.rds")

##Reading in all the genes that are present in all raw objects
all_genes_ra = rownames(alveri[[1]]@assays$RNA)

##Declaring and assigning variables
# dim=15

res1 = 1
res2 = 1
diff_res = 0
cluster_names = c("TREM2low", "TREM2high", "FOLR2+ID2+", "FOLR2highLYVE1+", "HLAhighCLEC10A+", "CD48highS100A12+", "CD48+SPP1+", "HLAhighISG15+", "FOLR2+ICAM1+")
fav_genes = c("CD48", "S100A9", "SPP1", "CLEC10A")
conditions = c("Healthy",  "UPA", "Naive RA", "Resistant RA", "Remission RA")
grps = unique(alveri[[1]]@meta.data$group)
pairwise <- combn(grps, 2)
conds = lapply(1:ncol(pairwise), function(x) paste(pairwise[,x], collapse = " VS ")) %>% unlist()
cluster.colours <-c("#6D67E8","dodgerblue2","deepskyblue","#6FEAFD","#95E949","#FF791A","#FF1E1A","#A42537", "#B625C8")
group.cols <- c("#66CC00", "orange","darkorange1","Red","#3399FF")
names(group.cols) <- levels(alveri[[1]]$group)
choice_gene = "TREM2"

if (!(all(exists("ra_macrophage_combined_clusters_tables_res"), exists("ra_macrophage_combined_de_tables"), exists("ra_macrophage_combined_de_ggplots_table")))) {
  
  ##Precomputing and saving the list of tables of cluster markers
  ra_macrophage_combined_clusters_tables_res = lapply(alveri, function(x) { 
    # clst = unique(Idents(x))
    DefaultAssay(x) = "RNA"
    lapply(0:(length(unique(x$FinalClusters))-1), function(y) {
      tryCatch(FindConservedMarkers(x, ident.1 = y, grouping.var = "group"), error=function(e) NULL)
      
    })
  })
  
  saveRDS(ra_macrophage_combined_clusters_tables_res, "ra_macrophage_combined_clusters_tables_res.rds")
  
  ##Generating pairwise list for all DE group comparisons per cluster
  grps = unique(alveri[[1]]@meta.data$group)
  pairwise <- combn(ra_grps, 2)
  
  ##Precomputing and saving the list of tables of DE genes per cluster
  ra_macrophage_combined_de_tables = lapply(alveri, function(x) { 
    # clst = unique(Idents(x))
    DefaultAssay(x) = "RNA"
    x$celltype.group <- paste(Idents(x), x$group, sep = "_")
    x$celltype <- Idents(x)
    Idents(x) <- "celltype.group"
    
    lapply(0:(length(unique(x$FinalClusters))-1), function(y) {
      lapply(1:ncol(pairwise), function(z) {
        tryCatch(FindMarkers(x, ident.1 = paste(y, pairwise[1,z], sep = "_"), ident.2 = paste(y, pairwise[2,z], sep = "_"), verbose = T, min.cells.group = 3), error=function(e) NULL)
        
      })
    })
  })
  
  saveRDS(ra_macrophage_combined_de_tables, "ra_macrophage_combined_de_tables.rds")
  
  ##Precomputing and saving the list of ggplot per cluster for all resolutions
  ra_macrophage_combined_de_ggplots_table = lapply(alveri, function(x) { 
    DefaultAssay(x) = "RNA"
    
    lapply(0:(length(unique(x$FinalClusters))-1), function(y) {
      # clst = unique(Idents(x))
      cells_type <- subset(x, idents = y)
      Idents(cells_type) <- "group"
      avg.cells <- log1p(AverageExpression(cells_type, verbose = FALSE)$RNA)
      avg.cells$gene <- rownames(avg.cells)
      avg.cells <- avg.cells %>% filter(!grepl("^mt-", gene))
      
    })
  })
  
  saveRDS(ra_macrophage_combined_de_ggplots_table, "ra_macrophage_combined_de_ggplots_table.rds")
  
}


##functions


which_numeric_cols = function(dat) {
  which(sapply(seq(ncol(dat)), function(i) {
    is.numeric(dat[,i])
  }))
}


