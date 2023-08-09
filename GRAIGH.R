library(Matrix)
library(dplyr)
library(tidyr)
library(Seurat)
library(Signac)
library(rtracklayer)
library(rhdf5)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(aricode)
library(caret)

#######

if (!(dir.exists("Results"))){
    dir.create("Results")
}
if (!(dir.exists("Results/IMAGES"))){
    dir.create("Results/IMAGES")
}
if (!(dir.exists("Results/Tables"))){
    dir.create("Results/Tables")
}

#### specificity functions

specificity_raw <- function(true_values, predicted_values){
    
    # Calculate True Positives (TP), True Negatives (TN), False Positives (FP), and False Negatives (FN)
    TP <- sum(true_values == 1 & predicted_values == 1)
    TN <- sum(true_values == 0 & predicted_values == 0)
    FP <- sum(true_values == 0 & predicted_values == 1)
    FN <- sum(true_values == 1 & predicted_values == 0)
    
    # Calculate Specificity
    specificity <- TN / (TN + FP)
    return(specificity)
}

specificity_analysis <- function(SEU, cell_type_table, cell_type_list, gene_selected){
    
    expression_vector <- data.frame(activity = SEU_obj@assays$ACTIVITY@counts[gene_selected,])
    cell_type_table$Activity <- expression_vector$activity[match(rownames(cell_type_table), rownames(expression_vector))]
    cell_type_table$Activity <- ifelse(cell_type_table$Activity != 0, 1,0 )
    
    gene_GHs <- GH_score %>%  dplyr::filter(symbol == gene_selected, is_elite == 1) #%>% top_n(n = -10, wt = p_val_adj) select(GHid)
    gh_list <- unique(gh_markers_ct %>% dplyr::filter(gene %in% gene_GHs$GHid,avg_log2FC > 0) %>%  select(gene))$gene# %>% top_n(n = -10, wt = avg_log2FC)
    
    acc_spec <- 0
    for (gh in gh_list){
        expression_vector <- data.frame(activity = SEU_obj@assays$GH@counts[gh,])
        cell_type_table$Accessibility <- expression_vector$activity[match(rownames(cell_type_table), rownames(expression_vector))]
        cell_type_table$Accessibility <- ifelse(cell_type_table$Accessibility != 0, 1,0 )
        
        cell_type_table <- cell_type_table %>% dplyr::mutate(bin_ct = ifelse(cell_type_table$predicted.id %in% cell_type_list, 1,0 ))
        
        acc_spec <- c(acc_spec ,specificity_raw(cell_type_table$bin_ct, cell_type_table$Accessibility))
    }
    acc_spec <- acc_spec[-1]
    result <- mean(acc_spec)
    
    specificity_raw(cell_type_table$bin_ct, cell_type_table$Activity)
    result <- c(result, specificity_raw(cell_type_table$bin_ct, cell_type_table$Activity), result-specificity_raw(cell_type_table$bin_ct, cell_type_table$Activity))
    return(result)
}


######## Loading Data #######

# GH data
gff <- readGFF("TMPDATA/Genehancer/GeneHancer_v5.15.gff")
GH_element <- read.csv(file = "TMPDATA/Genehancer/GeneHancer_AnnotSV_elements_v5.15.txt", sep = "\t")
GH_score <- read.csv(file = "TMPDATA/Genehancer/GeneHancer_AnnotSV_gene_association_scores_v5.15.txt", sep = "\t")

#scATAC-seq data
matrix <- readMM("TMPDATA/Human_PBMC/filtered_peak_bc_matrix/matrix.mtx")
cells <- read.table("TMPDATA/Human_PBMC/filtered_peak_bc_matrix/barcodes.tsv")
features <- read.delim("TMPDATA/Human_PBMC/filtered_peak_bc_matrix/peaks.bed", header=FALSE) %>% unite(features, sep = "-")
features <- features[1:165376,]
matrix <- matrix[1:165376,]
row.names(matrix) <- features
colnames(matrix) <- cells$V1

chrom_assay <- CreateChromatinAssay(
    counts = matrix,
    sep = c("-", "-"), 
    fragments = "TMPDATA/Human_PBMC/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(chrom_assay) <- annotations
SEU_obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
rm(chrom_assay)


#### Creation of the connection matrix ####

# peak coordinates
f <- read.delim("TMPDATA/Human_PBMC/filtered_peak_bc_matrix/peaks.bed", header=FALSE)
f <- f[1:165376,]
Peak_GRange <- GRanges(seqnames = f$V1, ranges =IRanges(start= f$V2 , end= f$V3))
metadata(Peak_GRange) <- as.list(features)

#GH coordinates
GH_GRange <- GRanges(seqnames = gff$seqid, ranges =IRanges(start= gff$start , end= gff$end) )
metadata(GH_GRange) <- gff[,-c(1:8)]

#overlapping them
overlaps <- findOverlaps(GH_GRange, Peak_GRange)

#creating matrix
connection_matrix <- sparseMatrix( 
  i = queryHits(overlaps),
  j = subjectHits(overlaps),
  dims = c(length(GH_GRange), length(Peak_GRange)),
  x = 1)

row.names(connection_matrix) <- GH_GRange@metadata[["genehancer_id"]]
colnames(connection_matrix) <- Peak_GRange@metadata[["features"]]

connection_matrix <- connection_matrix[rowSums(connection_matrix) !=0,]

#peak to GH distribution
peaks_per_GH <- as_tibble(rowSums(connection_matrix))
peaks_per_GH %>%  mutate(category = "1")%>% ggplot( aes(x=value, fill = category)) + geom_bar() + ggtitle("Peaks per GH element") + NoLegend() + labs(x = "Number of Peaks", y = "Gh elements") +
    theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"), axis.text = element_text(size=20), axis.title = element_text(size=30)) +
    scale_y_log10() + scale_fill_manual(values = ("#440156"))+ theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray90"),
                                                                     panel.grid.minor = element_line(color = "gray90")) #+ theme_minimal()
ggsave(path = "Results/IMAGES/", filename = "peaks_per_GH.pdf", width = 1080, height = 1080, units= "px",scale = 3.5)


#Create GH matrix
GH_cells_matrix <- connection_matrix %*% matrix

SEU_obj[["GH"]] <- CreateAssayObject(GH_cells_matrix)


##############

#processing ATAC matrix
DefaultAssay(SEU_obj) <- "ATAC"

SEU_obj <- RunTFIDF(SEU_obj)
SEU_obj <- FindTopFeatures(SEU_obj, min.cutoff = 'q0')
SEU_obj <- RunSVD(SEU_obj, reduction.name = "lsi.atac")
SEU_obj <- RunUMAP(object = SEU_obj, reduction = 'lsi.atac', dims = 2:30, reduction.name = "umap.ATAC")
SEU_obj <- FindNeighbors(object = SEU_obj, reduction = 'lsi.atac', dims = 2:30)
SEU_obj <- FindClusters(object = SEU_obj, verbose = FALSE, algorithm = 3)
DimPlot(SEU_obj, label = TRUE, reduction = "umap.ATAC", pt.size =2, label.size = 20) + ggtitle("ATAC embedding") + NoLegend() + 
    theme( axis.title = element_text(size=30), plot.title = element_text(size = 40, hjust = 0.5, face = "bold"), axis.text = element_text(size = 20)) 

ggsave(path = "Results/IMAGES/", filename = "ATAC_embedding.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

#Processing GH matrix
DefaultAssay(SEU_obj) <- "GH"
SEU_obj <- RunTFIDF(SEU_obj)
SEU_obj <- FindTopFeatures(SEU_obj, min.cutoff = 'q0')
SEU_obj <- RunSVD(SEU_obj, reduction.name = "lsi.gh")
SEU_obj <- RunUMAP(object = SEU_obj, reduction = 'lsi.gh', dims = 2:30, reduction.name = "umap.GH")
SEU_obj <- FindNeighbors(object = SEU_obj, reduction = 'lsi.gh', dims = 2:30)
SEU_obj <- FindClusters(object = SEU_obj, verbose = FALSE, algorithm = 3)
DimPlot(SEU_obj, label = TRUE, reduction = "umap.GH", pt.size =2, label.size = 20) + ggtitle("GH embedding") + NoLegend() + labs(x = "UMAP_1", y = "UMAP_2") +
    theme( axis.title = element_text(size=30), plot.title = element_text(size = 40, hjust = 0.5, face = "bold"), axis.text = element_text(size = 20)) 
ggsave(path = "Results/IMAGES/", filename = "GH_embedding.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

#superposing clusters
Idents(SEU_obj) <- "ATAC_snn_res.0.8"
DimPlot(object = SEU_obj, label = TRUE, reduction = "umap.GH") + NoLegend()
Idents(SEU_obj) <- "GH_snn_res.0.8"
DimPlot(SEU_obj, label = TRUE, reduction = "umap.ATAC") + NoLegend()

ARI <- ARI(SEU_obj@meta.data[["ATAC_snn_res.0.8"]], SEU_obj@meta.data[["GH_snn_res.0.8"]])
AMI <- AMI(SEU_obj@meta.data[["ATAC_snn_res.0.8"]], SEU_obj@meta.data[["GH_snn_res.0.8"]])

print("Clustering metrics:")
print(paste("ARI:", round(ARI, 3)))
print(paste("AMI:", round(AMI, 3)))

# creating and processsing GAM
DefaultAssay(SEU_obj) <- "ATAC"
gene.activities <- GeneActivity(SEU_obj)

SEU_obj[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(SEU_obj) <- "ACTIVITY"
SEU_obj <- FindVariableFeatures(SEU_obj, nfeatures = 3000)
SEU_obj <- NormalizeData(
    object = SEU_obj,
    assay = 'ACTIVITY',
    normalization.method = 'LogNormalize'
)
SEU_obj <- ScaleData(SEU_obj)
SEU_obj <- RunPCA(SEU_obj, npcs = 30)
SEU_obj <- RunUMAP(SEU_obj, dims = 1:30, reduction.name = "umap.activity")
SEU_obj <- FindNeighbors(SEU_obj, dims = 1:30)
SEU_obj <- FindClusters(SEU_obj, resolution = 0.5, algorithm = 3)

#Seurat Label transfer integration
DefaultAssay(SEU_obj) <- "ACTIVITY"
pbmc_rna <- readRDS("TMPDATA/pbmc_10k_v3.rds")

transfer.anchors <- FindTransferAnchors(
    reference = pbmc_rna,
    query = SEU_obj,
    reduction = 'cca'
)

predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = pbmc_rna$celltype,
    weight.reduction = SEU_obj[['lsi.atac']],
    dims = 2:30
)

SEU_obj <- AddMetaData(object = SEU_obj, metadata = predicted.labels)
DimPlot(SEU_obj, label = TRUE,  reduction = "umap.ATAC", group.by = "predicted.id", pt.size =2, label.size = 13, repel = TRUE) + ggtitle("Cell Types") + NoLegend() + labs(x = "UMAP_1", y = "UMAP_2") +
    theme( axis.title = element_text(size=30), plot.title = element_text(size = 40, hjust = 0.5, face = "bold"), axis.text = element_text(size = 20)) 
ggsave(path = "Results/IMAGES/", filename = "Cell_types.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)


#cell_type_table <- read.csv("cell_type_table.csv", col.names = "predicted.id")
#SEU_obj <- AddMetaData(object = SEU_obj, metadata = cell_type_table)

rm(matrix)

##### Differential Analysis #####

DefaultAssay(SEU_obj) <- "GH"
gh_markers <- FindAllMarkers(SEU_obj, assay = "GH")
write.CSV(gh_markers, file = "Results/gh_markers.csv")
#gh_markers <- read.csv(file = "Results/gh_markers.csv", row.names = 1)

gh_markers_ct <- FindAllMarkers(SEU_obj, assay = "GH", group.by = "predicted.id")
write.CSV(gh_markers_ct, file = "gh_markers_ct.csv")
#gh_markers_ct <- read.csv(file = "gh_markers_ct.csv", row.names = 1)


## CD4 +
CCR7_GHs <- GH_score %>%  dplyr::filter(symbol == "CCR7", is_elite == 1) #%>% top_n(n = -10, wt = p_val_adj) select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% CCR7_GHs$GHid,avg_log2FC > 0)# %>% top_n(n = -10, wt = avg_log2FC)

IL7R_GHs <- GH_score %>%  dplyr::filter(symbol == "IL7R",is_elite == 1) %>% dplyr::select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% IL7R_GHs$GHid,avg_log2FC > 0)# %>% top_n(n = -10, wt = avg_log2FC)

CD4_GHs <- GH_score %>%  dplyr::filter(symbol == "CD4",is_elite == 1) %>% dplyr::select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% CD4_GHs$GHid,avg_log2FC > 0)# %>% top_n(n = -10, wt = avg_log2FC)

# CD14 + 
CD14_GHs <- GH_score %>%  dplyr::filter(symbol == "CD14",is_elite == 1) %>% dplyr::select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% CD14_GHs$GHid,avg_log2FC > 0) %>% top_n(n = -10, wt = p_val_adj)

# CD8 +
CD8A_GHs <- GH_score %>%  dplyr::filter(symbol == "CD8A",is_elite == 1) %>% dplyr::select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% CD8A_GHs$GHid,avg_log2FC > 0) #%>% top_n(n = -10, wt = p_val_adj)

#B cells
MS4A1_GHs <- GH_score %>%  dplyr::filter(symbol == "MS4A1",is_elite == 1) %>% dplyr::select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% MS4A1_GHs$GHid, avg_log2FC > 0) #%>% top_n(n = -10, wt = avg_log2FC)

#NK cells
GNLY_GHs <- GH_score %>%  dplyr::filter(symbol == "GNLY", is_elite == 1) #%>% top_n(n = -10, wt = p_val_adj) select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% GNLY_GHs$GHid,avg_log2FC > 0)# %>% top_n(n = -10, wt = avg_log2FC)

# dendritic
FCER1A_GHs <- GH_score %>%  dplyr::filter(symbol == "FCER1A", is_elite == 1) #%>% top_n(n = -10, wt = p_val_adj) select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% FCER1A_GHs$GHid,avg_log2FC > 0)# %>% top_n(n = -10, wt = avg_log2FC)

CST3_GHs <- GH_score %>%  dplyr::filter(symbol == "CST3", is_elite == 1) #%>% top_n(n = -10, wt = p_val_adj) select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% CST3_GHs$GHid,avg_log2FC > 0)# %>% top_n(n = -10, wt = avg_log2FC)

# plat
PPBP_GHs <- GH_score %>%  dplyr::filter(symbol == "PPBP", is_elite == 1) #%>% top_n(n = -10, wt = p_val_adj) select(GHid)
gh_markers_ct %>% dplyr::filter(gene %in% PPBP_GHs$GHid,avg_log2FC > 0)# %>% top_n(n = -10, wt = avg_log2FC)


#saveRDS(SEU_obj, file = "SEU_obj.rds")
#SEU_obj <- readRDS("SEU_obj.rds")

#### specificty

cell_type_table <- read.csv("cell_type_table.csv", col.names = "predicted.id")

cell_types_names <- c("CD4+ T cells", "CD4+ T cells", "CD4+ T cells", "CD8+ T cells", "Monocytes","Monocytes", "NK cells", "Dendritic cells", " Dendritic cells", "Bcells")
gene_marker_names <- c("CCR7", "IL7R", "CD4","CD8A", "CD14", "MS4A7","GNLY", "FCER1A", "CST3", "MS4A1")

specificty_table <- data.frame(Cell_Type = cell_types_names, Marker = gene_marker_names)


## CD4 +
cell_type_list <-  c( "CD4 Naive")
CCR7_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "CCR7")
IL7R_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "IL7R")
CD4_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "CD4")
 

# CD8 +
cell_type_list <-  c( "CD8 Naive", "CD8 effector")
CD8A_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "CD8A")

# CD14 + 
cell_type_list <-  c( "CD14+ Monocytes","CD16+ Monocytes")
CD14_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "CD14")

cell_type_list <-  c("CD16+ Monocytes")
MS4A7_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "MS4A7")

#NK cells
cell_type_list <-  c( "NK dim", "NK bright")
GNLY_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "GNLY")

# dendritic
cell_type_list <-  c( "Dendritic cell")
FCER1A_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "FCER1A")
CST3_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "CST3")

#B cells
cell_type_list <-  c( "B cell progenitor", "pre-B cell")
MS4A1_spec <- specificity_analysis(SEU_obj, cell_type_table, cell_type_list, "MS4A1")

numbers <- rbind(CCR7_spec, IL7R_spec, CD4_spec, CD8A_spec, CD14_spec, MS4A7_spec, GNLY_spec, FCER1A_spec, CST3_spec, MS4A1_spec)
numbers <- as.data.frame(numbers)
colnames(numbers) <- c("Accessibility", "Activity", "Delta")

specificty_table <- cbind(specificty_table, numbers)


FeaturePlot(SEU_obj, features = c("CD14"), pt.size = 2.5 , reduction = "umap.ATAC") + ggtitle("CD14 activity") + theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"), axis.text = element_text(size=20), axis.title = element_text(size=30)) +
    scale_color_gradient(low = "grey85",high = "#440156")
ggsave(path = "Results/IMAGES/", filename = "CD14_activity.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

FeaturePlot(SEU_obj, features = c("GH05J140611"), pt.size = 2.5 , reduction = "umap.ATAC") + ggtitle("GH05J140611 accessibility") + 
    theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"), axis.text = element_text(size=20), axis.title = element_text(size=30)) +
    scale_color_gradient(low = "grey85",high = "#440156")
ggsave(path = "Results/IMAGES/", filename = "CD14_accessiblity.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)



FeaturePlot(SEU_obj, features = c("CD4"), pt.size = 2.5 , reduction = "umap.ATAC") + ggtitle("CD4 activity") + theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"), axis.text = element_text(size=20), axis.title = element_text(size=30)) +
    scale_color_gradient(low = "grey85",high = "#440156")
ggsave(path = "Results/IMAGES/", filename = "CD_activity.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

FeaturePlot(SEU_obj, features = c("GH12J006784"), pt.size = 2.5 , reduction = "umap.ATAC") + ggtitle("GH12J006784 accessibility") + 
    theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"), axis.text = element_text(size=20), axis.title = element_text(size=30)) +
    scale_color_gradient(low = "grey85",high = "#440156")
ggsave(path = "Results/IMAGES/", filename = "CD4_accessiblity.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)





