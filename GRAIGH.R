BiocManager::install("EnsDb.Hsapiens.v86")

library(Matrix)
library(dplyr)
library(tidyr)
library(Seurat)
library(Signac)
library(rtracklayer)
library(rhdf5)
library(ggplot2)
library(EnsDb.Hsapiens.v86)


gff <- readGFF("TMPDATA/Genehancer/GeneHancer_v5.15.gff")

#GH_element <- read.csv(file = "TMPDATA/Genehancer/GeneHancer_AnnotSV_hg19_v5.15.txt", sep = "\t")
GH_element <- read.csv(file = "TMPDATA/Genehancer/GeneHancer_AnnotSV_elements_v5.15.txt", sep = "\t")
GH_score <- read.csv(file = "TMPDATA/Genehancer/GeneHancer_AnnotSV_gene_association_scores_v5.15.txt", sep = "\t")

matrix <- readMM("TMPDATA/Human_PBMC/filtered_peak_bc_matrix/matrix.mtx")
cells <- read.table("TMPDATA/Human_PBMC/filtered_peak_bc_matrix/barcodes.tsv")
features <- read.delim("TMPDATA/Human_PBMC/filtered_peak_bc_matrix/peaks.bed", header=FALSE) %>% unite(features, sep = "-")
motif_mapping <- read.delim("TMPDATA/Human_PBMC/10k_pbmc_ATACv2_nextgem_Chromium_Controller_peak_motif_mapping.bed", header=FALSE) %>% unite(V1,V2,V3, col = "peaks", sep = "-")

row.names(matrix) <- features$features
colnames(matrix) <- cells$V1

chrom_assay <- CreateChromatinAssay(
    counts = matrix,
    sep = c("-", "-"), #indicates which are the separators for the peak notation of the file
    fragments = "TMPDATA/Human_PBMC/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(chrom_assay) <- annotations
SEU_obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

f <- read.delim("TMPDATA/Human_PBMC/filtered_peak_bc_matrix/peaks.bed", header=FALSE)
Peak_GRange <- GRanges(seqnames = f$V1, ranges =IRanges(start= f$V2 , end= f$V3))
metadata(Peak_GRange) <- features

GH_GRange <- GRanges(seqnames = gff$seqid, ranges =IRanges(start= gff$start , end= gff$end) )
metadata(GH_GRange) <- gff[,-c(1:8)]

overlaps <- findOverlaps(GH_GRange, Peak_GRange)

connection_matrix <- sparseMatrix( 
  i = queryHits(overlaps),
  j = subjectHits(overlaps),
  dims = c(length(GH_GRange), length(Peak_GRange)),
  x = 1)

row.names(connection_matrix) <- GH_GRange@metadata[["genehancer_id"]]
colnames(connection_matrix) <- Peak_GRange@metadata[["features"]]
            
            
connection_matrix <- connection_matrix[rowSums(connection_matrix) !=0,]

peaks_per_GH <- as_tibble(rowSums(connection_matrix)) #%>% group_by(value) %>% count()
#median(peaks_per_GH)
#mean(peaks_per_GH)
ggplot(peaks_per_GH, aes(x=value)) + geom_bar() + scale_y_log10()
GH_cells_matrix <- connection_matrix %*% matrix

SEU_obj[["GH"]] <- CreateAssayObject(GH_cells_matrix)

#######
unique_rows <- unique(motif_mapping$V4)
unique_cols <- unique(motif_mapping$peaks)
row_indices <- match(motif_mapping$V4, unique_rows)
col_indices <- match(motif_mapping$peaks, unique_cols)
motif_mapping_matrix <- sparseMatrix(
    i = row_indices,
    j = col_indices,
    x = 1,
    dims = c(length(unique_rows), length(unique_cols)),
    dimnames = list(unique_rows, unique_cols)
)
motif_cells_matrix <- motif_mapping_matrix %*% matrix[unique_cols,]

SEU_obj[["Motif"]] <- CreateAssayObject(motif_cells_matrix)

#####
DefaultAssay(SEU_obj) <- "ATAC"

SEU_obj <- RunTFIDF(SEU_obj)
SEU_obj <- FindTopFeatures(SEU_obj, min.cutoff = 'q0')
SEU_obj <- RunSVD(SEU_obj, reduction.name = "lsi.atac")
SEU_obj <- RunUMAP(object = SEU_obj, reduction = 'lsi.atac', dims = 2:30, reduction.name = "umap.ATAC")
SEU_obj <- FindNeighbors(object = SEU_obj, reduction = 'lsi.atac', dims = 2:30)
SEU_obj <- FindClusters(object = SEU_obj, verbose = FALSE, algorithm = 3)
DimPlot(SEU_obj, label = TRUE, reduction = "umap.ATAC") + NoLegend()

SEU_obj <- RunTFIDF(SEU_obj)
SEU_obj <- FindTopFeatures(SEU_obj, min.cutoff = 'q0')
SEU_obj <- RunSVD(SEU_obj, reduction.name = "lsi.gh")
SEU_obj <- RunUMAP(object = SEU_obj, reduction = 'lsi.gh', dims = 2:30, reduction.name = "umap.GH")
SEU_obj <- FindNeighbors(object = SEU_obj, reduction = 'lsi.gh', dims = 2:30)
SEU_obj <- FindClusters(object = SEU_obj, verbose = FALSE, algorithm = 3)
DimPlot(object = SEU_obj, label = TRUE, reduction = "umap.GH") + NoLegend()

Idents(SEU_obj) <- "ATAC_snn_res.0.8"
DimPlot(object = SEU_obj, label = TRUE, reduction = "umap.GH") + NoLegend()
Idents(SEU_obj) <- "GH_snn_res.0.8"
DimPlot(SEU_obj, label = TRUE, reduction = "umap.ATAC") + NoLegend()

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
DimPlot(SEU_obj)



DefaultAssay(SEU_obj) <- "GH"
FeaturePlot(SEU_obj, features = "GH22J046753")



DefaultAssay(SEU_obj) <- "GH"
gh_markers <- FindAllMarkers(SEU_obj, assay = "GH")
write.table(gh_markers, file = "gh_markers.csv")
gh_markers <- read.table(file = "gh")
    
top_gh <- gh_markers %>%  group_by(cluster) %>% top_n(n = 10)
GH_score %>%  filter(GHid == "GH22J046753")


CD8A_GHs <- GH_score %>%  filter(symbol == "CD8A",is_elite == 1) %>% select(GHid)
gh_markers %>% filter(gene %in% CD8A_GHs$GHid) # %>% top_n(n = -10, wt = p_val_adj)

MS4A1_GHs <- GH_score %>%  filter(symbol == "MS4A1",is_elite == 1) %>% select(GHid)
gh_markers %>% filter(gene %in% MS4A1_GHs$GHid) # %>% top_n(n = -10, wt = avg_log2FC)

CD14_GHs <- GH_score %>%  filter(symbol == "CD14",is_elite == 1) %>% select(GHid)
gh_markers %>% filter(gene %in% CD14_GHs$GHid) %>% top_n(n = -10, wt = p_val_adj)

IL7R_GHs <- GH_score %>%  filter(symbol == "IL7R",is_elite == 1) %>% select(GHid)
gh_markers %>% filter(gene %in% IL7R_GHs$GHid)# %>% top_n(n = -10, wt = avg_log2FC)

NKG7_GHs <- GH_score %>%  filter(symbol == "NKG7", is_elite == 1) #%>% top_n(n = -10, wt = p_val_adj) select(GHid)
gh_markers %>% filter(gene %in% NKG7_GHs$GHid)# %>% top_n(n = -10, wt = avg_log2FC)


#######

GH_score %>% dplyr::filter(symbol == "GTPBP6", is_elite == 1)

#########

pbmc_rna <- readRDS("../vignette_data/pbmc_10k_v3.rds")

transfer.anchors <- FindTransferAnchors(
    reference = pbmc_rna,
    query = pbmc,
    reduction = 'cca'
)

predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = pbmc_rna$celltype,
    weight.reduction = pbmc[['lsi']],
    dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)


