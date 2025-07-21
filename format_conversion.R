library(Seurat)
library(SeuratData)
#library(SeuratDisk)
#Convert("./h5_collection/he2020integrating_staahl2016visualization_ST_Luminal_B_BC23270_E2.h5ad", 
#        dest = "h5seurat", overwrite = FALSE)
#ts <- LoadH5Seurat("./h5_collection/he2020integrating_staahl2016visualization_ST_Luminal_B_BC23270_E2.h5seurat")
# ts
library(zellkonverter)
library(org.Hs.eg.db)

get.seurat <- function(path.name, file.name){
  sce1=readH5AD(path.name, verbose = TRUE)
  adata_Seurat <- as.Seurat(sce1, counts = "X", data = NULL)
  adata_Seurat <- RenameAssays(object = adata_Seurat, originalexp = "RNA")
  meta <- adata_Seurat@meta.data
  spatial <- adata_Seurat@reductions$spatial@cell.embeddings
  spatial <- as.data.frame(spatial)
  spatial[, 1] <- as.numeric(spatial[, 1])
  spatial[, 2] <- as.numeric(spatial[, 2])
  mtx <- adata_Seurat@assays$RNA@counts
  gvec <- as.vector(rownames(mtx))
  print(gvec)
  annots <- select(org.Hs.eg.db, keys=gvec, 
                   columns="SYMBOL", keytype="ENSEMBL")
  annots <- annots[complete.cases(annots), ]
  annots <- annots[!duplicated(annots$ENSEMBL), ]
  annots <- annots[!duplicated(annots$SYMBOL), ]
  
  # table(duplicated(annots$SYMBOL))
  mtx <- mtx[annots$ENSEMBL, ]
  rownames(mtx) <- annots$SYMBOL
  sr <- CreateSeuratObject(counts = mtx, meta.data = meta)
  sr[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(spatial), key = "spatial_")
  saveRDS(sr, file = paste0("./h5_collection/", file.name, ".rds"))
}

get.seurat.burglund <- function(path.name, file.name){
  sce1=readH5AD(path.name, verbose = TRUE)
  adata_Seurat <- as.Seurat(sce1, counts = "X", data = NULL)
  adata_Seurat <- RenameAssays(object = adata_Seurat, originalexp = "RNA")
  meta <- adata_Seurat@meta.data
  spatial <- adata_Seurat@reductions$spatial@cell.embeddings
  spatial <- as.data.frame(spatial)
  spatial[, 1] <- as.numeric(spatial[, 1])
  spatial[, 2] <- as.numeric(spatial[, 2])
  mtx <- adata_Seurat@assays$RNA@counts
  gvec <- as.vector(rownames(mtx))
#  print(gvec)
#  annots <- select(org.Hs.eg.db, keys=gvec, 
#                   columns="SYMBOL", keytype="ENSEMBL")
#  annots <- annots[complete.cases(annots), ]
#  annots <- annots[!duplicated(annots$ENSEMBL), ]
#  annots <- annots[!duplicated(annots$SYMBOL), ]
  
  # table(duplicated(annots$SYMBOL))
#  mtx <- mtx[annots$ENSEMBL, ]
  gvec1 <- gsub(" .*", "", gvec)
  rownames(mtx) <- gvec1
  print(mtx[1:5, 1:5])
  sr <- CreateSeuratObject(counts = mtx, meta.data = meta)
  sr[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(spatial), key = "spatial_")
  saveRDS(sr, file = paste0("./h5_collection/", file.name, ".rds"))
}

for (item in list.files("./h5_collection/", pattern = "^he.*h5ad")){
  print(item)
  get.seurat(path.name = paste0("./h5_collection/", item), 
             file.name = gsub("\\..*", "", item))
}
for (item in list.files("./h5_collection/", pattern = "^berg.*h5ad")){
  print(item)
  get.seurat.burglund(path.name = paste0("./h5_collection/", item), 
             file.name = gsub("\\..*", "", item))
}
