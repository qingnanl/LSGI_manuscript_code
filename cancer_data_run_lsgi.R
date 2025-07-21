library(dior)
library(Seurat)
library(tidyverse)
library(singlet)
library(RcppML)
library(Matrix)
library(anticlust)
library(copykat)
library(LSGI)
# source("./local.spatial.traj.R")

clean_glist <- function(glist){
  idx <- grep("^RPS|XIST|RP11|RPL|^MT|\\.", glist)
  glist <- glist[-idx]
  return(glist)
}
sr.process <- function(object){
  counts <- GetAssayData(object)
  counts <- counts[clean_glist(rownames(counts)),]
  object <- subset(object, features = rownames(counts))

  object <- object %>%
    NormalizeData() %>% 
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors() %>%
    FindClusters() 
  return(object)
}

#obj <- readRDS("./h5_collection/berglund2018spatial_berglund2018spatial_ST_mel2_rep1.rds")
#obj <- sr.process(obj)


find.ref <- function(object, gene.set, top.n = 0.05){
  object <- AddModuleScore(object, features = list(gene.set), name = "ref_")
  n.cells <- floor(ncol(object) * top.n)
  o <- order(object$ref_1, decreasing = T)[1:n.cells]
  cells <- colnames(object)[o]
  name.cluster <- names(which.max(prop.table(table(colnames(object) %in% cells, 
                                   object$seurat_clusters), margin = 2)[2, ]))
  ref.cells <- colnames(object)[object$seurat_clusters == name.cluster]
  return(ref.cells)
}


runcopykat <- function(object, sample.name, nc){
  exp.rawdata <- as.matrix(object@assays$RNA@counts)
  copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", 
                          ngene.chr=5, win.size=25, KS.cut=0.1, sam.name= sample.name, 
                          distance="euclidean", norm.cell.names=nc,
                          output.seg="FLASE", plot.genes="FALSE", genome="hg20",n.cores=24)
  pred.test <- data.frame(copykat.test$prediction)
  # pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
  CNA.test <- data.frame(copykat.test$CNAmat)
  return(list(CNA.test, pred.test))
}

gather.spatial.tumor <- function(obj, spatial_coords){
  spatial_coords$tumor_prediction <- obj$tumor
  return(spatial_coords)
}


imm.markers <- c("CD3E", "CD8A", "GZMK", "CD4", 
                 "CCR7", "GZMB", "FCER1G", "LHDB", 
                 "DUSP2", "IL7R", "S100A4")

# this function and code below (line 75 - line 145) is tailored for my file formats. 
# I will show how to run a single dataset below
run.panc.spatial.grad <- function(file.name, sample.name, gene.set, n.genes.factor = 50){
  if (length(grep("rds", file.name)) == 1){
          obj <- readRDS(file.name)
  }else{
          obj <- dior::read_h5(file = file.name, target.object = 'seurat')
  }
  #      obj <- dior::read_h5(file = file.name, target.object = 'seurat')
  # preprocess
  obj <- sr.process(obj)
  # copykat
  ref.cells <- find.ref(object = obj, gene.set = gene.set, top.n = 0.1)
  obj.cna <- runcopykat(object = obj, nc = ref.cells, sample.name = sample.name)
  obj.p <- read.delim(paste0(sample.name, "_copykat_prediction.txt"))
  rownames(obj.p) <- obj.p$cell.names
  obj$prediction <- obj.p[colnames(obj), "copykat.pred"]
  obj$tumor <- ifelse(obj$prediction == 'aneuploid', "tumor", "normal")
  # NMF (output NMF gene loading (top 50 and all))
  
  obj <- singlet::RunNMF(obj, k = c(6, 7, 8, 9, 10))
  saveRDS(obj, paste0("./output2/", sample.name, ".processed.obj.rds"))

  feature.loadings <- as.data.frame(obj@reductions$nmf@feature.loadings)

  top.gene.list <- list()
  for (i in 1:ncol(feature.loadings)){
    o <- order(feature.loadings[, i], decreasing = T)[1:n.genes.factor]
    features <- rownames(feature.loadings)[o]
    top.gene.list[[colnames(feature.loadings)[i]]] <- features
  }
  nmf.info <- list(feature.loadings = feature.loadings, top.genes = top.gene.list)


  # spatial grad
  nmfs <- obj@reductions$nmf@cell.embeddings
  if (!is.null(obj@reductions$spatial)){
          spatial_coords <- as.data.frame(obj@reductions$spatial@cell.embeddings)
  }else{
          spatial_coords <- as.data.frame(obj@images$slice1@coordinates[, c(4, 5)])
  }
  colnames(spatial_coords) <- c("X", "Y")
  local.traj <- local.traj.preprocessing(spatial_coords = spatial_coords, embeddings = nmfs)
  # grid.info <- local.traj$grid.info
  
  # output (spatial coords together with tumor cell identification)
  spatial_coords_tumor <- gather.spatial.tumor(obj = obj, spatial_coords = spatial_coords)
  out.data <- list(local.traj = local.traj, spatial.info = spatial_coords_tumor, nmf.info = nmf.info, 
              sample.info = sample.name, dataset = obj$dataset[1], tumor_type = obj$tumor_type[1], 
              experiment = obj$experiment[1], technology = obj$technology[1])
  saveRDS(out.data, paste0("./output2/", sample.name, ".out.info.rds"))
  return(out.data)
}

file.names <- list.files("./h5_collection")
file.names <- grep("ji2020|NYU|gracia|bergenstra|Gouin", value = T, file.names)
path.names <- paste0("./h5_collection/", file.names)
sample.names <- gsub("\\..*","",file.names)

for (i in 1:length(path.names)){
  tryCatch({
    file.name <- path.names[i]
    sample.name <- sample.names[i]
    run.panc.spatial.grad(file.name = file.name, 
                          sample.name = sample.name, 
                          gene.set = imm.markers)
  }, error=function(e){})
}
