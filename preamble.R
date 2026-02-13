# <!-- Start up loading and configuration.  -->


## General
project_name <- "WTx_BreastCancer"       

## Set up directory structure

parent_dir <- getwd()           # Location of analysis
if(!dir.exists(parent_dir)){
  dir.create(parent_dir, recursive=TRUE)
}
setwd(parent_dir)

analysis_dir <- file.path(parent_dir, "analysis")                # location of analysis_dir
objects_dir <- file.path(analysis_dir, "Robjects")                        # nested within analysis_dir
raw_dir <- file.path(parent_dir, "rawData")
qc_dir <- file.path(analysis_dir, "qc")
norm_dir <- file.path(analysis_dir, "normalization")
dimRed_dir <- file.path(analysis_dir, "dimensionReduction")
niche_dir <- file.path(analysis_dir, "niches")
celltyping_dir <- file.path(analysis_dir, "cellTyping")
tumorDE_dir <- file.path(analysis_dir, "tumorDE")
pathway_dir <- file.path(analysis_dir, "pathways")
immNeigh_dir <- file.path(analysis_dir, "immuneNeighborhoods")
vessel_dir <- file.path(analysis_dir, "vessels")
metabolism_dir <- file.path(analysis_dir, "metabolism")
prolif_dir <- file.path(analysis_dir, "prolif")

# Create directories
for(dir_to_create in c(analysis_dir, objects_dir, raw_dir)){
  dir.create(dir_to_create, recursive = TRUE, showWarnings = FALSE)
}


for (dir_path in c(tumorDE_dir, pathway_dir, vessel_dir, 
                   qc_dir, norm_dir, dimRed_dir, immNeigh_dir, 
                   metabolism_dir, niche_dir, celltyping_dir,
                   prolif_dir)){
  dir.create(file.path(dir_path, "Data"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_path, "Plots"), recursive = TRUE, showWarnings = FALSE)
}

# Key standard files
results_list_file <- file.path(analysis_dir, "results_list.rds")          # list object that keeps track of small results like number of cells.



#| label: libraries

# Set up environment

# List of CRAN packages for session
.cran_packages <- c("data.table", "plyr", "dplyr", "ggplot2", "knitr", "Seurat", "viridis", 
                    "yaml", "BiocManager", "devtools", "ggtangle",
                    "RColorBrewer", "msigdbr", "pheatmap", "GSEABase", "AUCell",
                    "doMC", "dbscan", "nsprcomp", "tidyr", "pals", "ggrepel",
                    "EnvStats", "mclust", "scales", "tibble", "grid", "Matrix",
                    "circlize", "dendextend", "writexl")

# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if(!("remotes" %in% installed.packages())){
  install.packages('remotes')
}
if(!("ggtangle" %in% installed.packages())){
  install.packages("https://cran.r-project.org/src/contrib/ggtangle_0.0.6.tar.gz", type = "source") # Default install is 0.0.3 which doesn't work with clusterProfiler
}

.inst <- .cran_packages %in% installed.packages()
if(length(.cran_packages[!.inst]) > 0) remotes::install_cran(.cran_packages[!.inst])

cran_loaded <- lapply(.cran_packages, require, character.only=TRUE)

## If Rcpp gives trouble, had success with the following:
# remotes::install_version("Rcpp", version = "1.0.12")

# Install Bioconductor packages (if not already installed)
.bioc_packages <- c('InteractiveComplexHeatmap', 'ComplexHeatmap', 
                    'SingleCellExperiment', 'clusterProfiler',
                    'ReactomePA', 'org.Hs.eg.db', 'fgsea',
                    'zellkonverter')
.bioc_inst <- .bioc_packages %in% installed.packages()
if(length(.bioc_packages[!.bioc_inst]) > 0){
  BiocManager::install(.bioc_packages[!.bioc_inst])
}
bioc_loaded <- lapply(.bioc_packages, require, character.only = TRUE)

# Load local packages
.local_packages <- c("smiDE", "CellularNeighborhoods", "HieraType")

if(!("smiDE" %in% installed.packages())){
  remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space",
                          subdir = "_code/smiDE", ref = "Main")
}

if(!("CellularNeighborhoods" %in% installed.packages())){
  devtools::install("local_packages/CellularNeighborhoods/")
}

if(!("HieraType" %in% installed.packages())){
  remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space", 
                          subdir = "_code/HieraType", ref = "Main")
}

other_loaded <- lapply(.local_packages, require, character.only=TRUE)

# Load in smiSmooth and scPearsonPCA directly
devtools::load_all("local_packages/scPearsonPCA/")

## Python packages
if(!reticulate::py_module_available("novae")){
  reticulate::py_install("novae")
}
reticulate::py_require('novae')

if(file.exists(results_list_file)){
  message("Reading in saved results list.")
  results_list <- readRDS(results_list_file)
} else {
  message("Creating results list from user-defined parameters.")
  results_list <- list(
    
  )
  saveRDS(results_list, results_list_file)
}


#### Set up functions for use in project ####
make_umap <- function(pcaobj_embeddings, min_dist=0.01, n_neighbors=30, metric="cosine",key ="UMAP_", spread = 2){
  ump <- 
    uwot::umap(pcaobj_embeddings
               ,n_neighbors = n_neighbors
               ,nn_method = "annoy"
               ,metric = metric
               ,min_dist = min_dist
               ,spread = spread
               ,ret_extra = c("fgraph","nn")
               ,verbose = TRUE)
  
  umpgraph <- ump$fgraph
  dimnames(umpgraph) <- list(rownames(ump$nn[[1]]$idx), rownames(ump$nn[[1]]$idx))
  colnames(ump$embedding) <- paste0(key, c(1,2)) 
  ump <- Seurat::CreateDimReducObject(embeddings = ump$embedding, key = key)
  return(list(grph = umpgraph
              ,ump = ump))
}

# simple function to put files in numeric order (not alphabetic order)
sortFiles <- function(x, i){
  cluster_number <- as.numeric(gsub(".png", "", unlist(lapply(strsplit(basename(x), split="_"), "[[", as.integer(i)))))
  out <- data.frame(file=x, cluster=cluster_number) %>% arrange(cluster) %>% pull(file)
  return(out)
}

# Add number of cells to a heatmap derived from Hieratype
#' show number of cells in heatmap plot
add_ncells <- function(heatmapp){
  newp <- copy(heatmapp)
  lvs <- levels(newp$data$cluster)
  havlvs <- intersect(lvs, unique(as.character(newp$data$cluster)))
  newlvs <- newp$data[,head(.SD,1),by=cluster]
  newlvs <- newlvs[match(havlvs, cluster),paste0(cluster, " (", scales::comma(ncells), ")")]
  newp$data[,cluster:=paste0(cluster, " (", scales::comma(newp$data$ncells), ")")]
  newp$data[,cluster:=factor(cluster, levels=newlvs)]
  return(newp)
}

## Set up scale bar function

make_scale_bar_r <- function(x_vals, y_vals) {
  
  # Adds a scale bar to a ggplot.
  
  #Parameters:
  #x_vals: vector of x coordinates in mm
  #y_vals: vector of y coordinates in mm
  
  # Example usage:
  # scale_bar = make_scale_bar_r(x_vals = cell_meta$x_slide_mm, y_vals = cell_meta$y_slide_mm)
  # ggplot() + scale_bar$bg + scale_bar$rect + scale_bar$label
  
  
  # Calculate x-axis range
  x_range <- range(x_vals, na.rm = TRUE)
  x_length_um <- diff(x_range)*1000 # um to mm conversion

  # Target scale length ~1/4 of the x-axis
  target <- x_length_um / 4
  
  # Compute order of magnitude
  order <- 10^floor(log10(target))
  mantissa <- target / order
  
  # Round mantissa to nearest 1, 2, or 5
  nice_mantissa <- if (mantissa < 1.5) {
    1
  } else if (mantissa < 3.5) {
    2
  } else if (mantissa < 7.5) {
    5
  } else {
    10
  }
  
  # Final scale length in mm
  scale_length_um <- nice_mantissa * order
  scale_length_mm <- scale_length_um/1000
  
  # Format label
  scale_label <- if (scale_length_um >= 1000) {
    paste0(scale_length_um / 1000, " mm")
  } else {
    paste0(scale_length_um, " µm")
  }
  
  # Set coordinates for the scale bar 
  x_start <- x_range[2] - scale_length_mm * 1.1
  x_end <- x_range[2] - scale_length_mm * 0.1
  y_pos <- min(y_vals, na.rm = TRUE) + scale_length_mm * 0.1
  
  # Generate scale bar background, scale bar and annotation to return
  list(
    bg = annotation_custom(
      grob = rectGrob(gp = gpar(fill = "white", alpha = 0.8, col = NA)),
      xmin = x_start- scale_length_mm*0.05, xmax = x_end+ scale_length_mm*0.05,
      ymin = y_pos- scale_length_mm*0.05, ymax = y_pos + scale_length_mm*0.3
    ),
    rect = annotation_custom(
      grob = rectGrob(gp = gpar(fill = "black")),
      xmin = x_start, xmax = x_end,
      ymin = y_pos, ymax = y_pos + scale_length_mm * 0.05
    ),
    label = annotation_custom(
      grob = textGrob(scale_label, gp = gpar(col = "black"), just = "center", vjust = 0),
      xmin = (x_start + x_end)/2, xmax = (x_start + x_end)/2,
      ymin = y_pos + scale_length_mm * 0.1, ymax = y_pos + scale_length_mm * 0.1
    )
  )
}
