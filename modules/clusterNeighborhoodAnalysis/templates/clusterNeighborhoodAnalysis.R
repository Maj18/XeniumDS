#!/opt/conda/bin/R

args = commandArgs(trailingOnly = TRUE)
# Parse key=value pairs
args_list = list()
i = 1
while (i <= length(args)) {
  if (startsWith(args[i], "--")) {
    key = substring(args[i], 3)           # remove leading "--"
    value = args[i + 1]                   # next argument is the value
    args_list[[key]] = value
    i = i + 2
  } else {
    i = i + 1
  }
}
OUTDIR = args_list[["OUTDIR"]]
INFILE = args_list[["INFILE"]]

print(OUTDIR)
print(INFILE)

# OUTDIR = "/mnt/results/TEST/AtoMx/Clustering/"
# dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
# INFILE = "/mnt/results/AtoMx_based/QC/All_filtered.RDS"


# https://satijalab.org/seurat/articles/visiumhd_analysis_vignette#identifying-spatially-defined-tissue-domains
# https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md
# https://github.com/prabhakarlab/Banksy
# https://www.sc-best-practices.org/spatial/neighborhood.html
# https://divingintogeneticsandgenomics.com/post/how-to-do-neighborhood-cellular-niches-analysis-with-spatial-transcriptome-data/
# https://prism-oncology.github.io/sopa/tutorials/spatial/

library(Seurat)
library(SeuratData)
library(ggplot2)
library(gridExtra)
library(pals)
library(BiocParallel)
library(dplyr)
library(patchwork)
library(tidyr)
library(tidyverse)
library(readr)
library(readxl)
library(openxlsx)
library(Seurat)
library(Matrix)
library(RANN)
library(dplyr)
library(Matrix)
library(igraph)


# Kelly pallette for visualization
mypal = kelly()[-1]

################# Import data ##################
# INFILE = "INFILE/merged_projected_k32_all_cells_refined_v1.rds"
dat = readRDS(INFILE)
print(dat)
print(colnames(dat@meta.data))


################# Celltype neighborhood analysis ##################
# https://prism-oncology.github.io/sopa/tutorials/spatial/
# library(Seurat)
# library(Matrix)
# library(RANN)
# library(dplyr)

#' Get spatial neighborhood adjacency per sample
getSpatialNeighborAdjacency = function(seurat_obj, k=50,
        coord_cols = c("x", "y"),  
        radius_max = 50,
        radius_min = 0  # same as sopa radius=[0,50] (we'll exclude self by using >0)
    ){
    img_id <- names(seurat_obj@images)[1]  # pick the correct image/slice
    print(img_id)
    cells <- colnames(seurat_obj)
    n <- length(cells)
    print(n)

    # k <- 50  # only used to speed up search; we'll still filter by radius
    # ---- extract coordinates aligned to Seurat cell names ----
    coords_df <- seurat_obj@images[[img_id]]@boundaries$centroids@coords %>%
            as.data.frame()
    rownames(coords_df) = seurat_obj@images[[img_id]]@boundaries$centroids@cells
    coords_df <- coords_df[colnames(seurat_obj), , drop = FALSE]
    print(dim(coords_df))

    # x <- coords_df[, coord_cols[1]]
    # y <- coords_df[, coord_cols[2]]
    # xy <- cbind(x, y)
    xy = coords_df[, coord_cols, drop = FALSE]
    print(head(xy))

    # ---- radius neighbors ----
    # RANN gives neighbors within a max distance; then we drop <= radius_min
    nn <- RANN::nn2(data = xy, query = xy, searchtype = "radius",
                    radius = radius_max, k = k)
    # nn$idx is a matrix of neighbor indices (with -1 for none)
    idx <- nn$nn.idx
    print(head(idx))
    # n <- ncol(seurat_obj)  # number of cells (matches colnames)
    # Build edges (i -> j) for all neighbors j where distance > radius_min
    # We'll compute distances only for found neighbors to filter by radius_min.
    # (To keep it efficient, use the returned neighbor distances if available.)
    # nn$dist returns distances corresponding to nn$idx.
    dist <- nn$nn.dist
    print(head(dist))
    # convert (row i, neighbor col idx[i, ] ) to edge list
    rows <- rep(seq_len(n), times = ncol(idx)) # it goes through cell list first
    cols <- as.vector(idx)
    keep <- (cols!=0) & (as.vector(dist)<radius_max) & (as.vector(dist)>0)   # drop padded "no neighbor"
    # as.vector(dist)[keep] %>% max()
    rows <- rows[keep]
    cols <- cols[keep]
    # Create adjacency (0/1). Then symmetrize to make it undirected connectivity.
    adj <- Matrix::sparseMatrix(i = rows, j = cols, x = 1, dims = c(n, n))
    # there is an edge between i and j if at least one direction exists. That is not “mutual neighbors”
    adj <- (adj + t(adj)) > 0
    adj <- Matrix::drop0(adj)   # optional cleanup
    print(adj[1:10, 1:10])
    # adj is now your “spatial_connectivities”-like graph
    # Row/col names (optional, but useful)
    rownames(adj) <- colnames(seurat_obj)
    colnames(adj) <- colnames(seurat_obj)
    # saveRDS(adj, "adj.RDS", compress=TRUE)

    return(adj)
} 

# library(Matrix)
# library(igraph)
# hop_distance: shortest-path length
mean_hop_distance <- function(adj, cell_type, title="Dynamic Distance Heatmap", figurePath) {
    # adj: sparse adjacency (n x n), binary or weighted (we only need edges)
    # cell_type: named vector length n
    # Ensure ordering consistent with adj
    if (!is.null(colnames(adj))) {
        cell_type <- cell_type[colnames(adj)]
    }
    print(head(cell_type))
    print("message1")
    g <- graph_from_adjacency_matrix(adj != 0, mode = "undirected", diag = FALSE)
    print(g)
    types <- sort(unique(cell_type))
    print(types)
    print("message2")
    # Order cells to align with cell_type
    cells_by_type = lapply(types, function(tp) which(cell_type==tp))
    names(cells_by_type) = types
    # Compute all-pairs distances once (|V|x|V|)
    dist_all = distances(g, v=V(g), to=V(g))
    # print(dist_all)
    # Aggregate:
    print("message3")
    D <- matrix(NA_real_, nrow = length(types), 
        ncol = length(types), dimnames = list(types, types))
    # print(D)
    print("message4")
    for (src in types) {
        src_cells = cells_by_type[[src]]
        # src_cells <- which(cell_type == src)
        # For each target type, we need distance to nearest target cell
        print(src)
        for (tgt in types) {
            tgt_cells = cells_by_type[[tgt]]
            # tgt_cells <- which(cell_type == tgt)
            # Single-source shortest paths from all tgt cells at once
            # Use multi-source distances by creating a super-source via union:
            # Approach: compute distances from tgt subset, then take min for src cells.
            # igraph doesn't directly do multi-source min hops, but we can do:
            # dist_mat <- distances(g, v = tgt_cells, to = src_cells)  # |tgt| x |src|
            dist_mat = dist_all[tgt_cells, src_cells, drop=FALSE]
            # For each src cell, nearest tgt distance:
            nearest_hops <- apply(dist_mat, 2, min)
            # mean over src cells
            D[src, tgt] <- mean(nearest_hops, na.rm = TRUE)
        }
    }
    print("message6")
    # print(D[1:5, 1:5])
    W <- 1 / (D + 1e-8)
    diag(W) <- 0
    # print(W[1:5, 1:5])
    base_width  = 5
    base_height = 5
    scale_factor = 0.3 # Inches to add per element
    pdf_width  = base_width + (ncol(D) * scale_factor)
    pdf_height = base_height + (nrow(D) * scale_factor)
    print(pdf_width)
    print(pdf_height)
    D_plot <- mean_matrix
    finite_vals <- D_plot[is.finite(D_plot)]
    max_finite <- max(finite_vals, na.rm = TRUE)
    D_plot[!is.finite(D_plot)] <- max_finite * 1.1
    print("message7")
    pdf(figurePath, h=pdf_height, w=pdf_width)
        print(pheatmap::pheatmap(
            mat = D_plot,
            # filename = "distance_heatmap.pdf",
            # width = pdf_width,
            # height = pdf_height,
            cluster_rows = TRUE,         # Automatically clusters similar rows
            cluster_cols = TRUE,         # Automatically clusters similar columns
            color = colorRampPalette(c("navy", "white", "firebrick3"))(50), # Custom color gradient
            display_numbers = TRUE,      # Set to FALSE if your matrix is too large for text numbers
            number_color = "black",
            main = title
        ))
    dev.off()
    print("message8")
    return(list(Weight=W, Distance=D))
}

clean_celltype <- function(x) {
  # Make a safe name: replace spaces/punctuation with _
  y <- gsub("[^A-Za-z0-9]+", "_", x)   # anything not alnum -> _
  y <- gsub("^_+|_+$", "", y)         # trim leading/trailing _
  y
}

# Calculate spatial connectivity, mean hop distance and plot separately for each sample
getAdj = function(celltypes, OUTDIR, dat, name="Celltype") {
    lapply(celltypes, function(celltype) {
        print(celltype)
        dir.create(paste0(OUTDIR, "/", clean_celltype(celltype)), recursive=T, showWarnings=F)
        adjs = lapply(unique(dat$Sample), function(sample){
            print(sample)
            sub = subset(dat, Sample==sample)
            Idents(sub) = as.data.frame(sub@meta.data)[[celltype]]
            sub = subset(x = sub, downsample = 1000)
            print(head(clean_celltype(sub[[celltype]])))
            print(head(sub[[celltype]]))
            sub[[celltype]] = clean_celltype(sub[[celltype]])
            print(sub)
            adj = getSpatialNeighborAdjacency(seurat_obj=sub, k=50,
                coord_cols = c("x", "y"), radius_max = 50,
                radius_min = 0)
            print(dim(adj))
            print(adj[1:4, 1:5])
            print(head(sub[[celltype]]))
            cell_type = setNames(as.data.frame(sub@meta.data)[[celltype]], colnames(sub))
            cell_type = cell_type[rownames(adj)]
            print(head(cell_type))
            meanDis = mean_hop_distance(adj, cell_type, 
                title=paste0("/HopDistance_", sample), 
                figurePath=paste0(OUTDIR, "/", celltype, 
                    "/", name, "HopDistance_", sample, "_", celltype, ".pdf"))
            return(meanDis)
        }) %>% setNames(unique(dat$Sample))
        adjs_D = lapply(adjs, function(adj){
            adj$Distance
        }) %>% setNames(paste0(unique(dat$Sample), "_Distance"))
        adjs_W = lapply(adjs, function(adj){
            adj$Weight
        }) %>% setNames(paste0(unique(dat$Sample), "_Weight"))
        saveRDS(adjs, paste0(OUTDIR, "/", clean_celltype(celltype), "/", 
            name, "AdjacentMatrices_", clean_celltype(celltype), ".RDS"), compress=TRUE)
        # Combine samples, calcualte mean hop distance, and plot
        ## Find the common intersecting Row names across ALL matrices
        common_rows = reduce(lapply(adjs_D, rownames), intersect)
        ## Find the common intersecting Column names across ALL matrices
        common_cols = reduce(lapply(adjs_D, colnames), intersect)
        ## Filter and reorder every matrix using the common intersecting names
        ## This guarantees they all match in size, naming, and sequential order
        cleaned_adjs_D <- lapply(adjs_D, function(mat) {
            mat[common_rows, common_cols, drop = FALSE]
        }) %>% as.data.frame()
        matrix_3d <- array(
            data = unlist(cleaned_adjs_D), 
            dim = c(nrow(cleaned_adjs_D[[1]]), ncol(cleaned_adjs_D[[1]]), length(cleaned_adjs_D)),
            dimnames = list(rownames(cleaned_adjs_D[[1]]), colnames(cleaned_adjs_D[[1]]), NULL)
        )
        # Calculate the mean for each intersecting spot
        ## margin = c(1, 2) instructs R to collapse the 3rd dimension (layers) by taking the mean
        mean_matrix <- apply(matrix_3d, MARGIN = c(1, 2), FUN = mean, na.rm = TRUE)
        base_width  = 5
        base_height = 5
        scale_factor = 0.3 # Inches to add per element
        pdf_width  = base_width + (ncol(mean_matrix) * scale_factor)
        pdf_height = base_height + (nrow(mean_matrix) * scale_factor)
        D_plot <- mean_matrix
        finite_vals <- D_plot[is.finite(D_plot)]
        max_finite <- max(finite_vals, na.rm = TRUE)
        D_plot[!is.finite(D_plot)] <- max_finite * 1.1
        pdf(paste0(OUTDIR, "/", clean_celltype(celltype), "/", name, "HopDistance_allSamples_", clean_celltype(celltype), ".pdf"), 
            h=pdf_height, w=pdf_width)
            print(pheatmap::pheatmap(
                mat = D_plot,
                # filename = "distance_heatmap.pdf",
                # width = pdf_width,
                # height = pdf_height,
                cluster_rows = TRUE,         # Automatically clusters similar rows
                cluster_cols = TRUE,         # Automatically clusters similar columns
                color = colorRampPalette(c("navy", "white", "firebrick3"))(50), # Custom color gradient
                display_numbers = TRUE,      # Set to FALSE if your matrix is too large for text numbers
                number_color = "black",
                main = "meanHopDistance_allSamples"
            ))
        dev.off()
        openxlsx::write_xlsx(c(adjs_D, adj_mean=mean_matrix),
            file = paste0(OUTDIR, "/", clean_celltype(celltype), "/", 
            name, "AdjacentMatrices_", clean_celltype(celltype), ".xlsx"),
            rowNames = TRUE
        )
    })
}

# Only use analysis_include == TRUE
dat$Sample = dat$sample_id
dat = subset(dat, subset = (analysis_include==TRUE)&(Sample!="ABX8and5"))
print(dat)
celltypes = c("final_cell_type") 
print(celltypes)
getAdj(celltypes, OUTDIR, dat=dat, name="Celltype")



print("Save sessionInfo ... ")
out = capture.output(sessionInfo())
writeLines(out, paste0(OUTDIR, "/sessionInfo_nh.txt"))
