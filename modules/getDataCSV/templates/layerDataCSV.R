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
print(args_list)
INFILE = args_list[["INFILE"]]
print(INFILE)


print("Load packages...")
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)


# DIR = "../"
# INFILE = paste0(DIR, "/data/merged_projected_k32_all_cells_refined_v1.rds")#
dat = readRDS(INFILE)
print(dat)
object.size(dat)


DefaultAssay(dat) <- "Xenium"
# dat <- JoinLayers(dat, assay = "Xenium")
# print(dat)
assay_name = "Xenium"
layer_names <- Layers(dat[[assay_name]])
dims_by_layer <- lapply(layer_names, function(l) {
  x <- GetAssayData(dat, assay = assay_name, layer = l)
  c(features = nrow(x), cells = ncol(x))
})
names(dims_by_layer) <- layer_names
print(dims_by_layer)


# Save layer normalized data to csv files:
assay_name <- "Xenium"  # change if needed
assay <- dat[[assay_name]]
dir.create("layer_csvs", showWarnings = FALSE)
layer_names <- Layers(assay_name)[grepl("^data\\.[0-9]+$", Layers(assay_name))]
dir.create(paste0(DIR, "/data/layer_csvs/"), recursive=T, showWarnings=F)
for (l in layer_names) {
  x <- GetAssayData(dat, assay = assay_name, layer = l)
  # Write sparse to long format: feature, cell, value
  # Works if x is dgCMatrix; if not, it will still run but may be huge.
  if (inherits(x, "dgCMatrix")) {
    tt <- summary(x)  # i, j, x
    df <- data.frame(
      feature = rownames(x)[tt$i],
      cell = colnames(x)[tt$j],
      value = tt$x
    )
  } else {
    # Dense fallback (can explode in size)
    df <- as.data.frame(as.matrix(x))
    df$feature <- rownames(x)
    df <- df[, c("feature", setdiff(names(df), "feature"))]
    # Note: this produces a wide CSV (often too large); adjust if needed.
  }
  out <- file.path("layer_csvs/", paste0(assay_name, "_", l, ".csv.gz"))
  write.table(df, out, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


print("Save sessionInfo ... ")
out = capture.output(sessionInfo())
writeLines(out, "layer_csvs/sessionInfo.txt")

