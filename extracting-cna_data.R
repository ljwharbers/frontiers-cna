  ## Author: Luuk Harbers
  ## Date: xxxx-xx-xx
  ## Script for 
  
  ## Load/install packages
  packages = c("data.table", "pbapply")
  sapply(packages, require, character.only = T)
  source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")
  
  nthreads = 38
  
  # Load in snp array probes, and select rows that have freqcnv = FALSE
  ref = fread("/home/luukharbers/Downloads/snp6.na35.remap.hg38.subset.txt.gz")
  ref = ref[freqcnv == FALSE, ]
  
  # Prepare for overlapping
  ref[, start := pos]
  ref[, end := pos + 1]
  setkey(ref, chr, start, end)
  
  # Load in TCGA data
  files = list.files("/mnt/AchTeraD/Documents/Projects/NucleAI/data/TCGA/gdc-datasets/segment_copy-ratios/",
                     pattern = ".txt$", recursive = TRUE, full.names = TRUE)
  files = files[!grepl("annotation", files)]
  files = togo
  # lapply through files 
  invisible(pblapply(files, function(x) {
    dt = fread(x)
    setnames(dt, c("aliquot", "chr", "start", "end", "nprobes", "segment_mean"))
    setkey(dt, chr, start, end)
    
    # Get overlaps
    overlaps = foverlaps(dt, ref)
    overlaps = overlaps[, .(probeid, chr, pos, segment_mean)]
    
    # Order 
    overlaps = overlaps[gtools::mixedorder(chr)]
    overlaps = merge(ref[, .(probeid, chr, pos)], overlaps, all.x = T)
    overlaps = unique(overlaps, by = "probeid")
    
    # Set aliquot name
    setnames(overlaps, "segment_mean", dt$aliquot[1])
    write.table(overlaps[, 4], paste0("/mnt/AchTeraD/Documents/Projects/ijms_review/data/tcga-cnv-log2/", dt$aliquot[1], ".tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    rm(c(dt, overlaps))
  
  }, cl = nthreads))

# CBIND IN BASH !!!
##### eval paste $(printf "<(cut -f4 %s) " *.tsv) > log2-combined.tsv ######
##### paste *.tsv > log2-combined.tsv ######  

# Load data
total = fread("/mnt/AchTeraD/Documents/Projects/ijms_review/data/tcga-cnv-log2/combined-log2-ratios.tsv.gz",
              header = TRUE)

  # Remove rows that contain NAs
total = na.omit(total)


# Make umap
config = umap.defaults
config$n_neighbors = 50
#config$min_dist = 0.1

total_umap = umap(t(total), method = "umap-learn", config = config)

umap_dt = data.table(x = total_umap$layout[, 1],
                     y = total_umap$layout[, 2],
                     sample = colnames(total[, 4:ncol(total)]))

plot = ggplot(umap_dt, aes(x = x, y = y)) +
  geom_point(size = 2) +
  scale_color_npg() +
  labs(x  = "UMAP 1", y = "UMAP 2")

