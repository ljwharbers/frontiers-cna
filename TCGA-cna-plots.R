## Author: Luuk Harbers
## Date: 2020-12-07
## Script for analysis of TCGA data

## Load/install packages
packages = c("data.table", "pbapply", "ggplot2", "TCGAutils", "RColorBrewer")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Set threads
nthreads = 40

# Load and list files
clinical = readRDS("/mnt/AchTeraD/Documents/Projects/ijms_review/data/clinical.rds")
tss = fread("/mnt/AchTeraD/Documents/Projects/ijms_review/data/tcga-tss-codetable-abb.tsv")
chr_annot = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/hg19_chromosomal_arms.bed")
cosmic = fread("/mnt/AchTeraD/Documents/Projects/ijms_review/data/cosmic/gene-census.tsv")
number = fread("/mnt/AchTeraD/Documents/Projects/ijms_review/data/number-of-samples.tsv")

files = list.files("/mnt/AchTeraD/Documents/Projects/ijms_review/data/tcga-cnv-log2/", 
                   recursive = T, full.names = TRUE, pattern = ".txt$")
files = files[!grepl("annotation", files)]

# Get chrom lengths
chrom_lengths = data.table(chr = c(1:22, "X"), length = c(249250621, 243199373, 198022430, 191154276, 180915260,
                                                          171115067, 159138663, 146364022, 141213431, 135534747,
                                                          135006516, 133851895, 115169878, 107349540, 102531392,
                                                          90354753, 81195210, 78077248, 59128983, 63025520, 48129895,
                                                          51304566, 155270560))
# Manual ordering function
# Get order function
position_orderdodge = function(
  width = NULL
) {
  ggproto(NULL, PositionOrderdodge, width = width)
}
PositionOrderdodge = ggproto(
  "PositionOrderdodge", PositionDodge,
  compute_panel = function(data, params, scales) {
    data = flip_data(data, params$flipped_aes)
    w = params$width
    data = data %>% group_by(x, group) %>%
      mutate(rank = scales::rescale(rank(y), to = c(-w, w)))
    data = as.data.frame(data)
    data$x = data$x + data$rank
    data$rank = NULL
    flip_data(data, params$flipped_aes)
  }
)

# Set amp/del threshold
amp = log2(2.5/2)
del = log2(1.5/2)

# Remove chr from chr annot file
chr_annot[, V1 := gsub("chr", "", V1)]
chr_annot[, arm_size := V3 - V2]

# apply through files and retrieve relevant information
result = pblapply(files, function(sample) {
  dt = fread(sample)
  
  # Annotate with chromosomal arms
  setkey(dt, Chromosome, Start, End)
  setkey(chr_annot, V1, V2, V3)
  
  dt = foverlaps(dt, chr_annot)
  dt[, Length := End - Start]
  
  # Annotate segments
  dt[!(Segment_Mean >= amp | Segment_Mean <= del), alteration_type := "neutral"]
  dt[Length < 1e4 & (Segment_Mean >= amp | Segment_Mean <= del), alteration_type := "indel"]
  dt[Length >= .75 * arm_size & (Segment_Mean >= amp | Segment_Mean <= del), alteration_type := "aneuploidy"]
  dt[Length < .75 * arm_size & Length > 1e4 & (Segment_Mean >= amp | Segment_Mean <= del), alteration_type := "cna"]
  
  # Reorder, so aneuploidy > cna > indel > neutral and then get unique rows
  setorder(dt, alteration_type)
  dt = unique(dt[, c(6, 1, 7:12), with = F], by = c("Chromosome", "Start", "End"))

  # Get genome length
  genome = dt[, sum(End - Start)]
  
  # Get number of events
  amp_events = dt[Segment_Mean >= amp & alteration_type == "cna", .N]
  del_events = dt[Segment_Mean <= del & alteration_type == "cna", .N]
  
  # Get percentage amp/del
  amp_percentage = dt[Segment_Mean >= amp & alteration_type == "cna", sum(End - Start)] / genome * 100
  del_percentage = dt[Segment_Mean <= del & alteration_type == "cna", sum(End - Start)] / genome * 100
  
  segments = dt[alteration_type != "neutral", .(Chromosome, Start, End, Num_Probes, Segment_Mean, Length, alteration_type)]
  segments[, CNA := ifelse(Segment_Mean >= amp, "AMP", "DEL")]
  
  res = list(aliquot = dt$GDC_Aliquot[1],
             genome = genome,
             n_events = c(amp = amp_events, del = del_events),
             percentage = c(amp = amp_percentage, del = del_percentage),
             segments = segments)
  return(res)
}, cl = nthreads)

# Extract segment
result_segments = pblapply(result, function(sample){
  cbind(data.table(aliquot = sample$aliquot),
        sample$segments)
})

result_segments = rbindlist(result_segments)
result_segments = result_segments[complete.cases(result_segments)]

# Per tumor type
# Get barcode, and match to tumor type
barcodes = UUIDtoBarcode(unique(result_segments$aliquot), from_type = "aliquot_ids")
setDT(barcodes)
barcodes[, tss := gsub("TCGA-|-.*", "", portions.analytes.aliquots.submitter_id)]
barcodes = merge(barcodes, tss[, c("TSS Code" , "Study Abbreviation")], by.x = "tss", by.y = "TSS Code")

# Merge with results
result_segments = merge(result_segments, barcodes[, c("portions.analytes.aliquots.aliquot_id", "Study Abbreviation")], 
                    by.x = "aliquot", by.y = "portions.analytes.aliquots.aliquot_id")

# Calculate percentage per subtype
alt_types_subtype = result_segments[, .N, by = c("alteration_type", "Study Abbreviation")]
alt_types_subtype = alt_types_subtype[, list(alteration_type = alteration_type, N = N, total = sum(N)), by = "Study Abbreviation"]
alt_types_subtype[, percentage := N / total * 100]

# Set factor levels
alt_types_subtype[, alteration_type := factor(alteration_type, levels = c("aneuploidy", "indel", "cna"))]
neworder = alt_types_subtype[alteration_type == "cna", ]
setorder(neworder, -percentage)
alt_types_subtype[, `Study Abbreviation` := factor(`Study Abbreviation`, levels = neworder$`Study Abbreviation`)]

# Plot per subtype stacked barplots
plt1 = ggplot(alt_types_subtype, aes(x = `Study Abbreviation`, y = percentage, fill = alteration_type)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(y = "Proportion of alteration type (%)", x = "Tumor type", fill = "Alteration type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/proportion-of-alteration-types",
              width = 14, height = 8)

# Plot distribution of segment lengths
segments_cnas = result_segments[alteration_type == "cna"]
segments_cnas_mean = segments_cnas[, .(Length = mean(Length)), by = .(CNA, `Study Abbreviation`, aliquot)]
setnames(segments_cnas_mean, "Study Abbreviation", "tumor_type")

# Get counts
neworder = segments_cnas_mean[, mean(Length), by = .(tumor_type)]
setorder(neworder, -V1)
segments_cnas_mean[, `Study Abbreviation` := factor(tumor_type, levels = neworder$tumor_type)]

number[, tumor_type := factor(tumor_type, levels = levels(segments_cnas_mean$tumor_type))]

plt0 = ggplot(segments_cnas_mean, aes(x = CNA, y = Length/1e6)) +
  facet_grid(~tumor_type, switch = "x") +
  geom_rect(aes(fill = as.factor(as.numeric(as.factor(tumor_type)) %% 2)),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(position = position_orderdodge(.3), size = 0.5, aes(color = CNA)) +
  geom_text(data=number, aes(y = -2, x = 1.35, label = paste0("n = ", V1)), hjust = 0.5) +
  stat_compare_means(paired = F, label = "p.signif") + 
  scale_color_brewer(palette = "Set1", labels = c("Amplified", "Deleted")) +
  scale_fill_manual(values = c("white", "gray90"), guide = guide_none()) +
  stat_summary(fun = mean, geom = "crossbar", color = "black") +
  labs(y = "Length of alteration (Mb)", x = "", color = "Alteration type") +
  theme(panel.spacing = unit(0, "points"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.placement = "outside")

save_and_plot(plt0, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/meanlength-per-tumortype",
              width = 22, height = 8)

segments_cnas = segments_cnas[, .(Length = mean(Length)), by = .(Chromosome, `Study Abbreviation`)]

# Set factors and normalize for chr length
segments_cnas = merge(segments_cnas, chrom_lengths, by.x = "Chromosome", by.y = "chr")
segments_cnas[, norm_length := Length/length]
segments_cnas[, Chromosome := factor(Chromosome, levels = c(1:22, "X"))]

plt10 = ggplot(segments_cnas, aes(x = Chromosome, y = `Study Abbreviation`, fill = norm_length * 1e2)) +
geom_tile(color = "black") +
scale_fill_viridis(begin = 0.1) +
labs(fill = "Average normalized alteration\nlength (Mb)")
  
save_and_plot(plt10, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/mean-alteration-length-norm",
              width = 14, height = 8)



# Extract Percentage of amp/del
result_perc = pblapply(result, function(sample){
  data.table(aliquot = sample$aliquot,
             amp = sample$percentage["amp"],
             del = sample$percentage["del"])
})

result_perc = rbindlist(result_perc)

# Merge with results
result_perc = merge(result_perc, barcodes[, c("portions.analytes.aliquots.aliquot_id", "Study Abbreviation")], 
                    by.x = "aliquot", by.y = "portions.analytes.aliquots.aliquot_id")

# Get correlation of mean amp/del per tumor type
result_perc_cor = result_perc[, .(mean_amp = mean(amp), mean_del = mean(del)), by = .(`Study Abbreviation`)]

plt2 = ggplot(result_perc_cor, aes(x = mean_amp, y = mean_del)) + 
  geom_point() +
  geom_smooth(method = lm, formula = y ~ x, se = F, color = "red", linetype = 2) +
  stat_cor() +
  labs(y = "Average percentage deleted", x = "Average percentage amplified")

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/correlation-mean-amp_del",
              width = 7, height = 7)

# # Plot scatter
# plt2 = ggplot(result_perc, aes(x = amp, y = del)) +
#   geom_point() +
#   stat_cor() +
#   geom_smooth(method = lm, formula = y ~ x, se = F, color = "red", linetype = 2) +
#   #facet_wrap(~ `Study Abbreviation`) +
#   labs(y = "Genome deleted (%)", x = "Genome amplified (%)")
# 
# save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/correlation-amp_del-all",
#               width = 14, height = 10)


# Melt
result_perc_m = melt(result_perc, id.vars = c("aliquot", "Study Abbreviation"))
setnames(result_perc_m, c("aliquot", "tumor_type", "cna", "value"))


# Set facet factors
neworder = result_perc_m[, mean(value), by = "tumor_type"]
setorder(neworder, -V1)
result_perc_m[, tumor_type := factor(tumor_type, levels = neworder$tumor_type)]

# Get counts
number = result_perc_m[, .N / 2, by = "tumor_type"]
number[, tumor_type := factor(tumor_type, levels = levels(result_perc_m$tumor_type))]

plt3 = ggplot(result_perc_m, aes(x = cna, y = value)) +
  facet_grid(~tumor_type, switch = "x") +
  geom_rect(aes(fill = as.factor(as.numeric(as.factor(tumor_type)) %% 2)),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(position = position_orderdodge(.3), size = 0.5, aes(color = cna)) +
  geom_text(data=number, aes(y = -2, x = 1.35, label = paste0("n = ", V1)), hjust = 0.5) +
  stat_compare_means(paired = F, label = "p.signif") + 
  scale_color_brewer(palette = "Set1", labels = c("Amplified", "Deleted")) +
  scale_fill_manual(values = c("white", "gray90"), guide = guide_none()) +
  stat_summary(fun = mean, geom = "crossbar", color = "black") +
  labs(y = "Genome altered (%)", x = "", color = "Alteration type") +
  theme(panel.spacing = unit(0, "points"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.placement = "outside")

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/Percentage-altered-dots",
              width = 22, height = 8)


# Extract segment counts
result_events = pblapply(result, function(sample){
  data.table(aliquot = sample$aliquot,
             amp = sample$n_events["amp"],
             del = sample$n_events["del"])
})

result_events = rbindlist(result_events)

# Merge with results
result_events = merge(result_events, barcodes[, c("portions.analytes.aliquots.aliquot_id", "Study Abbreviation")], 
                    by.x = "aliquot", by.y = "portions.analytes.aliquots.aliquot_id")

# Melt
result_events_m = melt(result_events, id.vars = c("aliquot", "Study Abbreviation"))
setnames(result_events_m, c("aliquot", "tumor_type", "cna", "value"))


# Set facet factors
neworder = result_events_m[, mean(value), by = "tumor_type"]
setorder(neworder, -V1)
result_events_m[, tumor_type := factor(tumor_type, levels = neworder$tumor_type)]

# Get counts
number = result_events_m[, .N / 2, by = "tumor_type"]
number[, tumor_type := factor(tumor_type, levels = levels(result_events_m$tumor_type))]

plt4 = ggplot(result_events_m, aes(x = cna, y = value)) +
  facet_grid(~tumor_type, switch = "x") +
  geom_rect(aes(fill = as.factor(as.numeric(as.factor(tumor_type)) %% 2)),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(position = position_orderdodge(.3), size = 0.5, aes(color = cna)) +
  geom_text(data=number, aes(y = -80, x = 1.35, label = paste0("n = ", V1)), hjust = 0.5) +
  scale_color_brewer(palette = "Set1", labels = c("Amplified", "Deleted")) +
  scale_fill_manual(values = c("white", "gray90"), guide = guide_none()) +
  stat_compare_means(paired = F, label = "p.signif") + 
  stat_summary(fun = mean, geom = "crossbar", color = "black") +
  labs(y = "Number of events", x = "", color = "Alteration type") +
  theme(panel.spacing = unit(0, "points"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.placement = "outside")

save_and_plot(plt4, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/num_events-altered-dots",
              width = 22, height = 8)

plt5 = ggplot(result_events_m, aes(x = cna, y = value)) +
  facet_grid(~tumor_type, switch = "x") +
  geom_rect(aes(fill = as.factor(as.numeric(as.factor(tumor_type)) %% 2)),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(position = position_orderdodge(.3), size = 0.5, aes(color = cna)) +
  geom_text(data=number, aes(y = -80, x = 1.35, label = paste0("n = ", V1)), hjust = 0.5) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_manual(values = c("white", "gray90"), guide = guide_none()) +
  stat_summary(fun = mean, geom = "crossbar", color = "black") +
  scale_y_continuous(limits = c(-80, 2000)) +
  labs(y = "Number of events", x = "", color = "Alteration type") +
  theme(panel.spacing = unit(0, "points"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.placement = "outside")

save_and_plot(plt5, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/num_events-altered-dots-2Kmax",
              width = 22, height = 8)

# Get info about lengths and locations
# Get counts per chromosome per tumor type get average counts per tumor type and chromosome
segment_counts = result_segments[, .N, by = c("Study Abbreviation", "Chromosome")]
segment_counts = segment_counts[, list(Chromosome = Chromosome, prop = N / sum(N)),  by = "Study Abbreviation"]
segment_counts = segment_counts[complete.cases(segment_counts)]



# Normalize for chr length
segment_counts = merge(segment_counts, chrom_lengths, by.x = "Chromosome", by.y = "chr")
segment_counts[, norm_prop := prop/length]

# Set factors
segment_counts[, Chromosome := factor(Chromosome, levels = c(1:22, "X"))]

plt6 = ggplot(segment_counts, aes(x = Chromosome, fill = Chromosome, y = norm_prop)) +
  facet_wrap(. ~ `Study Abbreviation`, ncol = 5) +
  geom_col() +
  labs(y = "Normalized proportion of total events") +
  scale_fill_viridis_d() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

save_and_plot(plt6, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/proportion_events-chr_tumortype",
              width = 18, height = 10)

# Plot as heatmap
plt7 = ggplot(segment_counts, aes(x = Chromosome, y = `Study Abbreviation`, fill = norm_prop *1e9)) +
  geom_tile(color = "black") +
  labs(fill = "Proportion of total events\n(proportion / chromosome length) * 1e9") +
  scale_fill_viridis(begin = 0.1, limits = c(0, 1.5), oob = scales::squish)

save_and_plot(plt7, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/proportion_events-chr_tumortype-heatmap-scale_limitedOOB",
              width = 16, height = 10)

# COSMIC gene analysis
# Set keys
setkey(result_segments, Chromosome, Start, End)
setkey(cosmic, chr, start, end)

result_genes = foverlaps(result_segments[alteration_type == "cna"], cosmic)
result_genes = result_genes[complete.cases(result_genes)]

# Merge with number of cancers for normalization
gene_counts = result_genes[, .N, by = c("CNA", "gene", "Study Abbreviation")]
gene_counts = merge(gene_counts, number, by.x = "Study Abbreviation", by.y = "tumor_type")
gene_counts[, alt_frequency := N / V1]
gene_counts[CNA == "DEL", alt_frequency := alt_frequency * -1]

gene_counts[, id := paste(gene, CNA, sep ="-")]

order = gene_counts[, sum(alt_frequency), by = c("id", "gene", "CNA")]
setorder(order, CNA, -V1)
select = c(head(order$id, 25), tail(order$id, 25))

gene_counts = gene_counts[id %in% select]
gene_counts[, gene := factor(gene, levels = c(head(order$gene, 25), tail(order$gene, 25)))]

# Select colors
n = 32
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Plot
plt8 = ggplot(gene_counts, aes(x = gene, y = alt_frequency, fill = `Study Abbreviation`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt8, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/cosmic_CNAgenes-subtype-stratified",
              width = 12, height = 8)

gene_counts_all = result_genes[, .N, by = c("CNA", "gene")]
gene_counts_all[, freq := N / sum(number$V1)]
gene_counts_all[CNA == "DEL", freq := freq * -1]
gene_counts_all[, id := paste(gene, CNA, sep = "-")]

setorder(gene_counts_all, -freq)
select = c(head(gene_counts_all$id, 25), tail(gene_counts_all$id, 25))

gene_counts_all = gene_counts_all[id %in% select]
gene_counts_all[, gene := factor(gene, levels = gene)]

plt9 = ggplot(gene_counts_all, aes(x = gene, y = freq, fill = CNA)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt9, "/mnt/AchTeraD/Documents/Projects/ijms_review/Plots/cosmic_CNAgenes-total",
              width = 12, height = 8)
