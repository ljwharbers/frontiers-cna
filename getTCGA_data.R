## Author: Luuk Harbers
## Date: 2020-12-07
## Script for downloading TCGA data

## Load/install packages
packages = c("data.table", "GenomicDataCommons", "TCGAutils")
sapply(packages, require, character.only = T)

# Set threads
nthreads = 40

# Query the GDC database for all TCGA copy number files
manifest_segments = files() %>% 
  filter(data_type == "Masked Copy Number Segment") %>% 
  filter(analysis.workflow_type == "DNAcopy") %>%
  filter(cases.samples.sample_type == "primary tumor") %>%
  filter(data_format == "txt") %>%
  manifest()

# Set cachce
gdc_set_cache("/mnt/AchTeraD/Documents/Projects/ijms_review/data/tcga-cnv-log2_new/")

# Download data to gdc_cache directory
download = transfer(manifest_segments$id)

# Get Clinical information
case_ids = UUIDtoUUID(manifest_segments$id, to_type = "case_id")
clinical = gdc_clinical(case_ids = case_ids$cases.case_id)

# Save object
saveRDS(clinical, "/mnt/AchTeraD/Documents/Projects/ijms_review/data/clinical.rds")
