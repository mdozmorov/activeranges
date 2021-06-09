### =========================================================================
### activeranges is an AnnotationHub package that stores genomic coordinates of
### DNAse I hypersensitive regions (open/active chromatin) as GRanges. 
### -------------------------------------------------------------------------
###

# activeranges package currently contains 16 Rds objects with cell/tissue-specific GRanges.
# The original data are available at 
# https://www.meuleman.org/research/dhsindex/

# The object names are structured as "<genome assembly>.<lab>.<original file name>", 
# e.g., "hg38.Meuleman.Cancer_epithelial".

# Get the full data
# wget https://www.meuleman.org/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz
  
# The following example demonstrate how the coordinates of cell/tissue-specific 
# DNAse I regions were converted into Rds objects

# Download results data from 

library(GenomicRanges)
library(readr)
library(ggplot2)
# Folder with results
dir_in <- "/Volumes/GoogleDrive/My Drive/activedata"
file <- "DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz"
mtx <- read_tsv(file.path(dir_in, file), col_types = c("cddcdddddc"))

# Unique cell/tissue types
cell_types <- unique(mtx[["component"]])

# Add seqinfo
genome_id <- "hg38"
# Get chromosome info and match it to the chromosome order in denyBED
chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
chrom_data <- chrom_data[chrom_data$chrom %in% seqlevels(activeGR), ]
chrom_data <- chrom_data[match(seqlevels(activeGR), chrom_data$chrom), ]

# Extract GRanges for each cell type
for (cell_type in cell_types) {
  # Subset to cell/tissue type
  mtx_subset <- mtx[mtx[["component"]] == cell_type, ]
  # Convert to GRanges object
  activeGR <- GenomicRanges::makeGRangesFromDataFrame(mtx_subset, keep.extra.columns = TRUE)
  # Check if chromosome order is the same
  if (!all.equal(seqlevels(activeGR), chrom_data$chrom)) {
    print(paste("Chromosome order does not match for", genome_id, "genome."))
    break
  }
  # Assign seqinfo data
  seqlengths(activeGR) <- chrom_data$size
  isCircular(activeGR) <- chrom_data$circular
  genome(activeGR)     <- genome_id
  
  # Construct file name, from "Organ devel. / renal"
  # to "hg38.Meuleman.Organ_devel_renal.rds"
  fileout <- paste0(genome_id, ".Meuleman.", gsub("[./ ]+", "_", cell_type, perl = TRUE), ".rds")
  # Save as Rds object
  saveRDS(object = activeGR, file = file.path(dir_in, fileout))
}






