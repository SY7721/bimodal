
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(AnnotationDbi)

nn <- 1 # 1, 2, 3, 4
extraCols_narrowPeak <- c(
  signalValue = "numeric",
  pValue      = "numeric",
  qValue      = "numeric",
  peak        = "integer"
)
if (nn == 1){
  group1_name <- "cross-talk"
  group2_name <- "two-state"
  
  h3k4me3_gr <- import(
    "peak/ENCFF357JNZ.bed.gz",
    format = "BED",
    extraCols = extraCols_narrowPeak
  )
sheet_c57 <- read_excel("gene name.xlsx", sheet = "F c57")
sheet_cast <- read_excel("gene name.xlsx", sheet = "F cast")
extract_genes <- function(sheet) {
  telegraph <- sheet[[1]][!is.na(sheet[[1]])]
  crosstalk <- sheet[[3]][!is.na(sheet[[3]])]
  return(list(telegraph = unique(telegraph), 
              crosstalk = unique(crosstalk)))
}
c57 <- extract_genes(sheet_c57)
cast <- extract_genes(sheet_cast)
my_genes <- unique(c(c57$crosstalk, cast$crosstalk))     # crosstalk
my_genes1 <- unique(c(c57$telegraph, cast$telegraph))  # telegraph
} else if (nn == 2){
  group1_name <- "cross-talk"
  group2_name <- "two-state"
  
  h3k4me3_gr <- import(
    "peak/ENCFF671UNN.bed",
    format = "BED",
    extraCols = extraCols_narrowPeak
  )
  sheet_c57 <- read_excel("gene name.xlsx", sheet = "E c57")
  sheet_cast <- read_excel("gene name.xlsx", sheet = "E cast")
extract_genes <- function(sheet) {
  telegraph <- sheet[[1]][!is.na(sheet[[1]])]
  crosstalk <- sheet[[3]][!is.na(sheet[[3]])]
  return(list(telegraph = unique(telegraph), 
              crosstalk = unique(crosstalk)))
}
c57 <- extract_genes(sheet_c57)
cast <- extract_genes(sheet_cast)
my_genes <- unique(c(c57$crosstalk, cast$crosstalk))     # crosstalk
my_genes1 <- unique(c(c57$telegraph, cast$telegraph))  # telegraph}
} else if (nn == 3){
  group1_name <- "bimodal"
  group2_name <- "unimodal"
  
  h3k4me3_gr <- import(
    "peak/ENCFF357JNZ.bed.gz",
    format = "BED",
    extraCols = extraCols_narrowPeak
  )
sheet_c57 <- read_excel("gene name.xlsx", sheet = "F c57")
sheet_cast <- read_excel("gene name.xlsx", sheet = "F cast")
extract_genes <- function(sheet) {
  bimodal <- c(sheet[[1]][!is.na(sheet[[1]])], 
               sheet[[3]][!is.na(sheet[[3]])])
  unimodal <- setdiff(sheet[[5]][!is.na(sheet[[5]])],bimodal)
  return(list(bimodal = unique(bimodal), 
              unimodal = unique(unimodal)))
}
c57 <- extract_genes(sheet_c57)
cast <- extract_genes(sheet_cast)
my_genes <- unique(c(c57$bimodal, cast$bimodal))     # bimodal
my_genes1 <- unique(c(c57$unimodal, cast$unimodal))  # unimodal
} else if (nn == 4){
  group1_name <- "bimodal"
  group2_name <- "unimodal"
  
  h3k4me3_gr <- import(
    "peak/ENCFF671UNN.bed",
    format = "BED",
    extraCols = extraCols_narrowPeak
  )
  sheet_c57 <- read_excel("gene name.xlsx", sheet = "E c57")
  sheet_cast <- read_excel("gene name.xlsx", sheet = "E cast")
  extract_genes <- function(sheet) {
    bimodal <- c(sheet[[1]][!is.na(sheet[[1]])], 
                 sheet[[3]][!is.na(sheet[[3]])])
    unimodal <- setdiff(sheet[[5]][!is.na(sheet[[5]])],bimodal)
    return(list(bimodal = unique(bimodal), 
                unimodal = unique(unimodal)))
  }
  c57 <- extract_genes(sheet_c57)
  cast <- extract_genes(sheet_cast)
  my_genes <- unique(c(c57$bimodal, cast$bimodal))     # bimodal
  my_genes1 <- unique(c(c57$unimodal, cast$unimodal))  # unimodal
}
rm(sheet_c57, sheet_cast, extract_genes, c57, cast )
all_genes <- unique(c(my_genes, my_genes1))

gene_group_df <- bind_rows(
  data.frame(
    SYMBOL = my_genes,
    group = group1_name,
    stringsAsFactors = FALSE
  ),
  data.frame(
    SYMBOL = my_genes1,
    group = group2_name,
    stringsAsFactors = FALSE
  )
) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  mutate(
    group = factor(group, levels = c(group1_name, group2_name))
  )


seqlevelsStyle(h3k4me3_gr) <- "UCSC"
h3k4me3_gr <- keepStandardChromosomes(
  h3k4me3_gr,
  pruning.mode = "coarse"
)
merge_adjacent_peaks <- TRUE
merge_gap <- 500

if (merge_adjacent_peaks) {
  h3k4me3_region_gr <- reduce(
    h3k4me3_gr,
    min.gapwidth = merge_gap + 1,
    ignore.strand = TRUE
  )
  
  mcols(h3k4me3_region_gr)$signalValue <- NA_real_
  
} else {
  h3k4me3_region_gr <- h3k4me3_gr
}

h3k4me3_region_gr <- keepStandardChromosomes(
  h3k4me3_region_gr,
  pruning.mode = "coarse"
)

gencode_gtf <- import("gencode.vM25.annotation.gtf.gz")
seqlevelsStyle(gencode_gtf) <- "UCSC"
tx_gr <- gencode_gtf[gencode_gtf$type == "transcript"]
tx_gr <- keepStandardChromosomes(
  tx_gr,
  pruning.mode = "coarse"
)
tx_df <- data.frame(
  gene_id = tx_gr$gene_id,
  transcript_id = tx_gr$transcript_id,
  SYMBOL = tx_gr$gene_name,
  gene_type = tx_gr$gene_type,
  chr = as.character(seqnames(tx_gr)),
  start = start(tx_gr),
  end = end(tx_gr),
  strand = as.character(strand(tx_gr)),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(SYMBOL)) %>%
  mutate(
    tss = ifelse(strand == "+", start, end)
  )

tss_df <- tx_df %>%
  left_join(
    gene_group_df,
    by = "SYMBOL"
  ) %>%
  filter(!is.na(group)) %>%
  distinct(SYMBOL, chr, tss, strand, group, .keep_all = TRUE) %>%
  mutate(
    tss_id = row_number()
  )

cat("TSS used for analysis:", nrow(tss_df), "\n")
cat("Genes with annotated TSS:", length(unique(tss_df$SYMBOL)), "\n")

missing_genes <- setdiff(gene_group_df$SYMBOL, unique(tss_df$SYMBOL))
tss_gr <- GRanges(
  seqnames = tss_df$chr,
  ranges = IRanges(start = tss_df$tss, width = 1),
  strand = tss_df$strand,
  tss_id = tss_df$tss_id,
  gene_symbol = tss_df$SYMBOL,
  group = tss_df$group
)

seqlevelsStyle(tss_gr) <- "UCSC"

tss_gr <- keepStandardChromosomes(
  tss_gr,
  pruning.mode = "coarse"
)

overlaps_tss <- findOverlaps(
  tss_gr,
  h3k4me3_region_gr,
  ignore.strand = TRUE
)

if (length(overlaps_tss) > 0) {
  
  tss_indices <- queryHits(overlaps_tss)
  peak_indices <- subjectHits(overlaps_tss)
  
  overlap_df <- data.frame(
    tss_id = tss_gr$tss_id[tss_indices],
    SYMBOL = tss_gr$gene_symbol[tss_indices],
    group = tss_gr$group[tss_indices],
    chr = as.character(seqnames(tss_gr[tss_indices])),
    strand = as.character(strand(tss_gr[tss_indices])),
    tss = start(tss_gr[tss_indices]),
    h3k4me3_start = start(h3k4me3_region_gr[peak_indices]),
    h3k4me3_end = end(h3k4me3_region_gr[peak_indices]),
    h3k4me3_width = width(h3k4me3_region_gr[peak_indices]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      left_length = tss - h3k4me3_start,
      right_length = h3k4me3_end - tss,
      upstream = ifelse(strand == "+", left_length, right_length),
      downstream = ifelse(strand == "+", right_length, left_length)
    )
  max_width_per_tss <- overlap_df %>%
    group_by(tss_id) %>%
    slice_max(
      order_by = h3k4me3_width,
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  tss_result <- tss_df %>%
    dplyr::select(tss_id, SYMBOL, group, chr, tss, strand) %>%
    left_join(
      max_width_per_tss %>%
        dplyr::select(
          tss_id,
          h3k4me3_start,
          h3k4me3_end,
          h3k4me3_width,
          upstream,
          downstream
        ),
      by = "tss_id"
    ) %>%
    mutate(
      h3k4me3_width = ifelse(is.na(h3k4me3_width), 0, h3k4me3_width),
      upstream = ifelse(is.na(upstream), 0, upstream),
      downstream = ifelse(is.na(downstream), 0, downstream),
      has_h3k4me3 = h3k4me3_width > 0
    )
  
} else {
  
  tss_result <- tss_df %>%
    dplyr::select(tss_id, SYMBOL, group, chr, tss, strand) %>%
    mutate(
      h3k4me3_start = NA_integer_,
      h3k4me3_end = NA_integer_,
      h3k4me3_width = 0,
      upstream = 0,
      downstream = 0,
      has_h3k4me3 = FALSE
    )
}
result <- tss_result %>%
  group_by(SYMBOL, group) %>%
  arrange(desc(h3k4me3_width), .by_group = TRUE) %>%
  summarise(
    n_tss = n(),
    selected_chr = dplyr::first(chr),
    selected_tss = dplyr::first(tss),
    selected_strand = dplyr::first(strand),
    h3k4me3_start = dplyr::first(h3k4me3_start),
    h3k4me3_end = dplyr::first(h3k4me3_end),
    h3k4me3_width = dplyr::first(h3k4me3_width),
    upstream = dplyr::first(upstream),
    downstream = dplyr::first(downstream),
    has_h3k4me3 = dplyr::first(has_h3k4me3),
    .groups = "drop"
  ) %>%
  mutate(
    group = factor(group, levels = c("cross-talk", "two-state")),
    has_h3k4me3 = as.logical(has_h3k4me3)
  )

wilcox_all <- wilcox.test(
  h3k4me3_width ~ group,
  data = result,
  exact = FALSE
)

cat("\n=== Wilcoxon test: all genes including zero values ===\n")
cat("p-value:", wilcox_all$p.value, "\n")
if (nn == 1 | nn == 2) {
  color_values <- c("cross-talk" = "#E41A1C", "two-state" = "#377EB8")
} else if (nn == 3 | nn == 4) {
  color_values <- c("bimodal" = "#E41A1C", "unimodal" = "#377EB8")
}


p_all <- ggplot(
  result,
  aes(x = group, y = h3k4me3_width, fill = group)
) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.8,
    outlier.shape = NA
  ) +
  scale_fill_manual(values = color_values) +
  coord_cartesian(
    ylim = c(0, 5000)
  ) +
  labs(
    title = "Overall H3K4me3 Modification Difference",
    subtitle = paste0(
      "All genes included; genes without TSS-overlapping H3K4me3 are assigned width = 0\n",
      "Wilcoxon p = ", signif(wilcox_all$p.value, 3)
    ),
    x = "Group",
    y = "H3K4me3 width (bp)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )
print(p_all)










