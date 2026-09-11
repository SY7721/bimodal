
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(IRanges)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)
library(stringr)

nn <- 4 # 1, 2, 3, 4

if (nn == 1){
  group1_name <- "cross-talk"
  group2_name <- "two-state"
  yy = c(3, 10)
  mef_loop <- read.delim(
    gzfile("loop/GSE113339_MEF-H3K27ac_high-confidence_interactions.tsv.gz"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
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
  yy = c(3, 11)
  mef_loop <- read.delim(
    gzfile("loop/GSE113339_ES-H3K27ac_high-confidence_interactions.tsv.gz"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
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
  yy = c(3, 10)
  mef_loop <- read.delim(
    gzfile("loop/GSE113339_MEF-H3K27ac_high-confidence_interactions.tsv.gz"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
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
  yy = c(3, 11)
  mef_loop <- read.delim(
    gzfile("loop/GSE113339_ES-H3K27ac_high-confidence_interactions.tsv.gz"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
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


gene_group <- data.frame(
  gene = all_genes,
  group = ifelse(
    all_genes %in% my_genes,
    group1_name,
    group2_name
  )
)

mef_loop$loop_strength <- rowMeans(
  mef_loop[, c("cpm.rep1", "cpm.rep2")],
  na.rm = TRUE
)
summary(mef_loop$loop_strength)
parse_bin <- function(x) {x_match <- str_match(x,
                                               "^(chr[^:.]+)[:.]([0-9]+)[-.]([0-9]+)$" )
data.frame(chr = x_match[, 2],
           start = as.numeric(x_match[, 3]),
           end = as.numeric(x_match[, 4]),
           stringsAsFactors = FALSE)
}
anchor1 <- parse_bin(
  mef_loop$bin.start
)
anchor2 <- parse_bin(
  mef_loop$bin.end
)
mef_loop$chr1 <- anchor1$chr
mef_loop$start1 <- anchor1$start
mef_loop$end1 <- anchor1$end
mef_loop$chr2 <- anchor2$chr
mef_loop$start2 <- anchor2$start
mef_loop$end2 <- anchor2$end
rm(anchor1, anchor2)
mef_loop$loop_id <- 1:nrow(mef_loop)
mef_loop <- mef_loop %>%
  filter(chr1 == chr2)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr <- genes(txdb)
gene_names <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = names(genes_gr),
  columns = "SYMBOL",
  keytype = "ENTREZID"
)
genes_df <- data.frame(
  gene_id = names(genes_gr),
  chr = as.character(seqnames(genes_gr)),
  start = start(genes_gr),
  end = end(genes_gr),
  strand = as.character(strand(genes_gr)),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    gene_names,
    by = c("gene_id" = "ENTREZID")
  ) %>%
  dplyr::filter(
    !is.na(SYMBOL)
  ) %>%
  dplyr::mutate(
    tss = ifelse(
      strand == "+",
      start,
      end
    ),
    promoter_start = pmax(
      1,
      tss - 1000
    ),
    promoter_end = tss + 1000
  )
all_promoters_gr <- GRanges(
  seqnames = genes_df$chr,
  ranges = IRanges(
    start = genes_df$promoter_start,
    end = genes_df$promoter_end
  ),
  gene = genes_df$SYMBOL
)
target_promoters_gr <- all_promoters_gr[
  all_promoters_gr$gene %in% all_genes
]

anchor1_gr <- GRanges(
  seqnames = mef_loop$chr1,
  ranges = IRanges(
    start = mef_loop$start1,
    end = mef_loop$end1
  ),
  loop_id = mef_loop$loop_id
)
anchor2_gr <- GRanges(
  seqnames = mef_loop$chr2,
  ranges = IRanges(
    start = mef_loop$start2,
    end = mef_loop$end2
  ),
  loop_id = mef_loop$loop_id
)

anchor1_in_promoter <- countOverlaps(
  anchor1_gr,
  all_promoters_gr,
  ignore.strand = TRUE
) > 0
anchor2_in_promoter <- countOverlaps(
  anchor2_gr,
  all_promoters_gr,
  ignore.strand = TRUE
) > 0

idx1 <- which(
  anchor1_in_promoter &
    !anchor2_in_promoter
)
overlap1 <- findOverlaps(
  anchor1_gr[idx1],
  target_promoters_gr,
  ignore.strand = TRUE
)
EP1 <- data.frame(loop_id = anchor1_gr$loop_id[idx1][queryHits(overlap1)],
                  gene = target_promoters_gr$gene[subjectHits(overlap1)],
                  promoter_anchor = "anchor1",
                  stringsAsFactors = FALSE)
idx2 <- which(
  !anchor1_in_promoter &
    anchor2_in_promoter
)
overlap2 <- findOverlaps(
  anchor2_gr[idx2],
  target_promoters_gr,
  ignore.strand = TRUE
)
EP2 <- data.frame( loop_id = anchor2_gr$loop_id[idx2][ queryHits(overlap2)],
                   gene = target_promoters_gr$gene[subjectHits(overlap2)],
                   promoter_anchor = "anchor2",
                   stringsAsFactors = FALSE)

EP_gene <- bind_rows( EP1,EP2) %>% distinct()
EP_gene <- EP_gene %>%
  left_join(
    mef_loop %>%
      dplyr::select(
        loop_id,
        loop_strength
      ),
    by = "loop_id"
  )
sum(is.na(EP_gene$loop_strength))

gene_EP_summary <- EP_gene %>% dplyr::group_by(gene) %>% dplyr::summarise(
  n_EP_loops = dplyr::n_distinct(loop_id),
  EP_contact_strength = sum( loop_strength,  na.rm = TRUE ),
  mean_loop_strength = mean( loop_strength,  na.rm = TRUE),
  median_loop_strength = median(loop_strength, na.rm = TRUE ),
  max_loop_strength = max( loop_strength, na.rm = TRUE ),
  .groups = "drop")
head(gene_EP_summary)
result <- gene_group %>% dplyr::left_join( gene_EP_summary, by = "gene" ) %>%
  dplyr::mutate(n_EP_loops = tidyr::replace_na(n_EP_loops, 0 ),
                EP_contact_strength = tidyr::replace_na( EP_contact_strength, 0 ),
                has_EP = n_EP_loops > 0 )
EP_presence_summary <- result %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n_genes = dplyr::n(),
                   n_EP_positive = sum(has_EP),
                   n_EP_negative = sum(!has_EP),
                   EP_positive_percent = mean(has_EP) * 100,
                   .groups = "drop")
EP_table <- table(result$group,result$has_EP)
fisher_EP <- fisher.test( EP_table)
result_EP <- result %>% dplyr::filter( has_EP)

EP_positive_summary <- result_EP %>%
  dplyr::group_by(group) %>%
  dplyr::summarise( n_genes = dplyr::n(),
   mean_EP_loops = mean( n_EP_loops ),
   median_EP_loops = median( n_EP_loops),
    Q1_EP_loops = quantile( n_EP_loops,0.25 ),
   Q3_EP_loops = quantile( n_EP_loops, 0.75 ),
   mean_loop_CPM = mean( mean_loop_strength, na.rm = TRUE ),
   median_loop_CPM = median( mean_loop_strength, na.rm = TRUE),
   mean_total_EP_strength = mean(EP_contact_strength),
   median_total_EP_strength = median( EP_contact_strength ),
   .groups = "drop"
  )

wilcox_mean_strength <- wilcox.test(
  mean_loop_strength ~ group,
  data = result_EP,
  exact = FALSE
)
wilcox_mean_strength
if (nn == 1 | nn == 2) {
  color_values <- c("cross-talk" = "#D55E00", "two-state" = "#0072B2")
} else if (nn == 3 | nn == 4) {
  color_values <- c("bimodal" = "#D55E00", "unimodal" = "#0072B2")
}

p_mean_strength <- ggplot(
  result_EP,
  aes(
    x = group,
    y = mean_loop_strength,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(color_values) ) +
  labs(
    title = "Mean H3K27ac-associated E-P Contact Strength",
    subtitle = paste0(
      "Wilcoxon p = ",
      signif(
        wilcox_mean_strength$p.value,
        3
      )
    ),
    x = "",
    y = " E-P intensity"
  ) +
  coord_cartesian(
    ylim = yy
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    axis.text = element_text(
      color = "black"
    ),
    
    panel.grid.minor = element_blank()
  )

p_mean_strength

