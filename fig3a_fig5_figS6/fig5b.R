

library(readxl)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(stringr)
library(dplyr)
library(tidyr)

nn <- 1 # 1, 2, 3, 4
if (nn == 1){
  group1_name <- "cross-talk"
  group2_name <- "two-state"
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
  my_genes1 <- unique(c(c57$telegraph, cast$telegraph))  # telegraph
} else if (nn == 3){
  group1_name <- "bimodal"
  group2_name <- "unimodal"
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

promoter_feature_c <-  data.frame(gene = my_genes)
promoter_feature_t <-  data.frame(gene = my_genes1)


promofile = list.files('EPD data',pattern="mouse.*.bed")
for (i in 1:length(promofile)){
  data <- read.table(paste0("EPD data/",promofile[i]),sep ='\t',header = FALSE)
  data$V4 = str_replace(data$V4, "_1","")
  data = data[c('V4',"V5")]
  promoter_feature_c = merge(promoter_feature_c,data,by.x = "gene",by.y = "V4",all.x = TRUE) 
}
for (i in 1:length(promofile)){
  data <- read.table(paste0("EPD data/",promofile[i]),sep ='\t',header = FALSE)
  data$V4 = str_replace(data$V4, "_1","")
  data = data[c('V4',"V5")]
  promoter_feature_t = merge(promoter_feature_t,data,by.x = "gene",by.y = "V4",all.x = TRUE) 
}

colnames(promoter_feature_c) <- c('genename','TATA','Inr','CCAAT','GC')
colnames(promoter_feature_t) <- c('genename','TATA','Inr','CCAAT','GC')

promoter_feature_c[is.na(promoter_feature_c)] <- 0
promoter_feature_t[is.na(promoter_feature_t)] <- 0
feature_name <- c('TATA', 'Inr', 'CCAAT', 'GC')

combined_df <- rbind(
  data.frame(promoter_feature_c[, c('TATA', 'Inr', 'CCAAT', 'GC')], Group = group1_name),  
  data.frame(promoter_feature_t[, c('TATA', 'Inr', 'CCAAT', 'GC')], Group = group2_name)  
)

prop_results <- combined_df %>%
  group_by(Group) %>%
  summarise(
    TATA_prop = mean(TATA == 1, na.rm = TRUE),
    Inr_prop = mean(Inr == 1, na.rm = TRUE),
    CCAAT_prop = mean(CCAAT == 1, na.rm = TRUE),
    GC_prop = mean(GC == 1, na.rm = TRUE)
  ) %>%
  mutate(across(ends_with("_prop"), ~ .x * 100)) 

prop_long <- prop_results %>%
  pivot_longer(cols = -Group, 
               names_to = "Feature", 
               values_to = "Percentage") %>%
  mutate(
    Feature = str_remove(Feature, "_prop"),
    Feature = factor(Feature, levels = c("TATA", "Inr", "GC", "CCAAT"))
  )


# 修改绘图部分 - 根据实际分组动态设置颜色
# 获取实际的分组名称
actual_groups <- unique(prop_long$Group)

# 根据分组名称设置颜色
if (nn == 1 | nn == 2) {
  color_values <- c("cross-talk" = "#4595D1", "two-state" = "#EB9A18")
} else if (nn == 3 | nn == 4) {
  color_values <- c("bimodal" = "#4582D1", "unimodal" = "#EB9A98")
}

# 只保留实际存在的分组
color_values <- color_values[names(color_values) %in% actual_groups]

ggplot(prop_long, aes(x = Feature, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_dodge(width = 0.8),
            vjust = -0.3,
            size = 3.5) +
  scale_fill_manual(values = color_values) +  # 使用动态颜色映射
  labs(title = "Promoter Element Enrichment",
       x = "Promoter Element",
       y = "Percentage of genes with promoter elements (%)") +
  coord_cartesian(ylim = c(0, 60)) + 
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))
