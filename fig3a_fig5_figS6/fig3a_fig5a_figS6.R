
#BiocManager::install("readxl")
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

df_BP <- read_excel("only_E.xlsx", 
                    sheet = "BP")  
df_CC <- read_excel("only_E.xlsx", 
                    sheet = "CC")
df_MF <- read_excel("only_E.xlsx", 
                    sheet = "MF")  

df_BP <- df_BP %>% mutate(Ontology = "BP")
df_CC <- df_CC %>% mutate(Ontology = "CC")
df_MF <- df_MF %>% mutate(Ontology = "MF")

df_BP$Benjamini <- as.numeric(df_BP$Benjamini)
df_CC$Benjamini <- as.numeric(df_CC$Benjamini)
df_MF$Benjamini <- as.numeric(df_MF$Benjamini)
df_BP$FDR <- as.numeric(df_BP$FDR)
df_CC$FDR <- as.numeric(df_CC$FDR)
df_MF$FDR <- as.numeric(df_MF$FDR)
df_BP$Bonferroni <- as.numeric(df_BP$Bonferroni)
df_CC$Bonferroni <- as.numeric(df_CC$Bonferroni)
df_MF$Bonferroni <- as.numeric(df_MF$Bonferroni)
go_df <- bind_rows(df_BP, df_CC, df_MF)


names(go_df)
go_df_clean <- go_df %>%
  mutate(
    PValue = as.numeric(`P-Value`),
    Benjamini = as.numeric(Benjamini),
    Count = as.numeric(Count),
    FoldEnrichment = as.numeric(`Fold Enrichment`),
    logP = -log10(PValue),
    Term_clean = ifelse(grepl("~", Term), 
                        sub(".*~", "", Term), 
                        Term),
    Term_clean = str_wrap(Term_clean, width = 35)
  ) %>%
  filter(!is.na(PValue), !is.na(Benjamini), !is.na(Count))

go_df_filtered <- go_df_clean %>%
  filter(Benjamini < 0.05) %>%
  group_by(Ontology) %>%
  arrange(Benjamini) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  arrange(Ontology, desc(logP)) %>%
  mutate(Term_clean = factor(Term_clean, levels = unique(Term_clean)))

print(table(go_df_filtered$Ontology))
ontology_colors <- c("BP" = "#E64B35",  
                     "CC" = "#4DBBD5", 
                     "MF" = "#00A087")  

print(go_df_filtered %>% 
        group_by(Ontology) %>% 
        summarise(
          number = n(),
          logP = paste(round(min(logP),1), "~", round(max(logP),1)),
          min_Benjamini = format(min(Benjamini), scientific = TRUE)
        ))

p_no_labels <- ggplot(go_df_filtered, aes(x = logP, y = Term_clean, fill = Ontology)) +
 geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", logP)), 
            hjust = -0.2, size = 3, color = "black") +
  facet_grid(factor(Ontology, levels = c("BP", "CC", "MF")) ~ ., 
             scales = "free_y", 
             space = "free_y") +
  scale_fill_manual(values = ontology_colors, 
                    name = NULL,  
                    labels = c("BP" = "BP",
                               "CC" = "CC",
                               "MF" = "MF")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-Log10(P value)", y = NULL,
       title = "cast_only_E") +
  theme_minimal() +
  theme(
    strip.text = element_blank(), 
    strip.background = element_blank(),  
    panel.spacing = unit(0, "lines"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8),
    legend.position = c(0.88, 0.08), 
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", color = "grey70", 
                                     linewidth = 0.3),
    legend.key.size = unit(0.4, "cm"),  
    legend.text = element_text(size = 7), 
    legend.margin = margin(2, 2, 2, 2) 
  )

print(p_no_labels)


### KEGG
kegg_df <- read_excel("only_E.xlsx", 
                      sheet = "KEGG") 
names(kegg_df)
kegg_clean <- kegg_df %>%
  mutate(
    PValue = as.numeric(`P-Value`),
    Benjamini = as.numeric(Benjamini),
    Count = as.numeric(Count),
    FoldEnrichment = as.numeric(`Fold Enrichment`),
    logP = -log10(PValue),
    Pathway = ifelse(grepl(":", Term), 
                     sub(".*:", "", Term), 
                     Term),
    Pathway = str_wrap(Pathway, width = 35)
  ) %>%
  filter(!is.na(PValue), !is.na(Benjamini), !is.na(Count)) %>%
  filter(Benjamini < 0.05) %>%
  arrange(Benjamini) %>%
  slice_head(n = 20) %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))

print(kegg_clean[, c("Pathway", "logP", "Count", "FoldEnrichment", "Benjamini")])
p <- ggplot(kegg_clean, aes(x = FoldEnrichment, y = Pathway)) +
  geom_point(aes(size = Count, color = logP), alpha = 0.8) +
  scale_color_gradient(low = "#4DBBD5",   
                       high = "#E64B35", 
                       name = "-Log10(P value)") +
  scale_size_continuous(name = "Gene Count", 
                        range = c(3, 12)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  labs(x = "Fold Enrichment", y = NULL,
       title = "cast_only_E") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.6, "cm"),
    legend.box = "vertical"
  )

print(p)

