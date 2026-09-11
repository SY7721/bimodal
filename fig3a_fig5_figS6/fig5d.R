library(readxl)
library(dplyr)
library(ggplot2)

nn <- 4 # 1, 2, 3, 4
if (nn == 1){
  load("promoter width/promoter_width_MEF.Rdata")
  promoter_width <- promoter_width_MEF
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
  y1 <- c(0,75)
} else if (nn == 2){
  load("promoter width/promoter_width_ESC.Rdata")
  promoter_width <- promoter_width_ESC
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
  y1 <- c(0,80)
} else if (nn == 3){
  load("promoter width/promoter_width_MEF.Rdata")
  promoter_width <- promoter_width_MEF
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
  y1 <- c(0,75)
} else if (nn == 4){
  load("promoter width/promoter_width_ESC.Rdata")
  promoter_width <- promoter_width_ESC
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
  y1 <- c(0,75)
}
rm(sheet_c57, sheet_cast, extract_genes, c57, cast )

group1_width <- promoter_width %>%
  filter(
    gene_name %in% my_genes
  ) %>%
  mutate(
    group = group1_name
  )
group2_width <- promoter_width %>%
  filter(
    gene_name %in% my_genes1
  ) %>%
  mutate(
    group = group2_name
  )
promoter_width_plot <- bind_rows(
  group1_width,
  group2_width
)

wilcox_all <- wilcox.test(
  width ~ group,
  data = promoter_width_plot
)
p_label <- paste0(
  "Wilcoxon p = ",
  signif(
    wilcox_all$p.value,
    3
  )
)
if (nn == 1 | nn == 2) {
  color_values <- c("cross-talk" = "#D55E00", "two-state" = "#009E73")
} else if (nn == 3 | nn == 4) {
  color_values <- c("bimodal" = "#E41A1C", "unimodal" = "#377EB8")
}

p <- ggplot(promoter_width_plot,aes(x = group, y = width, fill = group))+
     geom_boxplot(width = 0.6,alpha = 0.85,outlier.shape = NA,color = "black")+
     scale_fill_manual( values = color_values)+
     labs( x = NULL,y = "Promoter width", subtitle = p_label)+
     coord_cartesian( ylim = y1)+
     theme_bw()+ theme(legend.position = "none",
         panel.grid.major = element_line(   color = "grey85",   linewidth = 0.3 ),
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(   size = 12,   color = "black" ),
         axis.text.y = element_text(   size = 12,   color = "black" ),
         axis.title.y = element_text(   size = 14),
         plot.subtitle = element_text(   size = 12,   hjust = 0.5 ))

p
