library(DESeq2)
library(readxl)
library(pROC)
library(philentropy)

counts_data <- read_excel("tw_th_roc.xlsx", sheet = "tw_th")                      
counts_data <- as.data.frame(counts_data)
rownames(counts_data) <- counts_data[, 1]  
counts_data <- counts_data[, -1]                                  
counts_data[is.na(counts_data)] <- 0
counts_data1 <- read_excel("tw_th_roc.xlsx", sheet = "tw")                     
counts_data1 <- as.data.frame(counts_data1)
rownames(counts_data1) <- counts_data1[, 1] 
counts_data1 <- counts_data1[, -1]                              
labels <- c(rep(0, 500), rep(1, 500))

mat1 <- as.matrix(counts_data)
mat2 <- as.matrix(counts_data1)
distribution <- function(x,N){
  max_val <- max(x)
  y <- numeric(max_val+1)
  for (i in 1:(max_val+1)){
    y[i] <- sum(x == (i-1))/N
  }
  return(y)
}

diver <- function(g1, g2){
  N1 <- ncol(mat1)
  p <- distribution(g1, N1)
  q <- distribution(g2, N1)
  max_len <- max(length(p), length(q))
  length(p) <- max_len
  length(q) <- max_len
  p[is.na(p)] <- 0
  q[is.na(q)] <- 0
  p <- p+ 1e-10
  q <- q+ 1e-10
  p <- p / sum(p)
  q <- q / sum(q)
  k1_pq <- KL(rbind(p,q), unit = "log2")
  k1_qp <- KL(rbind(q,p), unit = "log2")
  sym_kl <- (k1_pq+k1_qp)/2
  js <- JSD(rbind(p,q), unit = "log2")
  
  return(c(k1_pq = k1_pq, k1_qp = k1_qp, sym_kl = sym_kl, JS = js))
}
result <- t(sapply(1:nrow(mat1), function(i) {
  diver(mat1[i, ], mat2[i, ])
}))

score <- result[1:nrow(mat1),4]

roc_obj <- roc(labels, score, 
               levels = c(0, 1), 
               direction = "<") 
print(roc_obj)

auc_value <- auc(roc_obj)
print(paste("AUC:", round(auc_value, 4)))

plot(roc_obj, 
     main = paste("ROC Curve (AUC =", round(auc_value, 4), ")"),
     col = "blue", 
     lwd = 2,
     legacy.axes = TRUE,         
     xlab = "False Positive Rate",  
     ylab = "True Positive Rate")       

