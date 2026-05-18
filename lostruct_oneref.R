#create personal library to write packages to
#dir.create("~/R/x86_64-pc-linux-gnu-library/4.4", recursive = TRUE, showWarnings = FALSE)
#set library paths
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")

library(data.table)
#devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)

# read the cpntest output (loci x individuals)
coded <- as.matrix(read.table("cpntest_FHA_all.txt"))

#make positions file 
positions <- read.table("positions.txt", col.names = c("chrom", "pos"))

#run pca windows
eigenstuff <- eigen_windows(coded, win=100, k=2)

#comput pairwise distances
windist <- pc_dist( eigenstuff, npc=2 )

#mds on window distances
mds_result <- cmdscale( windist, eig=TRUE, k=2 )
mds_df <- as.data.frame(mds_result$points)
colnames(mds_df) <- c("MDS1", "MDS2")

# check MDS distribution first
summary(mds_df$MDS1)
summary(mds_df$MDS2)

# check how many windows at different thresholds
cat("Outliers at 4 SD:", sum(abs(scale(mds_df$MDS1)) > 4 | abs(scale(mds_df$MDS2)) > 4), "\n")
cat("Outliers at 3 SD:", sum(abs(scale(mds_df$MDS1)) > 3 | abs(scale(mds_df$MDS2)) > 3), "\n")
cat("Outliers at 2 SD:", sum(abs(scale(mds_df$MDS1)) > 2 | abs(scale(mds_df$MDS2)) > 2), "\n")

# add window coordinates
mds_df$window <- 1:nrow(mds_df)

#identify outlier windows

sd_thresh<- 3
mds_df$outlier <- abs(scale(mds_df$MDS1)) > sd_thresh | 
                  abs(scale(mds_df$MDS2)) > sd_thresh

cat("Number of outlier windows:", sum(mds_df$outlier), "\n")

#plot
plot(mds_df$MDS1, mds_df$MDS2,
     col = ifelse(mds_df$outlier, "red", "grey50"),
     pch = 19,
     xlab = "MDS1", ylab = "MDS2",
     main = "Window PCA distances - outliers in red")

#extract SNPs for outlier windows and run PCA

outlier_wins <- which(mds_df$outlier)

outlier_table <- data.frame(
  window     = outlier_wins,
  snp_start  = (outlier_wins - 1) * 100 + 1,
  snp_end    = pmin(outlier_wins * 100, nrow(coded)),
  chrom_start = positions$chrom[(outlier_wins - 1) * 100 + 1],
  pos_start   = positions$pos[(outlier_wins - 1) * 100 + 1],
  chrom_end   = positions$chrom[pmin(outlier_wins * 100, nrow(coded))],
  pos_end     = positions$pos[pmin(outlier_wins * 100, nrow(coded))],
  MDS1       = mds_df$MDS1[outlier_wins],
  MDS2       = mds_df$MDS2[outlier_wins],
  MDS1_zscore = scale(mds_df$MDS1)[outlier_wins],
  MDS2_zscore = scale(mds_df$MDS2)[outlier_wins]
)

write.table(outlier_table, 
            file = "outlier_windows.txt", 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE)

cat("Outlier window table written to outlier_windows.txt\n")
print(outlier_table)


for (w in outlier_wins) {
  # get SNP indices for this window
  snp_start <- (w - 1) * 100 + 1
  snp_end   <- min(w * 100, nrow(coded))
  
  win_snps <- coded[snp_start:snp_end, ]
  
  # PCA on individuals for this window
  pca_out <- prcomp(t(win_snps), scale. = FALSE)
  
  # k-means clustering
  km <- kmeans(pca_out$x[, 1:2], centers = 3, nstart = 25)
  
  # plot
  plot(pca_out$x[, 1], pca_out$x[, 2],
       col = km$cluster,
       pch = 19,
       xlab = "PC1", ylab = "PC2",
       main = paste("Window", w, "|", 
                    positions$chrom[snp_start], ":", positions$pos[snp_start], "-",positions$pos[snp_end]))
  legend("topright", legend = paste("cluster", 1:3), 
         col = 1:3, pch = 19)
}

