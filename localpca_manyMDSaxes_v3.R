library(data.table)
library(lostruct)
library(ggplot2)
library(tidyr)
library(cluster)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript localpca_manyMDSaxes.R <input_file> <output_prefix> <window_size_snps> <n_axes>")
}
input_file    <- args[1]
output_prefix <- args[2]
win_size      <- as.numeric(args[3])  ## number of SNPs per window
n_axes        <- as.numeric(args[4])  ## number of MDS axes to compute

coded <- as.matrix(read.table(input_file))

# load positions
positions <- read.table("positions.txt", header = FALSE,
                        col.names = c("chrom", "pos"))
positions$pos <- as.numeric(positions$pos)

# filter zero-variance rows and keep positions in sync
row_vars <- apply(coded, 1, var)
keep      <- row_vars > 1e-10
coded     <- coded[keep, ]
positions <- positions[keep, ]

# standardize
g_scaled <- t(scale(t(coded)))

cat("Loaded:", nrow(coded), "SNPs x", ncol(coded), "individuals\n")
cat("Window size:", win_size, "SNPs\n")

# run lostruct eigen_windows (SNP windows)
eigenstuff <- eigen_windows(g_scaled, win = win_size, k = 2)

# track windows that were dropped (returned NA)
na_wins    <- which(is.na(eigenstuff[, 1]))
valid_wins <- setdiff(1:nrow(eigenstuff), na_wins)
cat("Windows retained:", length(valid_wins), "| Dropped:", length(na_wins), "\n")

eigenstuff_valid <- eigenstuff[valid_wins, , drop = FALSE]

# PC distance between axes 1 & 2 between windows
windist <- pc_dist(eigenstuff_valid, npc = 2)

# cap n_axes to max possible
n_axes <- min(n_axes, nrow(windist) - 1)
cat("Running MDS with", n_axes, "axes\n")

fit2d      <- cmdscale(windist, eig = TRUE, k = n_axes)
mds_points <- fit2d$points
colnames(mds_points) <- paste0("MDS", 1:n_axes)

# get genomic coordinates for each valid window
n_snps <- nrow(g_scaled)
n_wins <- floor(n_snps / win_size)

all_win_coords <- data.frame(
  win_index = 1:n_wins,
  chrom     = positions$chrom[ seq(1, n_wins * win_size, by = win_size) ],
  start_pos = positions$pos[   seq(1, n_wins * win_size, by = win_size) ],
  end_pos   = positions$pos[   seq(win_size, n_wins * win_size, by = win_size) ],
  mid_pos   = rowMeans(cbind(
    positions$pos[ seq(1,        n_wins * win_size, by = win_size) ],
    positions$pos[ seq(win_size, n_wins * win_size, by = win_size) ]
  ))
)

win_coords <- all_win_coords[valid_wins, ]

mds_df <- data.frame(win_index = valid_wins, mds_points)
mds_df <- merge(mds_df, win_coords, by = "win_index")
mds_df <- mds_df[order(mds_df$chrom, mds_df$start_pos), ]
mds_df$chrom_label <- sub("^(Scaffold_[^_]+)_.*", "\\1", mds_df$chrom)

### select best k for genome-wide window clustering ###
find_best_k <- function(data, max_k = 6, seed = 42) {
  set.seed(seed)
  max_k <- min(max_k, nrow(data) - 1)
  if (max_k < 2) {
    cat("Not enough windows to test k > 1, defaulting to k=1\n")
    return(1)
  }

  wss <- numeric(max_k)
  sil <- numeric(max_k)
  wss[1] <- (nrow(data) - 1) * sum(apply(data, 2, var))
  sil[1] <- NA

  for (k in 2:max_k) {
    km       <- kmeans(data, centers = k, nstart = 25)
    wss[k]   <- km$tot.withinss
    sil_vals <- silhouette(km$cluster, dist(data))
    sil[k]   <- mean(sil_vals[, 3])
  }

  scores_df <- data.frame(k = 1:max_k, wss = wss, mean_silhouette = sil)
  cat("\nK-means selection scores:\n")
  print(scores_df)

  write.table(scores_df,
              file      = paste0(output_prefix, "_", win_size, "snp_kmeans_scores.txt"),
              sep       = "\t",
              row.names = FALSE,
              quote     = FALSE)

  best_k <- which.max(sil[2:max_k]) + 1
  cat("Best k by silhouette:", best_k, "(score =", round(sil[best_k], 4), ")\n\n")

  p_elbow <- ggplot(scores_df[-1, ], aes(x = k, y = wss)) +
    geom_line() + geom_point(size = 3) +
    theme_classic() +
    labs(title = "K-means elbow plot (genome-wide windows)",
         x = "k", y = "Within-cluster SS")
  ggsave(paste0(output_prefix, "_", win_size, "snp_kmeans_elbow.png"), p_elbow)

  p_sil <- ggplot(scores_df[!is.na(scores_df$mean_silhouette), ],
                  aes(x = k, y = mean_silhouette)) +
    geom_line() + geom_point(size = 3) +
    geom_point(data = scores_df[best_k, ], color = "firebrick", size = 4) +
    theme_classic() +
    labs(title = "K-means silhouette scores (genome-wide windows)",
         x = "k", y = "Mean silhouette score")
  ggsave(paste0(output_prefix, "_", win_size, "snp_kmeans_silhouette.png"), p_sil)

  return(best_k)
}

best_k         <- find_best_k(mds_points, max_k = 6)
set.seed(42)
km_best        <- kmeans(mds_points, centers = best_k, nstart = 25)
mds_df$cluster <- km_best$cluster

# Outlier flagging 
flag_outliers <- function(x, thresh = 3) {
  abs(x - mean(x, na.rm = TRUE)) > thresh * sd(x, na.rm = TRUE)
}

# Find consecutive runs of outlier windows on a single chromosome 
# Returns a data.frame of qualifying runs (n_windows >= min_windows)
get_outlier_runs <- function(outlier_df, min_windows = 5) {
  if (nrow(outlier_df) == 0) return(NULL)

  runs <- list()
  for (chr in unique(outlier_df$chrom)) {
    sub  <- outlier_df[outlier_df$chrom == chr, ]
    sub  <- sub[order(sub$start_pos), ]
    gaps <- c(1, diff(sub$win_index))
    run_id <- cumsum(gaps != 1)

    for (r in unique(run_id)) {
      run_rows <- sub[run_id == r, ]
      if (nrow(run_rows) >= min_windows) {
        run_rows$run_id <- r          # attach run label while we have it
        runs[[length(runs) + 1]] <- run_rows
      } else {
        cat("    Skipping run of", nrow(run_rows), "window(s) on", chr,
            "(fewer than", min_windows, "consecutive windows)\n")
      }
    }
  }
  if (length(runs) == 0) return(NULL)
  do.call(rbind, runs)
}

### PCA across all SNPs in a region ###
run_region_pca <- function(run_df, coded_matrix, win_size) {
  all_snps <- unlist(lapply(run_df$win_index, function(wi) {
    snp_start <- (wi - 1) * win_size + 1
    snp_end   <- min(wi * win_size, nrow(coded_matrix))
    snp_start:snp_end
  }))
  all_snps <- unique(sort(all_snps))

  mat      <- coded_matrix[all_snps, , drop = FALSE]   # FIX: use all_snps, not snp_start:snp_end
  row_vars <- apply(mat, 1, var)
  mat      <- mat[row_vars > 1e-10, , drop = FALSE]
  if (nrow(mat) < 3) return(NULL)

  mat_scaled <- t(scale(t(mat)))
  prcomp(t(mat_scaled), center = TRUE, scale = FALSE)
}

####  Build one PCA plot (PC1v2 + PC2v3) per region ###
make_pca_plots <- function(outlier_runs, coded_matrix, win_size, mds_axis) {
  plots <- list()
  z_col <- paste0("z_", mds_axis)

  # outlier_runs already has run_id from get_outlier_runs()
  outlier_runs$region_key <- paste(outlier_runs$chrom, outlier_runs$run_id, sep = "__")

  for (rkey in unique(outlier_runs$region_key)) {
    run_df <- outlier_runs[outlier_runs$region_key == rkey, ]
    run_df <- run_df[order(run_df$start_pos), ]

    n_wins       <- nrow(run_df)
    region_chrom <- run_df$chrom[1]                    # FIX: was reg$chrom
    region_start <- run_df$start_pos[1]                # FIX: was reg$start_pos
    region_end   <- run_df$end_pos[nrow(run_df)]       # FIX: was reg$end_pos
    mean_z       <- round(mean(run_df[[z_col]]), 3)

    cat("    ", region_chrom,
        format(region_start, big.mark = ","), "-",
        format(region_end,   big.mark = ","),
        ":", n_wins, "windows\n")

    pca_out <- run_region_pca(run_df, coded_matrix, win_size)
    if (is.null(pca_out)) next

    has_pc3 <- ncol(pca_out$x) >= 3
    var_exp <- round(100 * pca_out$sdev^2 / sum(pca_out$sdev^2), 1)

    pca_df <- data.frame(
      sample = colnames(coded_matrix),
      PC1    = pca_out$x[, 1],
      PC2    = pca_out$x[, 2]
    )
    if (has_pc3) pca_df$PC3 <- pca_out$x[, 3]

    title_str    <- paste0("Region PCA (", mds_axis, "): ", region_chrom,
                           " ", format(region_start, big.mark = ","),
                           "-",  format(region_end,   big.mark = ","))
    subtitle_str <- paste0(n_wins, " windows | mean z = ", mean_z)

    # PC1 vs PC2
    p12 <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
      geom_point(size = 2.5, alpha = 0.8, color = "steelblue") +
      theme_classic() +
      labs(title    = title_str,
           subtitle = subtitle_str,
           x        = paste0("PC1 (", var_exp[1], "%)"),
           y        = paste0("PC2 (", var_exp[2], "%)"))
    plots[[length(plots) + 1]] <- p12

    # PC2 vs PC3 — only if enough PCs exist          FIX: was always attempted
    if (has_pc3) {
      p23 <- ggplot(pca_df, aes(x = PC2, y = PC3)) +
        geom_point(size = 2.5, alpha = 0.8, color = "darkorange") +
        theme_classic() +
        labs(title    = paste0(title_str, " [PC2 vs PC3]"),
             subtitle = subtitle_str,
             x        = paste0("PC2 (", var_exp[2], "%)"),
             y        = paste0("PC3 (", var_exp[3], "%)"))
      plots[[length(plots) + 1]] <- p23
    }
  }
  plots
}

# Main loop over MDS axes
axes_with_outliers <- c()
all_outlier_rows   <- list()

for (ax in 1:n_axes) {
  mds_col <- paste0("MDS", ax)
  out_col <- paste0("out_", mds_col)
  z_col   <- paste0("z_", mds_col)

  mds_df[[z_col]]   <- (mds_df[[mds_col]] - mean(mds_df[[mds_col]])) / sd(mds_df[[mds_col]])
  mds_df[[out_col]] <- flag_outliers(mds_df[[mds_col]])

  n_out <- sum(mds_df[[out_col]])
  cat("Outlier windows on", mds_col, ":", n_out, "\n")

  if (n_out == 0) next
  axes_with_outliers <- c(axes_with_outliers, ax)

  # Manhattan-style MDS plot
  p <- ggplot(mds_df, aes(x = mid_pos / 1e6, y = .data[[mds_col]],
                           color = .data[[out_col]])) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_manual(values = c("grey60", "firebrick"),
                       labels = c("normal", "outlier")) +
    geom_hline(yintercept = mean(mds_df[[mds_col]]) + 3 * sd(mds_df[[mds_col]]),
               linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_hline(yintercept = mean(mds_df[[mds_col]]) - 3 * sd(mds_df[[mds_col]]),
               linetype = "dashed", color = "red", linewidth = 0.5) +
    facet_wrap(~chrom_label, scales = "free_x") +
    theme_classic() +
    labs(x = "Position (Mb)", y = mds_col, color = "",
         title = paste0("Outlier windows — ", mds_col))

  ggsave(paste0(output_prefix, "_", win_size, "snp_", mds_col, "_plot.png"), p)

  # Collect outlier windows for summary table
  outlier_wins <- mds_df[mds_df[[out_col]] == TRUE, ]

  all_outlier_rows[[ax]] <- data.frame(
    mds_axis  = mds_col,
    chrom     = outlier_wins$chrom,
    start_pos = outlier_wins$start_pos,
    end_pos   = outlier_wins$end_pos,
    mid_pos   = outlier_wins$mid_pos,
    win_index = outlier_wins$win_index,
    z_score   = outlier_wins[[z_col]],
    cluster   = outlier_wins$cluster
  )

  # Filter to runs of >= 5 consecutive windows
  cat("  Filtering to runs of >= 5 consecutive windows...\n")
  outlier_runs <- get_outlier_runs(outlier_wins, min_windows = 5)

  if (is.null(outlier_runs)) {
    cat("  No runs of >= 5 consecutive outlier windows for", mds_col, "- skipping PCAs\n")
    next
  }

  # Report candidate regions
  outlier_runs$region_key <- paste(outlier_runs$chrom, outlier_runs$run_id, sep = "__")
  cat("  Candidate regions for", mds_col, ":\n")
  for (rkey in unique(outlier_runs$region_key)) {
    rd <- outlier_runs[outlier_runs$region_key == rkey, ]
    cat("    ", rd$chrom[1],
        format(min(rd$start_pos), big.mark = ","), "-",
        format(max(rd$end_pos),   big.mark = ","),
        ":", nrow(rd), "windows\n")
  }

  # Make and save PCA plots (one per region)
  pca_plots <- make_pca_plots(outlier_runs, coded, win_size, mds_col)

  if (length(pca_plots) > 0) {
    pdf(paste0(output_prefix, "_", win_size, "snp_outlier_pca_", mds_col, ".pdf"),
        width = 7, height = 6)
    for (pp in pca_plots) print(pp)
    dev.off()
    n_regions <- length(unique(outlier_runs$region_key))
    cat("  Saved", length(pca_plots), "PCA plots for", n_regions, "regions on", mds_col, "\n")
  }
}

# Write summary table
if (length(all_outlier_rows) > 0) {
  outlier_table <- do.call(rbind, all_outlier_rows)
  outlier_table <- outlier_table[order(outlier_table$mds_axis,
                                       outlier_table$chrom,
                                       outlier_table$start_pos), ]

  outlier_table$run_id             <- NA
  outlier_table$n_consecutive_windows <- NA

  for (ax in unique(outlier_table$mds_axis)) {
    for (chr in unique(outlier_table$chrom[outlier_table$mds_axis == ax])) {
      idx  <- which(outlier_table$mds_axis == ax & outlier_table$chrom == chr)
      sub  <- outlier_table[idx, ]
      gaps <- c(1, diff(sub$win_index))
      runs <- cumsum(gaps != 1)
      outlier_table$run_id[idx] <- runs
      for (r in unique(runs)) {
        ridx <- idx[runs == r]
        outlier_table$n_consecutive_windows[ridx] <- length(ridx)
      }
    }
  }

  write.table(outlier_table,
              file      = paste0(output_prefix, "_", win_size, "snp_outlier_windows.txt"),
              sep       = "\t",
              row.names = FALSE,
              quote     = FALSE)
  cat("\nOutlier window summary table written with", nrow(outlier_table), "rows\n")
}

cat("\nAxes with outliers:   ", paste(paste0("MDS", axes_with_outliers), collapse = ", "), "\n")
cat("Axes without outliers:", paste(paste0("MDS", setdiff(1:n_axes, axes_with_outliers)), collapse = ", "), "\n")
