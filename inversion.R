
# load libraries ----------------------------------------------------------


library(tidyverse)
library(dplyr)
library(fuzzyjoin)
library(RColorBrewer)
library(ggvenn)
library(ggplot2)
library(cowplot)
library(scales)
library(GenomicRanges)
library(ggpubr)
library(purrr)
library(dunn.test)
setwd("/Users/a02499139//Desktop/Gompert_Lab_Research/TimemaSVmethods/analyses")


### PANGENOME ###-----------------------------------------------------------------
  # load pantree inversion coordinates --------------------------------------
df <- read_tsv("/Users/a02499139//Desktop/Gompert_Lab_Research/TimemaSVmethods/PantreeSummaries/inversion_coordinates.tsv")  

scaffold_info <- df %>%
  distinct(scaffold) %>%
  mutate(
    scaffold_num = as.integer(str_extract(scaffold, "(?<=Scaffold_)\\d+")),
    scaf_length  = as.numeric(str_extract(scaffold, "(?<=length_)\\d+"))
  ) %>%
  arrange(scaffold_num) %>%
  mutate(scaffold = str_extract(scaffold, "Scaffold_\\d+"))

#add chromosome name from science paper
chromosome_map <- c(
  "1"  = 4,
  "2"  = 3,
  "3"  = 13,
  "4"  = 8,
  "5" = 6,
  "6" = 2,
  "7" = 10,
  "8" = 7,
  "9" = 9,
  "10" = 11,
  "11" = 12,
  "13" = 1)
scaffold_info$scaffold_num <- trimws(scaffold_info$scaffold_num)
scaffold_info$Chromosome <- chromosome_map[scaffold_info$scaffold_num]

#add chromosome 5
chr5<- data.frame(
  scaffold    = "Scaffold_12",
  scaffold_num = 12,
  scaf_length  = 47609450,
  Chromosome   = 5
)

scaffold_info <- rbind(scaffold_info, chr5)

scaffold_info <-scaffold_info %>%
  arrange(Chromosome) %>%
  dplyr::mutate(y = row_number())

df <- df %>%
  mutate(scaffold = str_extract(scaffold, "Scaffold_\\d+")) %>%
  mutate(
    start_pos = as.numeric(start_pos),
    end_pos   = as.numeric(end_pos)
  )


df <- df  %>% 
  separate(genome,
           into = c("genome_id", "haplotype", "scaf_from_genome", "extra"),
           sep  = "#",
           extra = "merge",
           fill  = "right") %>%
  mutate(
    haplotype = as.integer(haplotype),              # 1 or 2
    # scaffold to plot on (use the 'scaffold' column, not scaf_from_genome)
    scaffold  = factor(scaffold,
                       levels = paste0("Scaffold_", 1:13))  # fix order
  ) 

genomes      <- unique(df$genome_id)
n_genomes    <- length(genomes)

df <- df %>%
  mutate(
    genome_idx = as.integer(factor(genome_id, levels = genomes)),
    hap_offset = (genome_idx - 1) * 0.4 + (haplotype - 1) * 0.18,
    y_pos      = as.integer(scaffold) + hap_offset - 0.3   # centred around scaffold integer
  )

  # Load pantree inversion sizes ----------------------------------------------------

sizes <- read.csv("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/PantreeSummaries/all_scaffolds_inversions_tcrist_genotypes_NRbp.csv")  

colnames(sizes)<-c("scaffold", "variant_id", "RC", "AC", "TP", "NIA", "AN", "NR_bp", 
                   "Hap2_hwy154_cen4119", "hwy154_cen4119.1", "hwy154_cen4280.1", "hwy154_cen4280.2",
                   "refug_cen4120.1", "refug_cen4120.2", "refug_cen4122.1", "refug_cen4122.2")
sizes <- sizes %>%
  mutate(scaffold = str_extract(scaffold, "Scaffold_\\d+"))

df_sizes<-left_join(df,sizes, by=c("scaffold","variant_id"))

df_sizes  <- df_sizes  %>%
  mutate(
    start_pos = as.numeric(start_pos),
    end_pos   = as.numeric(end_pos)
  )

df_sizes$pos_diff<- abs(df_sizes$end_pos - df_sizes$start_pos)
df_sizes<-left_join(df_sizes, scaffold_info, by="scaffold")
df_sizes <- df_sizes %>%
  mutate(size = pmax(pos_diff, NR_bp, na.rm = TRUE))

df_sizes_small<-df_sizes %>%
  dplyr::select(c("genome_id","haplotype","variant_id", "Chromosome", "start_pos","end_pos", "NR_bp", "pos_diff", "size"))

# Parse genome column to remove rows with missing start_pos or end_pos
df_parsed_sizes <- df_sizes %>%
  filter(start_pos != ".") 

#filter it down to just one unique Chromosome & variant_id combo, keep the row with the lowest available genome_idx &
# create a new column called position, where if the genome_idx=1 it says "exact" and if it is another number it says "approx"

df_parsed_sizes_oneofeach <- df_parsed_sizes %>%
  arrange(Chromosome, variant_id, genome_idx) %>%
  group_by(Chromosome, variant_id) %>%
  slice_min(genome_idx, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(position = if_else(genome_idx == 1, "exact", "approx")) 

#this makes start and stop end/start_pos first and NR_BP second
df_parsed_sizes_oneofeach_startend <- df_parsed_sizes_oneofeach %>%
  mutate(
    start = case_when(
      !is.na(end_pos) ~ pmin(start_pos, end_pos),
      is.na(end_pos) ~ start_pos),
    stop = case_when(
      !is.na(end_pos) ~ pmax(start_pos, end_pos),
      is.na(end_pos) & NR_bp != 0 ~ start_pos + NR_bp
    )
  )

df_parsed_sizes_oneofeach_small_startend<-df_parsed_sizes_oneofeach_startend%>%
  dplyr::select(c("genome_id","haplotype", "Chromosome", "start_pos","end_pos", "NR_bp", "pos_diff", "start", "stop"))

  # Parse ref/alt ------------------------

#  parse a GT:CR:CA string 
parse_genotype <- function(x) {
  # Check before as.character() so is.na() works correctly
  if (is.na(x) || x == "") return(list(GT = NA, CR = NA, CA = NA, call = "missing"))
  
  x <- as.character(x)
  
  if (x == ".") return(list(GT = NA, CR = NA, CA = NA, call = "missing"))
  
  parts <- strsplit(x, ":")[[1]]
  
  if (length(parts) < 3) return(list(GT = NA, CR = NA, CA = NA, call = "missing"))
  
  GT <- suppressWarnings(as.integer(parts[1])) #genotype
  CR <- suppressWarnings(as.integer(parts[2])) #times visiting REF
  CA <- suppressWarnings(as.integer(parts[3]))# times visiting ALT
  
  call <- dplyr::case_when(
    is.na(GT)                                    ~ "missing",
    GT == 1                                      ~ "ALT",
    GT == 0                                      ~ "REF",
    !is.na(CA) & CA > 0 & (is.na(CR) | CR == 0) ~ "ALT",
    !is.na(CR) & CR > 0 & (is.na(CA) | CA == 0) ~ "REF",
    !is.na(CR) & !is.na(CA) & CR > 0 & CA > 0   ~ "ambiguous",
    TRUE                                         ~ "missing"
  )
  
  list(GT = GT, CR = CR, CA = CA, call = call)
}

# Identify haplotype columns 
hap_cols <- names(df_sizes)[20:27]

# Find rows where start_pos or end_pos is missing
target_rows <- seq_len(nrow(df_sizes))

# Parse and classify 
results <- lapply(target_rows, function(i) {
  row_result <- lapply(hap_cols, function(col) {
    parsed <- parse_genotype(df_sizes[[col]][i])  # [[col]] not [col]
    data.frame(
      row_index = i,
      haplotype = col,
      GT        = parsed$GT,
      CR        = parsed$CR,
      CA        = parsed$CA,
      call      = parsed$call,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(row_result)
})

result_long <- dplyr::bind_rows(results)

# pivot wide: one row per inversion, one col per haplotype
result_wide <- result_long %>%
  dplyr::select(row_index, haplotype, call) %>%
  tidyr::pivot_wider(names_from = haplotype, values_from = call)

result_wide <- result_wide %>%
  dplyr::mutate(
    n_genomes = rowSums(
      dplyr::across(all_of(hap_cols), ~ .x == "ALT"),
      na.rm = TRUE
    )
  )
# Merge back to original df 
df_sizes_annotated <- df_sizes[target_rows, ] %>%
  dplyr::mutate(row_index = target_rows) %>%
  dplyr::left_join(result_wide, by = "row_index") %>%
  dplyr::select(-row_index)
df_sizes_annotated$haplotype<-as.numeric(df_sizes_annotated$haplotype)

#filter to just ones that match Chromosome/genome_id/variant_id in df_parsed_sizes_oneofeach
df_sizes_annotated_oneofeach<- df_sizes_annotated %>%
  dplyr::semi_join(
    df_parsed_sizes_oneofeach_startend,
    by = c("Chromosome", "genome_id", "variant_id", "haplotype")
  )%>%
  mutate(
    start = case_when(
      !is.na(end_pos) ~ pmin(start_pos, end_pos),
      is.na(end_pos) ~ start_pos),
    stop = case_when(
      !is.na(end_pos) ~ pmax(start_pos, end_pos),
      is.na(end_pos) & NR_bp != 0 ~ start_pos + NR_bp
    )
  )

#population level frequency
df_sizes_annotated_oneofeach<- df_sizes_annotated_oneofeach%>%
  dplyr::mutate(
    has_refugio = rowSums(across(starts_with("refug"), ~ . == "ALT")) > 0,
    has_hwy154  = rowSums(across(starts_with(c("hwy154", "Hap2")), ~ . == "ALT")) > 0,
    Population  = case_when(
      has_refugio & !has_hwy154 ~ "Refugio",
      has_hwy154  & !has_refugio ~ "Hwy154",
      has_refugio & has_hwy154  ~ "Both",
      TRUE                      ~ NA_character_
    )
  ) %>%
  dplyr::select(-has_refugio, -has_hwy154)

sum(df_sizes_annotated_oneofeach$size, na.rm = TRUE)

  # Fig3c. pantree inversions across chromosomes -----------------------------

cluster_summary <- df_sizes_annotated_oneofeach %>%
  transmute(
    ref_chr       = Chromosome,  
    cluster_start = start,         
    cluster_end   = stop,  
    n_genomes = n_genomes,
    Population = Population
  )


 chr_order <- unique(scaffold_info$Chromosome) %>%
  .[order(as.numeric(gsub("[^0-9]", "", .)),
          na.last = TRUE,
          method  = "radix")]

cluster_summary <- cluster_summary%>%
  mutate(ref_chr = factor(ref_chr, levels = chr_order))


chr_lengths <- scaffold_info%>%
  dplyr::group_by(Chromosome) %>%
  dplyr::summarise(chr_len = max(scaf_length), .groups = "drop") %>%
  dplyr::mutate(ref_chr = factor(Chromosome, levels = chr_order)) %>%
  dplyr::arrange(ref_chr)

chr_extent <- chr_lengths %>%
  transmute(ref_chr = ref_chr,
            cluster_start = 0,
            cluster_end   = chr_len)

chr_layout <- chr_lengths %>%
  arrange(as.numeric(ref_chr)) 

shared_min <- 1
shared_max <- 8

SEG_Y    <- -0.45
SEG_YEND <-  0.45

cluster_summary <- cluster_summary %>%
  dplyr::mutate(ref_chr = factor(ref_chr, levels = chr_order))%>%
  dplyr::mutate(ref_chr = as.character(gsub("^Chr", "", ref_chr))) %>%
  left_join(chr_layout %>% 
  dplyr::select(ref_chr, Chromosome), by = c("ref_chr"))  %>%
  dplyr::mutate(n_genomes = as.character(n_genomes))


genome_colors <- c(
  "1" = "#E69F00",
  "2" = "#56B4E9",
  "3" = "#009E73",
  "4" = "#CC79A7",
  "6" = "#56B4E9"
)


p4 <- ggplot() +
  # Chromosome background blocks
  geom_rect(
    data = chr_layout,
    aes(xmin = 0, xmax = chr_len,
        ymin =Chromosome + SEG_Y, ymax =Chromosome + SEG_YEND),
    fill = "grey92", colour = "grey70", linewidth = 0.3
  ) +
  geom_rect(
    data = cluster_summary,
    aes(xmin = cluster_start, xmax = cluster_end,
        ymin =  Chromosome + SEG_Y,     ymax =  Chromosome + SEG_YEND,
        fill = n_genomes),
    color = NA,
    alpha = 0.92
  ) +
  scale_fill_manual(
    name   = "# genomes\nsharing",
    values = genome_colors,
    labels = c("1", "2", "3", "4", "6")
  ) +
  scale_x_continuous(
    name   = "Reference position (Mb)",
    labels = function(x) comma(x / 1e6, accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = chr_layout$Chromosome,
    labels = chr_layout$Chromosome
  ) +
  labs(title = "Pantree Inversions across genomes", y = NULL) +
  theme_cowplot(11) +
  theme(
    strip.text         = element_text(face = "bold", size = 9),
    panel.spacing      = unit(0.4, "lines"),
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    legend.position    = "right",
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
print(p4)

ggsave("fig3c_shared_inversions_bychromosome_pantree.svg",
       p4, width = 12, height = 10)

p6a <- ggplot() +
  geom_rect(data = scaffold_info_gbs_full,
            aes(xmin = 0, xmax = scaf_length,
                ymin = y - 0.45, ymax = y + 0.45),
            fill = "grey92", colour = "grey70", linewidth = 0.3) +
  geom_rect(
    data = cluster_summary,
    aes(xmin = cluster_start, xmax = cluster_end,
        ymin = y + SEG_Y,     ymax = y + SEG_YEND,
        fill = Population),
    color = NA,
    alpha = 0.92
  ) +
  scale_fill_manual(values = c("#E69F00","#61D04F","#2297E6"), name = "Population") +
  scale_x_continuous(
    name   = "Reference position (Mb)",
    labels = function(x) comma(x / 1e6, accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = scaffold_info_gbs_full$y,
    labels = scaffold_info_gbs_full$Chromosome
  ) +
  labs(title = "Pantree Inversions across genomes", y = NULL) +
  theme_cowplot(11) +
  theme(
    strip.text         = element_text(face = "bold", size = 9),
    panel.spacing      = unit(0.4, "lines"),
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    legend.position    = "right",
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
print(p6a)

ggsave("fig_shared_inversions_bypop_pantree.svg",
       p6a, width = 12, height = 10)

# histogram of population frequency

p6b <- ggplot(df_sizes_annotated_oneofeach,
              aes(x=Population,fill = Population)) +
  geom_bar() +
  scale_y_continuous(name = "Number of inversion clusters") +
  scale_fill_manual(values = c("#E69F00","#61D04F","#2297E6"), name = "Population") +
  labs(title = "Inversion frequency by population") +
  theme_cowplot(11) +
  theme(legend.position  = "none",
        strip.background = element_rect(fill = "grey92"),
        plot.subtitle    = element_text(size = 9, color = "grey40"))

print(p6b)
ggsave("fig_pantree_freqeuncyofpop.svg",
       p6b, width = 8, height = 8)

  # Fig. Pantree Size by Frequency ---------------------------------------

df_sizes_annotated_oneofeach<-df_sizes_annotated_oneofeach%>%
  mutate(
    frequency = case_when(
      n_genomes == 1               ~ "Unique (1 genome)",
      n_genomes == 8 ~ sprintf("Core (%d genomes)", n_genomes),
      TRUE                         ~ sprintf("Partial (%d genomes)", n_genomes)
    )
  )


frequency_levels <- c(
  "Unique (1 genome)",
  paste0("Partial (", 2:7, " genomes)")
)

frequency_levels <- intersect(frequency_levels, unique(df_sizes_annotated_oneofeach$frequency))
n_levels <- length(frequency_levels)

size_palette <- setNames( c(
  "#E69F00",
  colorRampPalette(c("#56B4E9", "#009E73"))(max(1, n_levels - 2)),
  "#CC79A7"
), nm = frequency_levels )

cluster_plot <- df_sizes_annotated_oneofeach %>%
  mutate(frequency = factor(frequency, levels = frequency_levels))



p2a <- ggplot(cluster_plot,
              aes(x = abs(pos_diff) / 1000, fill = frequency)) +
  geom_histogram(bins = 60, color = "white", linewidth = 0.2, alpha = 0.9) +
  scale_x_log10(
    name   = "Inversion size — max per cluster (kb, log scale)",
    labels = comma_format(accuracy = 0.1)
  ) +
  scale_y_continuous(name = "Number of inversion clusters") +
  scale_fill_manual(values = size_palette) +
  facet_wrap(~ frequency, scales = "free_y", ncol = 1) +
  labs(title = "Inversion size distribution by frequency") +
  theme_cowplot(11) +
  theme(legend.position  = "none",
        strip.background = element_rect(fill = "grey92"),
        plot.subtitle    = element_text(size = 9, color = "grey40"))

print(p2a)
ggsave("fig2a_inversion_size_histogram_pantree.svg",
       p2a, width = 8, height = 3 * n_levels)

  # Fig. Inversion length vs. frequency ####

p2b <- ggplot(cluster_plot,
              aes(x = frequency, y = cluster_size / 1000,
                  fill = frequency, color = frequency)) +
  geom_violin(alpha = 0.35, linewidth = 0.7,
              quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.18, alpha = 0.5, size = 1.2, shape = 16) +
  scale_y_log10(
    name   = "Inversion size — max per cluster (kb, log scale)",
    labels = comma_format(accuracy = 0.1)
  ) +
  scale_x_discrete(name = "Frequency across genomes") +
  scale_fill_manual(values = size_palette) +
  scale_color_manual(values = size_palette) +
  labs(title    = "Inversion length vs. frequency") +
  theme_cowplot(12) +
  theme(legend.position = "none",
        axis.text.x     = element_text(angle = 30, hjust = 1),
        plot.subtitle   = element_text(size = 9, color = "grey40"))

ggsave(file.path(OUT_DIR, "fig_length_vs_sharedness_RO.svg"),
       p2b, width = 8, height = 6)



  # Load INVPG inversions ---------------------------------------------------
invpg <- read_tsv("/Users/a02499139//Desktop/Gompert_Lab_Research/TimemaSVmethods/INVPG_annot/invpg_HWY154_REF_4119Hap2_sumbp.tsv")  

invpg <- invpg  %>%
  mutate(
   scaffold_num = as.character(str_extract(CHROM, "(?<=Scaffold_)\\d+")),
   pos_end=as.numeric(POS+REF_bp))

invpg<-left_join(invpg,scaffold_info, by=c("scaffold_num"))

#plot
invpg_plot <- ggplot() +
  # Chromosome background blocks
  geom_rect(
    data = chr_layout,
    aes(xmin = 0, xmax = chr_len,
        ymin = y + SEG_Y, ymax = y + SEG_YEND),
    fill = "grey92", colour = "grey70", linewidth = 0.3
  ) +
  geom_rect(
    data = invpg,
    aes(xmin = POS, xmax = pos_end,
        ymin = y + SEG_Y,     ymax = y + SEG_YEND)
  ) +
  scale_x_continuous(
    name   = "Reference position (Mb)",
    labels = function(x) comma(x / 1e6, accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = chr_layout$y,
    labels = chr_layout$y
  ) +
  labs(title = "INVPG_annot Inversions across genomes", y = NULL) +
  theme_cowplot(11) +
  theme(
    strip.text         = element_text(face = "bold", size = 9),
    panel.spacing      = unit(0.4, "lines"),
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    legend.position    = "right",
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
print(invpg_plot)

ggsave("invpg_plot.svg",invpg_plot, width = 12, height = 10)

### GBS ### ###---------------------------------------------------------------------
  # load GBS local pca windows prefiltering -------------------------------------

#load GBS windows file
ref_gbs <- read_delim("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/GBS/REF/REF_all_100snp_outlier_windows_40mds.txt")
hwy_gbs <- read_delim("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/GBS/HWY154(FHA)/FHA_all_100snp_outlier_windows_40mds.txt")

ref_gbs$pop<-"REF"
hwy_gbs$pop<-"HWY"



  # load and view post-filtering GBS windows -----------------------------------

#ref
ref_merged_100snp <- read_delim("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/GBS/REF/REF_all_100snp_putative_inversions.txt")
ref_merged_50snp <- read_delim("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/GBS/REF/REF_all_50snp_putative_inversions.txt")

#merge positive/negative direction
ref_directionmerged_100snp<-ref_merged_100snp%>%
  group_by(chrom, start_pos, end_pos) %>%
  summarise(
    start_pos = min(start_pos),
    end_pos   = max(end_pos),
    n_windows = sum(n_windows),
    perm_p    = min(perm_p),  # keep most significant p value
    direction = if (n_distinct(direction) > 1) "both" else direction[1],
    .groups   = "drop"
  )

ref_MDSmerged_100snp<-ref_merged_100snp%>%
  dplyr::arrange(chrom, start_pos) %>%
  dplyr::group_by(chrom) %>%
  dplyr::mutate(overlap_group = cumsum(cummax(lag(end_pos, default = 0)) < start_pos)) %>%
  dplyr::group_by(chrom, overlap_group) %>%
  dplyr::summarise(
    start_pos = min(start_pos),
    end_pos   = max(end_pos),
    n_windows = sum(n_windows),
    perm_p    = min(perm_p),
    direction = if (n_distinct(direction) > 1) "both" else direction[1],
    mds_axes  = paste(sort(unique(mds_axis)), collapse = ","),
    .groups   = "drop"
  ) %>%
  dplyr::select(-overlap_group)

#plot ref
ref_MDSmerged_100snp  <- ref_MDSmerged_100snp  %>%
  dplyr::mutate(Chromosome = as.numeric(str_extract(chrom, "\\d+")))%>%
  dplyr::left_join(scaffold_info %>% dplyr::select(Chromosome, y,scaf_length), by = "Chromosome")%>%
  dplyr::mutate(Chromosome = factor(Chromosome, 
                             levels = unique(Chromosome[order(as.numeric(Chromosome))])))

ref_MDSmerged_100snp_plot <- ggplot(ref_MDSmerged_100snp) +
  # Scaffold rectangles
  geom_rect(data = scaffold_info,
            aes(xmin = 0, xmax = scaf_length,
                ymin = y - 0.45, ymax = y + 0.45),
            fill = "grey92", colour = "grey70", linewidth = 0.3) +
  # One segment per window, coloured by mds_axis, offset by pop
  geom_segment(aes(x     = start_pos,
                   xend  = end_pos,
                   y     = y,
                   yend  = y),
               linewidth = 1.8,
               lineend   = "round") +
  scale_y_continuous(
    breaks = scaffold_info$y,
    labels = scaffold_info$Chromosome
  ) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
  labs(title = "GBS outlier windows across scaffolds",
       x = "Position (Mb)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    axis.text.y        = element_text(hjust = 1)
  )

print(ref_MDSmerged_100snp_plot)
ggsave("ref_MDSmerged_100snp_plot.svg", ref_MDSmerged_100snp_plot, width = 14, height = 8)

# Colour palette: one hue per mds_axis
mds_axes   <- unique(ref_merged_100snp$mds_axis)
n_axes     <- length(mds_axes)
axis_colors <- setNames(
  hcl.colors(n_axes, palette = "Dark 2"),
  mds_axes
)


ref_merged_100snp_plot <- ggplot(ref_merged_100snp) +
  # Scaffold rectangles
  geom_rect(data = scaffold_info,
            aes(xmin = 0, xmax = scaf_length,
                ymin = y - 0.45, ymax = y + 0.45),
            fill = "grey92", colour = "grey70", linewidth = 0.3) +
  # One segment per window, coloured by mds_axis, offset by pop
  geom_segment(aes(x     = start_pos,
                   xend  = end_pos,
                   y     = y,
                   yend  = y,
                   colour = mds_axis),
               linewidth = 1.8,
               lineend   = "round") +
  scale_colour_manual(values = axis_colors) +
  scale_y_continuous(
    breaks = scaffold_info$y,
    labels = scaffold_info$Chromosome
  ) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
  labs(title = "GBS outlier windows across scaffolds",
       x = "Position (Mb)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    axis.text.y        = element_text(hjust = 1)
  )

print(ref_merged_100snp_plot)
ggsave("ref_merged_100snp_plot.svg", ref_merged_100snp_plot, width = 14, height = 8)

# FHA
fha_merged_100snp <- read_delim("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/GBS/HWY154(FHA)/FHA_all_100snp_putative_inversions.txt")

#merge positive/negative direction
fha_directionmerged_100snp<-fha_merged_100snp%>%
  group_by(chrom, start_pos, end_pos) %>%
  summarise(
    start_pos = min(start_pos),
    end_pos   = max(end_pos),
    n_windows = sum(n_windows),
    perm_p    = min(perm_p),  # keep most significant p value
    direction = if (n_distinct(direction) > 1) "both" else direction[1],
    .groups   = "drop"
  )

fha_MDSmerged_100snp<-fha_merged_100snp%>%
  dplyr::arrange(chrom, start_pos) %>%
  dplyr::group_by(chrom) %>%
  dplyr::mutate(overlap_group = cumsum(cummax(lag(end_pos, default = 0)) < start_pos)) %>%
  dplyr::group_by(chrom, overlap_group) %>%
  dplyr::summarise(
    start_pos = min(start_pos),
    end_pos   = max(end_pos),
    n_windows = sum(n_windows),
    perm_p    = min(perm_p),
    direction = if (n_distinct(direction) > 1) "both" else direction[1],
    mds_axes  = paste(sort(unique(mds_axis)), collapse = ","),
    .groups   = "drop"
  ) %>%
  dplyr::select(-overlap_group)


fha_MDSmerged_100snp  <- fha_MDSmerged_100snp  %>%
  dplyr::mutate(Chromosome = as.numeric(str_extract(chrom, "\\d+")))%>%
  dplyr::left_join(scaffold_info %>% 
  dplyr::select(Chromosome, y,scaf_length), by = "Chromosome")%>%
  dplyr::mutate(Chromosome = factor(Chromosome, 
                             levels = unique(Chromosome[order(as.numeric(Chromosome))])))

#plot fha

fha_MDSmerged_100snp_plot <- ggplot(fha_MDSmerged_100snp) +
  # Scaffold rectangles
  geom_rect(data = scaffold_info,
            aes(xmin = 0, xmax = scaf_length,
                ymin = y - 0.45, ymax = y + 0.45),
            fill = "grey92", colour = "grey70", linewidth = 0.3) +
  # One segment per window, coloured by mds_axis, offset by pop
  geom_segment(aes(x     = start_pos,
                   xend  = end_pos,
                   y     = y,
                   yend  = y),
               linewidth = 1.8,
               lineend   = "round") +
  scale_y_continuous(
    breaks = scaffold_info$y,
    labels = scaffold_info$Chromosome
  ) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
  labs(title = "GBS outlier windows across scaffolds",
       x = "Position (Mb)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    axis.text.y        = element_text(hjust = 1)
  )

print(fha_MDSmerged_100snp_plot)
ggsave("fha_MDSmerged_100snp_plot.svg", fha_MDSmerged_100snp_plot, width = 14, height = 8)

# Colour palette: one hue per mds_axis
mds_axes   <- unique(fha_merged_100snp$mds_axis)
n_axes     <- length(mds_axes)
axis_colors <- setNames(
  hcl.colors(n_axes, palette = "Dark 2"),
  mds_axes
)


fha_merged_100snp_plot <- ggplot(fha_merged_100snp) +
  # Scaffold rectangles
  geom_rect(data = scaffold_info,
            aes(xmin = 0, xmax = scaf_length,
                ymin = y - 0.45, ymax = y + 0.45),
            fill = "grey92", colour = "grey70", linewidth = 0.3) +
  # One segment per window, coloured by mds_axis, offset by pop
  geom_segment(aes(x     = start_pos,
                   xend  = end_pos,
                   y     = y,
                   yend  = y,
                   colour = mds_axis),
               linewidth = 1.8,
               lineend   = "round") +
  scale_colour_manual(values = axis_colors) +
  scale_y_continuous(
    breaks = scaffold_info$y,
    labels = scaffold_info$Chromosome
  ) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
  labs(title = "GBS outlier windows across scaffolds",
       x = "Position (Mb)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    axis.text.y        = element_text(hjust = 1)
  )

print(fha_merged_100snp_plot)
ggsave("fha_merged_100snp_plot.svg", fha_merged_100snp_plot, width = 14, height = 8)


  # merge and sort GBS data -------------------------------------------------
all_gbs<-rbind(ref_gbs, hwy_gbs)
scaffold_info$y <- rank(-scaffold_info$Chromosome)

all_gbs <- all_gbs %>%
  dplyr::mutate(Chromosome = as.numeric(str_extract(chrom, "\\d+")))%>%
  dplyr::left_join(scaffold_info %>% 
                     dplyr::select(Chromosome, y,scaf_length), by = "Chromosome")%>%
  dplyr::select(-c("mid_pos", "cluster","run_id")) %>%
  dplyr::mutate(Chromosome = factor(Chromosome, 
                             levels = unique(Chromosome[order(as.numeric(Chromosome))])))


# filter to just runs of more than 5 outlier windows
valid_runs <- all_gbs %>%
  group_by(Chromosome, n_consecutive_windows,pop,y, mds_axis) %>%
  filter(n_consecutive_windows >= 5) %>%
  ungroup()

#collapse consecutive windows
merged_runs <- valid_runs %>%
  dplyr::group_by(Chromosome, n_consecutive_windows,pop,y, mds_axis) %>%
  dplyr::summarise(
    start_pos = min(start_pos),
    end_pos   = max(end_pos),
    n_windows = n(),
    .groups = "drop"
  )

#merge filtered datasets
ref_MDSmerged_100snp$pop<-"REF"
fha_MDSmerged_100snp$pop<-"HWY"

all_gbs_filtered<-rbind(ref_MDSmerged_100snp, fha_MDSmerged_100snp)

#merge if any overlap, results in 26 inversions
all_gbs_filtered_mergedany <- all_gbs_filtered %>%
  dplyr::arrange(chrom, start_pos) %>%
  dplyr::group_by(chrom) %>%
  dplyr::mutate(overlap_group = cumsum(cummax(lag(end_pos, default = 0)) < start_pos)) %>%
  dplyr::group_by(chrom, overlap_group) %>%
  dplyr::summarise(
    start_pos = min(start_pos),
    end_pos   = max(end_pos),
    n_windows = sum(n_windows),
    perm_p    = min(perm_p),
    direction = if (n_distinct(direction) > 1) "both" else direction[1],
    pops  = paste(sort(unique(pop)), collapse = ","),
    .groups   = "drop"
  ) %>%
  dplyr::select(-overlap_group)

#merge if 2/3 80% rule is followed
reciprocal_overlap <- function(s1, e1, s2, e2, size_thresh = 2/3, overlap_thresh = 0.8) {
  overlap <- max(0, min(e1, e2) - max(s1, s2))
  len1    <- e1 - s1
  len2    <- e2 - s2
  
  smaller <- min(len1, len2)
  larger  <- max(len1, len2)
  
  # smaller must be at least 2/3 the size of larger
  if (smaller / larger < size_thresh) return(0)
  
  # smaller region must be 80% overlapping with itself
  overlap / smaller >= overlap_thresh
}
merge_across_pops <- function(inv_df, overlap_thresh = 0.8) {
  inv_df  <- inv_df %>% arrange(chrom, start_pos)
  n       <- nrow(inv_df)
  group   <- seq_len(n)  # each row starts in its own group
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (inv_df$chrom[i] != inv_df$chrom[j]) next
      
      ro <- reciprocal_overlap(inv_df$start_pos[i], inv_df$end_pos[i],
                               inv_df$start_pos[j], inv_df$end_pos[j])
      if (ro >= overlap_thresh) {
        # merge j into i's group
        group[group == group[j]] <- group[i]
      }
    }
  }
  
  inv_df$group <- group
  
  inv_df %>%
    dplyr::group_by(chrom, group) %>%
    dplyr::summarise(
      start_pos = min(start_pos),
      end_pos   = max(end_pos),
      n_windows = sum(n_windows),
      perm_p    = min(perm_p),
      direction = if (n_distinct(direction) > 1) "both" else direction[1],
      mds_axes  = paste(sort(unique(mds_axes)), collapse = ","),
      pops      = paste(sort(unique(pop)), collapse = ","),
      .groups   = "drop"
    ) %>%
    dplyr::select(-group)
}

all_gbs_filtered_merged80 <- merge_across_pops( all_gbs_filtered, overlap_thresh = 0.8)
  # See how many new inversions MDS axes contribute -------------------------

axes <- paste0("MDS", 1:40)

# For each axis in order, find windows not overlapping any previous axis
seen <- data.frame(Chromosome = character(), start_pos = numeric(), end_pos = numeric())
new_counts <- numeric(length(axes))

for (i in seq_along(axes)) {
  current <- merged_runs %>% filter(mds_axis == axes[i])
  
  if (nrow(seen) == 0) {
    new_wins <- current
  } else {
    # For each window, check if it overlaps anything in seen
    new_wins <- current %>%
      rowwise() %>%
      filter(!any(seen$Chromosome == Chromosome &
                    seen$start_pos <= end_pos &
                    seen$end_pos   >= start_pos)) %>%
      ungroup()
  }
  
  new_counts[i] <- nrow(new_wins)
  seen <- bind_rows(seen, new_wins %>% select(Chromosome, start_pos, end_pos))
}

# Build plotting dataframe
cumulative_df <- data.frame(
  mds_axis  = axes,
  new_wins  = new_counts,
  cumulative = cumsum(new_counts)
)


newinversions_mds<-ggplot(cumulative_df, aes(x = 1:40, y = new_wins)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, colour = "steelblue") +
  scale_x_continuous(breaks = 1:40, labels = axes) +
  labs(x = "MDS Axis", y = "New outlier windows",
       title = "New outlier windows contributed by each MDS axis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("gbs_newinversions_by_mdsaxes.svg", newinversions_mds, width = 10, height = 8)

  # plot all GBS inversions -------------------------------------------------

# Assign numeric y offset per scaffold × pop combination
pop_levels <- unique(all_gbs_filtered$pop)
n_pops     <- length(pop_levels)

pop_offsets <- tibble(
  pop      = pop_levels,
  pop_idx  = seq_along(pop_levels),
  y_offset = scales::rescale(seq_along(pop_levels), to = c(-0.35, 0.35))
)

# Join offsets onto all_gbs
all_gbs_filtered_plot <- all_gbs_filtered %>%
  left_join(pop_offsets, by = "pop") %>%
  mutate(y_pos = y + y_offset)

# Join offsets onto merged_runs
merged_runs_plot <- merged_runs %>%
  left_join(pop_offsets, by = "pop") %>%
  mutate(y_pos = y + y_offset)





p_gbs2 <- ggplot(merged_runs_plot) +
  geom_rect(data = scaffold_info,
            aes(xmin = 0, xmax = scaf_length,
                ymin = y - 0.45, ymax = y + 0.45),
            fill = "grey92", colour = "grey70", linewidth = 0.3) +
  geom_segment(aes(x      = start_pos,
                   xend   = end_pos,
                   y      = y_pos,
                   yend   = y_pos,
                   colour = pop),
               linewidth = 1.8,
               lineend   = "round") +
  scale_colour_manual(values = c("#61D04F","#2297E6"), name = "Population") +
  scale_y_continuous(
    breaks = scaffold_info$y,
    labels = scaffold_info$Chromosome
  ) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
  labs(title = "GBS outlier windows across scaffolds",
       x = "Position (Mb)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    axis.text.y        = element_text(hjust = 1)
  )
print(p_gbs2)
ggsave("gbs_min5outlierwindows_by_pop_mds_sciencepapernames.pdf", p_gbs2, width = 14, height = 8)
ggsave("gbs_min5outlierwindows_by_pop_mds_sciencepapernames.svg", p_gbs2, width = 14, height = 8)




### COMPARATIVE ALIGNMENT ### ###---------------------------------------------------
  # load syri inversions ----------------------------------------------------
all_inv_raw <- read_tsv("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/syri/all_inversions_raw.tsv")

#split up 3&4 before clustering based off of orientation
chr3_length <- scaffold_info$scaf_length[scaffold_info$Chromosome == "3"]
chr4_length <- scaffold_info$scaf_length[scaffold_info$Chromosome == "4"]

remap_3_4 <- function(df, version) {
  df %>%
    dplyr::filter(ref_chr == "Chr3_4") %>%
    dplyr::mutate(
      ref_chr = case_when(
        genome == "TcrRGUS1" & ref_start <  chr3_length ~ "Chr3",
        genome == "TcrRGUS1" & ref_start >= chr3_length ~ "Chr4",
        genome %in% c("TcrRGUS2", "TcrRGS2") & ref_start <  chr4_length ~ "Chr4",
        genome %in% c("TcrRGUS2", "TcrRGS2") & ref_start >= chr4_length ~ "Chr3"
      ),
      new_start = case_when(
        # RGUS1: chr3 side — unchanged
        genome == "TcrRGUS1" & ref_start <  chr3_length ~ ref_start,
        # RGUS1: chr4 side — revcomp flip
        genome == "TcrRGUS1" & ref_start >= chr3_length ~ chr4_length - (ref_end   - chr3_length),
        # RGUS2/RGS2: chr4 side — unchanged
        genome %in% c("TcrRGUS2", "TcrRGS2") & ref_start <  chr4_length ~ ref_start,
        # RGUS2/RGS2: chr3 side — revcomp flip
        genome %in% c("TcrRGUS2", "TcrRGS2") & ref_start >= chr4_length ~ chr3_length - (ref_end   - chr4_length)
      ),
      new_end = case_when(
        genome == "TcrRGUS1" & ref_start <  chr3_length ~ ref_end,
        genome == "TcrRGUS1" & ref_start >= chr3_length ~ chr4_length - (ref_start - chr3_length),
        genome %in% c("TcrRGUS2", "TcrRGS2") & ref_start <  chr4_length ~ ref_end,
        genome %in% c("TcrRGUS2", "TcrRGS2") & ref_start >= chr4_length ~ chr3_length - (ref_start - chr4_length)
      ),
      ref_start = new_start,
      ref_end   = new_end
    ) %>%
    dplyr::select(-new_start, -new_end)
}

all_inv_raw_remap <- all_inv_raw  %>%
  filter(ref_chr != "Chr3_4") %>%
  bind_rows(remap_3_4(all_inv_raw)) 


  # RECIPROCAL OVERLAP ANALYSIS ####
# For each genome G and RO threshold t, count how many of G's
# inversions have RO >= t with at least one inversion in any
# other genome. RO(A, B) = overlap_length / min(len_A, len_B)
# This is the workflow:
#   1. Filter by candidate pairs in which the smallest one is at least two thirds the size of the biggest
#   2. Compute true RO on coordinates
#   3. Count focal inversions whose best RO >= threshold t

compute_best_ro <- function(focal_df, other_df) {
  best_ro <- rep(0.0, nrow(focal_df))
  
  focal_gr <- GRanges(
    seqnames = focal_df$ref_chr,
    ranges   = IRanges(
      start = pmax(1L, focal_df$ref_start),
      end   = focal_df$ref_end
    )
  )
  other_gr <- GRanges(
    seqnames = other_df$ref_chr,
    ranges   = IRanges(start = other_df$ref_start,
                       end   = other_df$ref_end)
  )
  
  hits <- findOverlaps(focal_gr, other_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(best_ro)
  
  fi <- queryHits(hits)
  oi <- subjectHits(hits)
  
  # True overlap on original coordinates
  ov_len   <- pmax(0L,
                   pmin(focal_df$ref_end[fi], other_df$ref_end[oi]) -
                     pmax(focal_df$ref_start[fi], other_df$ref_start[oi]) + 1L)
  min_size <- pmin(focal_df$inv_size[fi], other_df$inv_size[oi])
  max_size <-  pmax(focal_df$inv_size[fi], other_df$inv_size[oi])
  SIZE_RATIO_MIN <- 2/3
  size_ratio <- min_size / max_size
  ro         <- ifelse(size_ratio >= SIZE_RATIO_MIN, ov_len / min_size, 0)
  
  # Keep best RO per focal inversion
  for (i in seq_along(fi)) {
    if (ro[i] > best_ro[fi[i]]) best_ro[fi[i]] <- ro[i]
  }
  best_ro
}

ro_thresholds <- seq(RO_MIN, RO_MAX, by = RO_STEP)

# compute best RO for each genome (against all others combined)
best_ro_list <- map(names(inv_list), function(gname) {
  focal  <- inv_list[[gname]]
  others <- bind_rows(inv_list[setdiff(names(inv_list), gname)])
  compute_best_ro(focal, others)
})
names(best_ro_list) <- names(inv_list)

sweep_results <- map_dfr(names(inv_list), function(gname) {
  bro <- best_ro_list[[gname]]
  map_dfr(ro_thresholds, function(t) {
    tibble(genome = gname, ro_thresh = t, n_shared = sum(bro >= t))
  })
})

write.table(sweep_results,
            file.path(OUT_DIR, "ro_sweep.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

  # SELECT IDEAL RO THRESHOLD ####

#normalized
sweep_normalized <- sweep_results %>%
  group_by(genome) %>%
  mutate(prop_shared = n_shared / max(n_shared)) %>%
  ungroup()

p0<-ggplot(sweep_normalized, aes(x = ro_thresh, y = prop_shared,
                                 color = genome, group = genome)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     name   = "Minimum reciprocal overlap threshold") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     name   = "Proportion of inversions with a partner") +
  scale_color_viridis_d(option = "turbo") +
  labs(title    = "Proportion of inversions shared vs. RO threshold",
       subtitle = "Curves normalized to RO=0") +
  theme_cowplot(12)

ggsave(file.path(OUT_DIR, "fig1_ro_normalized.pdf"),  p0, width = 10, height = 5.5)

  # Cluster Analysis ---------------------------------------------------------
#start with the largest inversions and work in

cluster_inversions <- function(all_inv_raw_remap, ideal_ro = 0.80, size_ratio_min = 2/3) {
  
  # make a copy, adding cluster assignment column
  df <- all_inv_raw_remap %>% mutate(cluster_id = NA_integer_)
  
  # Build GRanges
  gr <- GRanges(
    seqnames = df$ref_chr,
    ranges   = IRanges(start = df$ref_start, end = df$ref_end)
  )
  
  next_cluster <- 1L
  unassigned   <- which(is.na(df$cluster_id))
  
  while (length(unassigned) > 0) {
    
    # Take the largest unassigned inversion as the anchor
    anchor_idx <- unassigned[which.max(df$inv_size[unassigned])]
    anchor     <- df[anchor_idx, ]
    
    # Find all unassigned inversions that positionally overlap the anchor
    anchor_gr <- GRanges(
      seqnames = anchor$ref_chr,
      ranges   = IRanges(start = anchor$ref_start, end = anchor$ref_end)
    )
    candidates_gr <- gr[unassigned]
    hits <- findOverlaps(anchor_gr, candidates_gr, ignore.strand = TRUE)
    
    if (length(hits) == 0) {
      # No overlaps — anchor is its own cluster
      df$cluster_id[anchor_idx] <- next_cluster
      next_cluster <- next_cluster + 1L
      unassigned   <- which(is.na(df$cluster_id))
      next
    }
    
    candidate_local_idx <- subjectHits(hits)
    candidate_global_idx <- unassigned[candidate_local_idx]
    
    # For each candidate, compute RO with the anchor
    cands <- df[candidate_global_idx, ]
    
    ov_len <- pmax(0L,
                   pmin(anchor$ref_end,   cands$ref_end) -
                     pmax(anchor$ref_start, cands$ref_start) + 1L)
    
    min_size   <- pmin(anchor$inv_size, cands$inv_size)
    max_size   <- pmax(anchor$inv_size, cands$inv_size)
    size_ratio <- min_size / max_size
    ro         <- ov_len / min_size  
    passes <- size_ratio >= size_ratio_min & ro >= ideal_ro
    members <- c(anchor_idx, candidate_global_idx[passes])
    
    # All-pairs RO check among members
    # Keep iterating until no members are dropped
    repeat {
      member_df   <- df[members, ]
      n           <- length(members)
      if (n <= 1) break
      
      # Compute RO for every pair
      all_pass <- rep(TRUE, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          if (i == j) next
          ov <- max(0L,
                    min(member_df$ref_end[i],   member_df$ref_end[j]) -
                      max(member_df$ref_start[i], member_df$ref_start[j]) + 1L)
          min_s <- min(member_df$inv_size[i], member_df$inv_size[j])
          max_s <- max(member_df$inv_size[i], member_df$inv_size[j])
          ro_ij <- ifelse(min_s / max_s >= size_ratio_min, ov / min_s, 0)
          if (ro_ij < ideal_ro) {
            all_pass[i] <- FALSE
            break
          }
        }
      }
      
      if (all(all_pass)) break  # all members pass, done
      
      # Drop failing members and re-check
      members <- members[all_pass]
    }
    # Check that the smallest member is >= 2/3 the size of the largest
    member_sizes <- df$inv_size[members]
    cluster_size_ratio <- min(member_sizes) / max(member_sizes)
    
    if (cluster_size_ratio < size_ratio_min) {
      # Drop the smallest members until constraint is satisfied
      size_order <- order(member_sizes, decreasing = TRUE)
      members_sorted <- members[size_order]
      sizes_sorted   <- member_sizes[size_order]
      keep <- sizes_sorted >= sizes_sorted[1] * size_ratio_min
      members <- members_sorted[keep]
    }
    
    df$cluster_id[members] <- next_cluster
    next_cluster <- next_cluster + 1L
    unassigned   <- which(is.na(df$cluster_id))
  }
  
  df
}

all_inv_clustered <- cluster_inversions(all_inv_raw_remap, ideal_ro = 0.80, size_ratio_min = 2/3)

### Summaries ###
n_genomes_total <- 7

cluster_summary <- all_inv_clustered %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(
    ref_chr       = ref_chr[1],
    cluster_start = min(ref_start),
    cluster_end   = max(ref_end),
    cluster_size  = cluster_end - cluster_start + 1L,
    n_genomes     = n_distinct(genome),
    genomes_list  = paste(sort(unique(genome)), collapse = ","),
    ann_types     = paste(sort(unique(ann_type)), collapse = ","),
    mean_inv_size = mean(inv_size, na.rm = TRUE),
    max_inv_size  = max(inv_size, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  dplyr::mutate(
    frequency = case_when(
      n_genomes == 1 | n_genomes == n_genomes_total              ~ sprintf("Unique (1 genome)"),
      TRUE                         ~ sprintf("Shared (%d genomes)", n_genomes)
    )
  )

per_genome_summary <- all_inv_clustered %>%
  dplyr::left_join(cluster_summary %>% dplyr::select(cluster_id, n_genomes), by = "cluster_id") %>%
  dplyr::group_by(genome, ref_chr) %>%
  dplyr::summarise(
    n_total  = n_distinct(cluster_id),
    n_unique = n_distinct(cluster_id[n_genomes == 1 | n_genomes == 7]),
    n_shared = n_distinct(cluster_id[n_genomes  > 1]),
    .groups  = "drop"
  )

per_genome_all <- all_inv_clustered %>%
  dplyr::left_join(cluster_summary %>% dplyr::select(cluster_id, n_genomes), by = "cluster_id") %>%
  dplyr::group_by(genome) %>%
  dplyr::summarise(
    ref_chr  = "ALL",
    n_total  = n_distinct(cluster_id),
    n_unique = n_distinct(cluster_id[n_genomes == 1 | n_genomes == 7]),
    n_shared = n_distinct(cluster_id[n_genomes  > 1]),
    .groups  = "drop"
  )

per_genome_summary <- all_inv_clustered %>%
  left_join(cluster_summary %>% dplyr::select(cluster_id, n_genomes), by = "cluster_id") %>%
  dplyr::group_by(genome, ref_chr) %>%
  dplyr::summarise(
    n_total  = n_distinct(cluster_id),
    n_unique = n_distinct(cluster_id[n_genomes == 1 | n_genomes == 7]),
    n_shared = n_distinct(cluster_id[n_genomes  > 1]),
    .groups  = "drop"
  )

genome_summary <- bind_rows(per_genome_summary, per_genome_all)

#population level frequency
cluster_summary<- cluster_summary%>%
  dplyr::mutate(
    has_refugio = grepl("TcrR", genomes_list),
    has_hwy154  = grepl("TcrH", genomes_list),
    Population  = case_when(
      has_refugio & !has_hwy154 ~ "Refugio",
      has_hwy154 & !has_refugio ~ "Hwy154",
      has_refugio & has_hwy154  ~ "Both",
      TRUE                      ~ NA_character_
    )
  ) %>%
  dplyr::select(-has_refugio, -has_hwy154)


cluster_summary<-cluster_summary%>%
  filter(cluster_size>50)
table(cluster_summary$ann_types)

  # load and cluster translocations and duplications ------------------------------------
all_trans_raw <- read_tsv("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/syri/all_translocations_raw.tsv")
all_dup_raw <- read_tsv("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/syri/all_duplications_raw.tsv")

all_trans_raw_remap <- all_trans_raw  %>%
  filter(ref_chr != "Chr3_4") %>%
  bind_rows(remap_3_4(all_trans_raw)) 

all_dup_raw_remap <- all_dup_raw  %>%
  filter(ref_chr != "Chr3_4") %>%
  bind_rows(remap_3_4(all_dup_raw)) 

all_trans_clustered <- cluster_inversions(all_trans_raw_remap, ideal_ro = 0.80, size_ratio_min = 2/3)
all_dup_clustered <- cluster_inversions(all_dup_raw_remap, ideal_ro = 0.80, size_ratio_min = 2/3)

trans_cluster_summary <- all_trans_clustered %>%
  group_by(cluster_id) %>%
  summarise(
    ref_chr       = ref_chr[1],
    cluster_start = min(ref_start),
    cluster_end   = max(ref_end),
    cluster_size  = cluster_end - cluster_start + 1L,
    n_genomes     = n_distinct(genome),
    genomes_list  = paste(sort(unique(genome)), collapse = ","),
    ann_types     = paste(sort(unique(ann_type)), collapse = ","),
    mean_inv_size = mean(inv_size, na.rm = TRUE),
    max_inv_size  = max(inv_size, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  mutate(
    frequency = case_when(
      n_genomes == 1 | n_genomes == n_genomes_total              ~ sprintf("Unique (1 genome)"),
      TRUE                         ~ sprintf("Shared (%d genomes)", n_genomes)
    )
  )

trans_cluster_summary<-trans_cluster_summary%>%
  filter(cluster_size>50)

dup_cluster_summary <- all_dup_clustered %>%
  group_by(cluster_id) %>%
  summarise(
    ref_chr       = ref_chr[1],
    cluster_start = min(ref_start),
    cluster_end   = max(ref_end),
    cluster_size  = cluster_end - cluster_start + 1L,
    n_genomes     = n_distinct(genome),
    genomes_list  = paste(sort(unique(genome)), collapse = ","),
    ann_types     = paste(sort(unique(ann_type)), collapse = ","),
    mean_inv_size = mean(inv_size, na.rm = TRUE),
    max_inv_size  = max(inv_size, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  mutate(
    frequency = case_when(
      n_genomes == 1 | n_genomes == n_genomes_total              ~ sprintf("Unique (1 genome)"),
      TRUE                         ~ sprintf("Shared (%d genomes)", n_genomes)
    )
  )

dup_cluster_summary<-dup_cluster_summary%>%
  filter(cluster_size>50)

  # Figures -----------------------------------------------------------------

## FIGURE 3b inversions across chromosomes colored by shared number of genomes ##

frequency_levels <- c(
  "Unique (1 genome)",
  paste0("Shared (", 2:6, " genomes)")
)
frequency_levels <- intersect(frequency_levels, unique(cluster_summary$frequency))

shared_clusters <- cluster_summary %>%
  dplyr::filter(n_genomes > 1 | n_genomes < 7 ) %>%
  dplyr::mutate(
    fill_val = n_genomes,
    alpha_val = 0.9
  )

unique_clusters <- cluster_summary %>%
  dplyr::filter(n_genomes == 1 | n_genomes ==7) %>%
  dplyr::mutate(
    fill_val  = NA_real_,
    alpha_val = 0.4
  )

chr_order <- unique(cluster_summary$ref_chr) %>%
  .[order(as.numeric(gsub("[^0-9]", "", .)),
          na.last = TRUE,
          method  = "radix")]

chr_lengths <- cluster_summary %>%
  dplyr::group_by(ref_chr) %>%
  dplyr::summarise(chr_len = max(cluster_end), .groups = "drop") %>%
  dplyr::mutate(ref_chr = factor(ref_chr, levels = chr_order)) %>%
  dplyr::arrange(ref_chr)%>%
  dplyr::mutate(ref_chr = gsub("^Chr", "", ref_chr))


chr_layout <- chr_lengths %>%
  dplyr::arrange(as.numeric(ref_chr)) %>%
  dplyr::mutate(y = rev(row_number()))

SEG_Y    <- -0.45
SEG_YEND <-  0.45

shared_clusters <- shared_clusters %>%
  dplyr::mutate(ref_chr = factor(ref_chr, levels = chr_order))%>%
  dplyr::mutate(ref_chr = as.character(gsub("^Chr", "", ref_chr))) %>%
  dplyr::left_join(chr_layout %>% dplyr::select(ref_chr, y), 
            by = c("ref_chr")) 

unique_clusters <- unique_clusters %>%
  dplyr::mutate(ref_chr = factor(ref_chr, levels = chr_order))%>%
  dplyr::mutate(ref_chr = as.character(gsub("^Chr", "", ref_chr))) %>%
  dplyr::left_join(chr_layout %>% dplyr::select(ref_chr, y), 
            by = c("ref_chr")) 

shared_min <- 1
shared_max <- 7

genome_colors <- c(
  "1" = "#E69F00",
  "2" = "#56B4E9",
  "3" = "#009E73",
  "4" = "#CC79A7",
  "5" = "#009E73",
  "6" = "#56B4E9",
  "7" = "#E69F00"
)

shared_clusters <- shared_clusters %>%
  mutate(n_genomes = as.character(n_genomes))

p4 <- ggplot() +
  # Chromosome background blocks
  geom_rect(
    data = chr_layout,
    aes(xmin = 0, xmax = chr_len,
        ymin = y + SEG_Y, ymax = y + SEG_YEND),
    fill = "grey92", colour = "grey70", linewidth = 0.3
  ) +
  # Unique inversions
  geom_rect(
    data = unique_clusters,
    aes(xmin = cluster_start, xmax = cluster_end,
        ymin = y + SEG_Y,     ymax = y + SEG_YEND),
    fill  = "grey75",
    color = NA,
    alpha = 0.5
  ) +
  # Shared inversions
  geom_rect(
    data = shared_clusters,
    aes(xmin = cluster_start, xmax = cluster_end,
        ymin = y + SEG_Y,     ymax = y + SEG_YEND,
        fill = n_genomes),
    color = NA,
    alpha = 0.92
  ) +
  scale_fill_manual(
    name   = "# genomes\nsharing",
    values = genome_colors,
    labels = c("1", "2", "3", "4", "5", "6", "7")
  ) +
  scale_x_continuous(
    name   = "Reference position (Mb)",
    labels = function(x) comma(x / 1e6, accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = chr_layout$y,
    labels = chr_layout$ref_chr
  ) +
  labs(title = "Comparative Alignment Inversions >50 bp across genomes", y = NULL) +
  theme_cowplot(11) +
  theme(
    strip.text         = element_text(face = "bold", size = 9),
    panel.spacing      = unit(0.4, "lines"),
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    legend.position    = "right",
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
print(p4)

ggsave("fig4_shared_inversions_bychromosome_min50bp.svg",
       p4, width = 12, height = 10)



# chromosome lebel by population
all_clusters <- cluster_summary %>%
  dplyr::mutate(ref_chr = factor(ref_chr, levels = chr_order))%>%
  dplyr::mutate(ref_chr = as.character(gsub("^Chr", "", ref_chr))) %>%
  dplyr::left_join(chr_layout %>% dplyr::select(ref_chr, y), 
            by = c("ref_chr")) 

p6 <- ggplot() +
  # Chromosome background blocks
  geom_rect(
    data = chr_layout,
    aes(xmin = 0, xmax = chr_len,
        ymin = y + SEG_Y, ymax = y + SEG_YEND),
    fill = "grey92", colour = "grey70", linewidth = 0.3
  ) +
  geom_rect(
    data = all_clusters,
    aes(xmin = cluster_start, xmax = cluster_end,
        ymin = y + SEG_Y,     ymax = y + SEG_YEND,
        fill = Population),
    color = NA,
    alpha = 0.92
  ) +
  scale_fill_manual(values = c("#E69F00","#61D04F","#2297E6"), name = "Population") +
  scale_x_continuous(
    name   = "Reference position (Mb)",
    labels = function(x) comma(x / 1e6, accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = chr_layout$y,
    labels = chr_layout$ref_chr
  ) +
  labs(title = "Comparative Alignment Inversions >50 bp across genomes", y = NULL) +
  theme_cowplot(11) +
  theme(
    strip.text         = element_text(face = "bold", size = 9),
    panel.spacing      = unit(0.4, "lines"),
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    legend.position    = "right",
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
print(p6)


ggsave("fig_compalign_bypop_min50bp.svg",
       p6, width = 12, height = 10)

# histogram of population frequency

p6b <- ggplot(all_clusters,
              aes(x=Population,fill = Population)) +
  geom_bar() +
  scale_y_continuous(name = "Number of inversion clusters") +
  scale_fill_manual(values = c("#E69F00","#61D04F","#2297E6"), name = "Population") +
  labs(title = "Inversion frequency by population") +
  theme_cowplot(11) +
  theme(legend.position  = "none",
        strip.background = element_rect(fill = "grey92"),
        plot.subtitle    = element_text(size = 9, color = "grey40"))

print(p6b)
ggsave("fig_compalign_freqeuncyofpop.svg",
       p6b, width = 8, height = 8)


### COMPARE METHODS ###---------------------------------------------------------
#Compare my three inversion calling methods:
#1) method=pangenome; df name=df_parsed_sizes_oneofeach; columns= Chromosome, start_pos, end_pos
#2) method=localPCA; df name=merged_runs; columns=Chromosome, start_pos, end_pos
#3) method=comparativealignment; df name=compalign; columns=ref_chr, cluster_start, cluster_end. 

  # Standarize Dataframes --------------------------------------------------------

pg <- df_sizes_annotated_oneofeach %>%
  transmute(
    method    = "pangenome",
    chr       = as.character(Chromosome),
    start_pos = start,
    end_pos   = stop,
    size      = size
  )

lp <- all_gbs_filtered_merged80  %>%
  mutate(Chromosome = as.character(gsub("^Chr", "", chrom))) %>%
  transmute(
    method    = "localPCA",
    chr       = as.character(Chromosome),
    start_pos = pmin(as.integer(start_pos), as.integer(end_pos)),
    end_pos   = pmax(as.integer(start_pos), as.integer(end_pos)),
    size      = abs(end_pos - start_pos)
  ) 

ca <- cluster_summary %>%
  mutate(Chromosome = as.character(gsub("^Chr", "", ref_chr))) %>%
  transmute(
    method    = "comparativealignment",
    chr       = Chromosome,
    start_pos = pmin(as.integer(cluster_start), as.integer(cluster_end)),
    end_pos   = pmax(as.integer(cluster_start), as.integer(cluster_end)),
    size      = abs(end_pos - start_pos)
  )

pg2 <- invpg %>%
  transmute(
    method    = "pangenome2",
    chr       = as.character(Chromosome),
    start_pos = POS,
    end_pos   = pos_end,
    size      = REF_bp
  )

message(sprintf("  pangenome            : %d inversions", nrow(pg))) #51 inversions
message(sprintf("  localPCA             : %d inversions", nrow(lp))) #32 inversions
message(sprintf("  comparativealignment : %d inversions", nrow(ca))) #4904 inversions
message(sprintf("  pangenome2            : %d inversions", nrow(pg2))) #403 inversions


all_inv <- bind_rows(
  pg %>% dplyr::select(method, chr, start_pos, end_pos, size),
  lp %>% dplyr::select(method, chr, start_pos, end_pos, size),
  ca %>% dplyr::select(method, chr, start_pos, end_pos, size)
) %>% dplyr::mutate(global_idx = row_number())

all_inv_summary<-all_inv  %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    mean   = mean(size, na.rm = TRUE),
    median = median(size, na.rm = TRUE),
    min    = min(size, na.rm = TRUE),
    max    = max(size, na.rm = TRUE),
    n      = n()
  )

tapply(all_inv$size, all_inv$method, shapiro.test)
#size is not normally distributed

kruskal.test(size ~ method, data = all_inv) #sig
dunn.test(all_inv$size, all_inv$method, method = "bh")

#redo with new pangenome inversion annotation
all_inv2 <- bind_rows(
  pg %>% select(method, chr, start_pos, end_pos, size),
  pg2 %>% select(method, chr, start_pos, end_pos, size),
  lp %>% select(method, chr, start_pos, end_pos, size),
  ca %>% select(method, chr, start_pos, end_pos, size)
) %>% mutate(global_idx = row_number())

all_inv_summary2<-all_inv2  %>%
  group_by(method) %>%
  summarise(
    mean   = mean(size, na.rm = TRUE),
    median = median(size, na.rm = TRUE),
    min    = min(size, na.rm = TRUE),
    max    = max(size, na.rm = TRUE),
    n      = n()
  )


  # Identify Reciprocal Overlap, Any ----------------------------------------
cluster_methods_any <- function(all_inv) {
  
  df <- all_inv %>% mutate(cluster_id = NA_integer_)
  
  gr <- GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start_pos, end = df$end_pos)
  )
  
  next_cluster <- 1L
  unassigned   <- which(is.na(df$cluster_id))
  
  while (length(unassigned) > 0) {
    
    anchor_idx    <- unassigned[which.max(df$size[unassigned])]
    anchor        <- df[anchor_idx, ]
    
    anchor_gr     <- GRanges(
      seqnames = anchor$chr,
      ranges   = IRanges(start = anchor$start_pos, end = anchor$end_pos)
    )
    candidates_gr <- gr[unassigned]
    hits          <- findOverlaps(anchor_gr, candidates_gr, ignore.strand = TRUE)
    
    if (length(hits) == 0) {
      df$cluster_id[anchor_idx] <- next_cluster
      next_cluster <- next_cluster + 1L
      unassigned   <- which(is.na(df$cluster_id))
      next
    }
    
    candidate_local_idx  <- subjectHits(hits)
    candidate_global_idx <- unassigned[candidate_local_idx]
    cands                <- df[candidate_global_idx, ]
    
    # Only requirement: different method from anchor
    passes  <- cands$method != anchor$method
    members <- c(anchor_idx, candidate_global_idx[passes])
    
    df$cluster_id[members] <- next_cluster
    next_cluster <- next_cluster + 1L
    unassigned   <- which(is.na(df$cluster_id))
  }
  
  df
}

all_inv_anyoverlap <- cluster_methods_any(all_inv)

all_inv_anyoverlap %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(n_methods = n_distinct(method)) %>%
  dplyr::count(n_methods)

all_inv_anyoverlap %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::filter(n_distinct(method) > 1) %>%
  nrow() #1714

all_inv_anyoverlap %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(
    has_pg = "pangenome"            %in% method,
    has_lp = "localPCA"             %in% method,
    has_ca = "comparativealignment" %in% method,
    .groups = "drop"
  ) %>%
  dplyr::count(has_pg, has_lp, has_ca)

venn_list <- list(
  pangenome            = all_inv_anyoverlap %>% filter(method=="pangenome")%>% pull(cluster_id),
  localPCA             = all_inv_anyoverlap %>% filter(method=="localPCA")%>% pull(cluster_id),
  comparativealignment = all_inv_anyoverlap %>% filter(method=="comparativealignment")%>% pull(cluster_id)
)


p_venn <- ggvenn(
  venn_list,
  fill_color   = c("#E69F00", "#56B4E9", "#CC79A7"),
  fill_alpha   = 0.4,
  stroke_size  = 0.6,
  text_size    = 4,
  set_name_size = 4
) +
  labs(
    title    = "Overlap of inversion calls across SV methods",
    subtitle = "Any Overlap") +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))
print(p_venn)
ggsave("fig_venn_sv_methods_anyoverlap.svg", p_venn, width = 7, height = 6)

#redo with new pg inversion dataset

all_inv_anyoverlap2 <- cluster_methods_any(all_inv2)

all_inv_anyoverlap2 %>%
  group_by(cluster_id) %>%
  summarise(n_methods = n_distinct(method)) %>%
  count(n_methods)

all_inv_anyoverlap2 %>%
  group_by(cluster_id) %>%
  filter(n_distinct(method) > 1) %>%
  nrow() #901

all_inv_anyoverlap2 %>%
  group_by(cluster_id) %>%
  summarise(
    has_pg = "pangenome"            %in% method,
    has_pg2 = "pangenome2"            %in% method,
    has_lp = "localPCA"             %in% method,
    has_ca = "comparativealignment" %in% method,
    .groups = "drop"
  ) %>%
  count(has_pg,  has_pg2, has_lp, has_ca)

venn_list2 <- list(
  pangenome            = all_inv_anyoverlap2 %>% filter(method=="pangenome")%>% pull(cluster_id),
  pangenome2            = all_inv_anyoverlap2 %>% filter(method=="pangenome2")%>% pull(cluster_id),
  localPCA             = all_inv_anyoverlap2 %>% filter(method=="localPCA")%>% pull(cluster_id),
  comparativealignment = all_inv_anyoverlap2 %>% filter(method=="comparativealignment")%>% pull(cluster_id)
)


p_venn2 <- ggvenn(
  venn_list2,
  fill_color   = c("#E69F00", "#009E73","#56B4E9", "#CC79A7"),
  fill_alpha   = 0.4,
  stroke_size  = 0.6,
  text_size    = 4,
  set_name_size = 4
) +
  labs(
    title    = "Overlap of inversion calls across SV methods",
    subtitle = "Any Overlap") +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))
print(p_venn2)
ggsave("fig_venn_sv_methods_anyoverlap_withnewpginversion.svg", p_venn2, width = 7, height = 6)

  # Identify reciprocal Overlap, 2/3 and 80% --------------------------------

# Overlap defined as reciprocal overlap >= RO_THRESHOLD
# RO(A,B) = overlap_length / min(size_A, size_B)
RO_THRESHOLD <- 0.80

has_ro_match <- function(focal, other, threshold = RO_THRESHOLD, size_ratio_min = 2/3) {
  if (nrow(focal) == 0 || nrow(other) == 0) return(rep(FALSE, nrow(focal)))
  
  focal_gr <- GRanges(
    seqnames = focal$chr,
    ranges   = IRanges(start = focal$start_pos, end = focal$end_pos)
  )
  other_gr <- GRanges(
    seqnames = other$chr,
    ranges   = IRanges(start = other$start_pos, end = other$end_pos)
  )
  
  hits <- findOverlaps(focal_gr, other_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(rep(FALSE, nrow(focal)))
  
  fi <- queryHits(hits)
  oi <- subjectHits(hits)
  
  ov_len     <- pmax(0L,
                     pmin(focal$end_pos[fi],   other$end_pos[oi]) -
                       pmax(focal$start_pos[fi], other$start_pos[oi]) + 1L)
  min_size   <- pmin(focal$size[fi], other$size[oi])
  max_size   <- pmax(focal$size[fi], other$size[oi])
  size_ratio <- min_size / max_size
  ro         <- ifelse(size_ratio >= size_ratio_min, ov_len / min_size, 0)
  
  matched <- rep(FALSE, nrow(focal))
  for (i in seq_along(fi)) {
    if (ro[i] >= threshold) matched[fi[i]] <- TRUE
  }
  matched
}

## compute overlap ##
pg <- pg %>% mutate(
  in_localPCA  = has_ro_match(pg, lp),
  in_compalign = has_ro_match(pg, ca)
)

lp <- lp %>% mutate(
  in_pangenome = has_ro_match(lp, pg),
  in_compalign = has_ro_match(lp, ca)
)

ca <- ca %>% mutate(
  in_pangenome = has_ro_match(ca, pg),
  in_localPCA  = has_ro_match(ca, lp)
)

#compute overlap with new pangnoem inversions
pg3 <- pg %>% mutate(
  in_localPCA  = has_ro_match(pg, lp),
  in_compalign = has_ro_match(pg, ca),
  in_pangenome2 = has_ro_match(pg, pg2)
)

lp2 <- lp %>% mutate(
  in_pangenome = has_ro_match(lp, pg),
  in_compalign = has_ro_match(lp, ca),
  in_pangenome2 = has_ro_match(lp, pg2)
)

ca2 <- ca %>% mutate(
  in_pangenome = has_ro_match(ca, pg),
  in_localPCA  = has_ro_match(ca, lp),
  in_pangenome2 = has_ro_match(ca, pg2)
)

pg2 <- pg2 %>% mutate(
  in_localPCA  = has_ro_match(pg2, lp),
  in_compalign = has_ro_match(pg2, ca),
  in_pangenome2 = has_ro_match(pg2, pg)
)
### CLUSTER ACROSS METHODS ###
cluster_methods_fn <- function(all_inv, ideal_ro = 0.80, size_ratio_min = 2/3) {
  
  df <- all_inv %>% mutate(cluster_id = NA_integer_)
  
  gr <- GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start_pos, end = df$end_pos)
  )
  
  next_cluster <- 1L
  unassigned   <- which(is.na(df$cluster_id))
  
  while (length(unassigned) > 0) {
    
    anchor_idx <- unassigned[which.max(df$size[unassigned])]
    anchor     <- df[anchor_idx, ]
    
    anchor_gr     <- GRanges(
      seqnames = anchor$chr,
      ranges   = IRanges(start = anchor$start_pos, end = anchor$end_pos)
    )
    candidates_gr <- gr[unassigned]
    hits <- findOverlaps(anchor_gr, candidates_gr, ignore.strand = TRUE)
    
    if (length(hits) == 0) {
      df$cluster_id[anchor_idx] <- next_cluster
      next_cluster <- next_cluster + 1L
      unassigned   <- which(is.na(df$cluster_id))
      next
    }
    
    candidate_local_idx  <- subjectHits(hits)
    candidate_global_idx <- unassigned[candidate_local_idx]
    
    cands <- df[candidate_global_idx, ]
    
    ov_len <- pmax(0L,
                   pmin(anchor$end_pos,   cands$end_pos) -
                     pmax(anchor$start_pos, cands$start_pos) + 1L)
    
    min_size   <- pmin(anchor$size, cands$size)
    max_size   <- pmax(anchor$size, cands$size)
    size_ratio <- min_size / max_size
    ro         <- ifelse(size_ratio >= size_ratio_min, ov_len / min_size, 0)
    
    # Only allow cross-method joins
    diff_method <- cands$method != anchor$method
    passes      <- diff_method & size_ratio >= size_ratio_min & ro >= ideal_ro
    
    members <- c(anchor_idx, candidate_global_idx[passes])
    
    # All-pairs check 
    repeat {
      member_df <- df[members, ]
      n         <- length(members)
      if (n <= 1) break
      
      all_pass <- rep(TRUE, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          if (i == j) next
          if (member_df$method[i] == member_df$method[j]) next
          ov    <- max(0L,
                       min(member_df$end_pos[i],   member_df$end_pos[j]) -
                         max(member_df$start_pos[i], member_df$start_pos[j]) + 1L)
          min_s <- min(member_df$size[i], member_df$size[j])
          max_s <- max(member_df$size[i], member_df$size[j])
          ro_ij <- ifelse(min_s / max_s >= size_ratio_min, ov / min_s, 0)
          if (ro_ij < ideal_ro) {
            all_pass[i] <- FALSE
            break
          }
        }
      }
      
      if (all(all_pass)) break
      members <- members[all_pass]
    }
    
    # Final size ratio check
    if (length(members) == 0) {
      df$cluster_id[anchor_idx] <- next_cluster
      next_cluster <- next_cluster + 1L
      unassigned   <- which(is.na(df$cluster_id))
      next
    }
    
    member_sizes       <- df$size[members]
    if (any(is.na(member_sizes)) || max(member_sizes) == 0) {
      df$cluster_id[members] <- next_cluster
      next_cluster <- next_cluster + 1L
      unassigned   <- which(is.na(df$cluster_id))
      next
    }
    
    cluster_size_ratio <- min(member_sizes) / max(member_sizes)
    if (cluster_size_ratio < size_ratio_min) {
      size_order     <- order(member_sizes, decreasing = TRUE)
      members_sorted <- members[size_order]
      sizes_sorted   <- member_sizes[size_order]
      keep           <- sizes_sorted >= sizes_sorted[1] * size_ratio_min
      members        <- members_sorted[keep]
    }
    
    df$cluster_id[members] <- next_cluster
    next_cluster <- next_cluster + 1L
    unassigned   <- which(is.na(df$cluster_id))
  }
  
  df
}

all_inv_RO <- cluster_methods_fn(all_inv, ideal_ro = 0.80, size_ratio_min = 2/3)

all_inv_RO2 <- cluster_methods_fn(all_inv2, ideal_ro = 0.80, size_ratio_min = 2/3)


### venn diagram ###
cluster_methods_RO <- all_inv_RO %>%
  group_by(cluster_id) %>%
  summarise(
    in_pangenome = "pangenome"            %in% method,
    in_localPCA  = "localPCA"             %in% method,
    in_compalign = "comparativealignment" %in% method,
    n_methods    = n_distinct(method),
    max_size     = max(size),
    .groups      = "drop"
  )



venn_list <- list(
  pangenome            = cluster_methods_RO %>% filter(in_pangenome)%>% pull(cluster_id),
  localPCA             = cluster_methods_RO %>% filter(in_localPCA)%>% pull(cluster_id),
  comparativealignment = cluster_methods_RO %>% filter(in_compalign)%>% pull(cluster_id)
)

p_venn_RO <- ggvenn(
  venn_list,
  fill_color   = c("#E69F00", "#56B4E9", "#CC79A7"),
  fill_alpha   = 0.4,
  stroke_size  = 0.6,
  text_size    = 4,
  set_name_size = 4
) +
  labs(
    title    = "Overlap of inversion calls across SV methods",
    subtitle = sprintf("Reciprocal overlap >= %.0f%%", RO_THRESHOLD * 100)
  ) +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))
print(p_venn_RO)
ggsave("fig_venn_sv_methods_RO80.svg", p_venn_RO, width = 7, height = 6)

### venn diagram with new inversions from pangenome###
cluster_methods_RO2 <- all_inv_RO2 %>%
  group_by(cluster_id) %>%
  summarise(
    in_pangenome = "pangenome"            %in% method,
    in_pangenome2 = "pangenome2"            %in% method,
    in_localPCA  = "localPCA"             %in% method,
    in_compalign = "comparativealignment" %in% method,
    n_methods    = n_distinct(method),
    max_size     = max(size),
    .groups      = "drop"
  )



venn_list2 <- list(
  pangenome            = cluster_methods_RO2 %>% filter(in_pangenome)%>% pull(cluster_id),
  pangenome2            = cluster_methods_RO2 %>% filter(in_pangenome2)%>% pull(cluster_id),
  localPCA             = cluster_methods_RO2 %>% filter(in_localPCA)%>% pull(cluster_id),
  comparativealignment = cluster_methods_RO2 %>% filter(in_compalign)%>% pull(cluster_id)
)

p_venn_RO2 <- ggvenn(
  venn_list2,
  fill_color   = c("#E69F00", "#009E73","#56B4E9", "#CC79A7"),
  fill_alpha   = 0.4,
  stroke_size  = 0.6,
  text_size    = 4,
  set_name_size = 4
) +
  labs(
    title    = "Overlap of inversion calls across SV methods",
    subtitle = sprintf("Reciprocal overlap >= %.0f%%", RO_THRESHOLD * 100)
  ) +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))
print(p_venn_RO2)
ggsave("fig_venn_sv_methods_RO80_withnewpginversions.svg", p_venn_RO2, width = 7, height = 6)


  # Size/Method analysis, 2/3 and 80% ----------------------------------------------------

# Size histogram faceted by method
p_hist <- ggplot(all_inv_RO, aes(x = size / 1000, fill = method)) +
  geom_histogram(bins = 50, color = "white", linewidth = 0.2, alpha = 0.85) +
  scale_x_log10(
    name   = "Inversion size (kb, log scale)",
    labels = comma_format(accuracy = 0.1)
  ) +
  scale_y_continuous(name = "Count") +
  scale_fill_manual(
    values = c("pangenome"            = "#E69F00",
               "localPCA"             = "#56B4E9",
               "comparativealignment" = "#CC79A7"),
    guide  = "none"
  ) +
  facet_wrap(~ method, ncol = 1, scales = "free_y") +
  labs(title    = "Inversion size distribution per method") +
  theme_cowplot(11) +
  theme(strip.background = element_rect(fill = "grey92"),
        plot.subtitle    = element_text(size = 9, color = "grey40"))
print(p_hist)
ggsave("fig_size_hist_by_method.pdf", p_hist, width = 8, height = 8)


p_hist2 <- ggplot(all_inv_RO2, aes(x = size / 1000, fill = method)) +
  geom_histogram(bins = 50, color = "white", linewidth = 0.2, alpha = 0.85) +
  scale_x_log10(
    name   = "Inversion size (kb, log scale)",
    labels = comma_format(accuracy = 0.1)
  ) +
  scale_y_continuous(name = "Count") +
  scale_fill_manual(
    values = c("pangenome"            = "#E69F00",
               "pangenome2"           = "#009E73",
               "localPCA"             = "#56B4E9",
               "comparativealignment" = "#CC79A7"),
    guide  = "none"
  ) +
  facet_wrap(~ method, ncol = 1, scales = "free_y") +
  labs(title    = "Inversion size distribution per method") +
  theme_cowplot(11) +
  theme(strip.background = element_rect(fill = "grey92"),
        plot.subtitle    = element_text(size = 9, color = "grey40"))
print(p_hist2)
ggsave("fig_size_hist_by_method_newinverionsmethod.pdf", p_hist2, width = 8, height = 8)


#compare size versus shared method
all_inv_RO <- all_inv_RO  %>%
  group_by(cluster_id) %>%
  mutate(shared_methods = n_distinct(method)) %>%
  ungroup()

all_inv_RO$shared_methods<-as.factor(all_inv_RO$shared_methods)
# Fig: size vs overlap group (across all methods combined)
p_size_overlap_RO <- ggplot(all_inv_RO,
                         aes(x = shared_methods, y = size / 1000)) +
  geom_boxplot(alpha = 0.35, linewidth = 0.7) +
  geom_jitter(width = 0.18, alpha = 0.5, size = 1.2, shape = 16) +
  scale_y_log10(
    name   = "Inversion size (kb, log scale)",
    labels = comma_format(accuracy = 0.1)
  ) +
  scale_x_discrete(name = NULL) +
  labs(
    title    = "Inversion size vs. number of methods detecting it",
    subtitle = "Reciprocal overlap >= 80%") +
  theme_cowplot(12) +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))+
  stat_compare_means(label = "p.signif",comparisons = list(c("1", "2")))
print(p_size_overlap_RO)
ggsave("fig_size_vs_overlap_RO.pdf", p_size_overlap_RO, width = 8, height = 6)
ggsave("fig_size_vs_overlap_RO.svg", p_size_overlap_RO, width = 8, height = 6)

# Fig: per-method detected vs missed
all_inv_RO$shared_methods<-as.numeric(all_inv_RO$shared_methods)
p_detect <- all_inv_RO %>%
  mutate(detected = if_else(shared_methods > 1, "Detected by\nanother method",
                            "Not detected by\nany other method")) %>%
  ggplot(aes(x = detected, y = size / 1000,
             fill = detected, color = detected)) +
  geom_violin(alpha = 0.35, linewidth = 0.7,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.0, shape = 16) +
  scale_y_log10(
    name   = "Inversion size (kb, log scale)",
    labels = comma_format(accuracy = 0.1)
  ) +
  scale_x_discrete(name = NULL) +
  scale_fill_manual(
    values = c("Detected by\nanother method"      = "#009E73",
               "Not detected by\nany other method" = "#D55E00"),
    guide  = "none"
  ) +
  scale_color_manual(
    values = c("Detected by\nanother method"      = "#009E73",
               "Not detected by\nany other method" = "#D55E00"),
    guide  = "none"
  ) +
  facet_wrap(~ method, ncol = 3) +
  labs(
    title    = "Does inversion size predict whether a method is corroborated?",
    subtitle = sprintf("RO >= %.0f%%  |  each point = one raw inversion call",
                       RO_THRESHOLD * 100)
  ) +
  theme_cowplot(11) +
  theme(strip.background = element_rect(fill = "grey92"),
        plot.subtitle    = element_text(size = 9, color = "grey40"))
print(p_detect)
ggsave("fig_detection_by_size_RO.pdf", p_detect, width = 11, height = 5)

# Stats
wt_z <-wilcox.test(size ~ shared_methods, data = all_inv_RO)
library(rstatix)
all_inv_RO %>% wilcox_effsize(size ~ shared_methods)
all_inv_RO %>%
  group_by(shared_methods) %>%
  summarise(
    median_size = median(size, na.rm = TRUE),
    IQR_size = IQR(size, na.rm = TRUE),
    n = n()
  )
# Size/Method analysis, any overlap ----------------------------------------------------


#compare size versus shared method
all_inv_anyoverlap <- all_inv_anyoverlap  %>%
  group_by(cluster_id) %>%
  mutate(shared_methods = n_distinct(method)) %>%
  ungroup()

all_inv_anyoverlap$shared_methods<-as.factor(all_inv_anyoverlap$shared_methods)
# Fig: size vs overlap group (across all methods combined)
p_size_overlap_anyoverlap <- ggplot(all_inv_anyoverlap,
                            aes(x = shared_methods, y = size / 1000)) +
  geom_boxplot(alpha = 0.35, linewidth = 0.7) +
  geom_jitter(width = 0.18, alpha = 0.5, size = 1.2, shape = 16) +
  scale_y_log10(
    name   = "Inversion size (kb, log scale)",
    labels = comma_format(accuracy = 0.1)
  ) +
  scale_x_discrete(name = NULL) +
  labs(
    title    = "Inversion size vs. number of methods detecting it",
    subtitle = "Any overlap") +
  theme_cowplot(12) +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))+
  stat_compare_means(label = "p.signif",comparisons = list(c("1", "2"), c("2", "3")))
print(p_size_overlap_anyoverlap)
ggsave("fig_size_vs_overlap_anyoverlap.pdf", p_size_overlap_anyoverlap, width = 8, height = 6)
ggsave("fig_size_vs_overlap_anyoverlap.svg", p_size_overlap_anyoverlap, width = 8, height = 6)

# Fig: per-method detected vs missed
all_inv_anyoverlap$shared_methods<-as.numeric(all_inv_anyoverlap$shared_methods)
p_detect <- all_inv_anyoverlap %>%
  mutate(detected = if_else(shared_methods > 1, "Detected by\nanother method",
                            "Not detected by\nany other method")) %>%
  ggplot(aes(x = detected, y = size / 1000,
             fill = detected, color = detected)) +
  geom_violin(alpha = 0.35, linewidth = 0.7,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.0, shape = 16) +
  scale_y_log10(
    name   = "Inversion size (kb, log scale)",
    labels = comma_format(accuracy = 0.1)
  ) +
  scale_x_discrete(name = NULL) +
  scale_fill_manual(
    values = c("Detected by\nanother method"      = "#009E73",
               "Not detected by\nany other method" = "#D55E00"),
    guide  = "none"
  ) +
  scale_color_manual(
    values = c("Detected by\nanother method"      = "#009E73",
               "Not detected by\nany other method" = "#D55E00"),
    guide  = "none"
  ) +
  facet_wrap(~ method, ncol = 3) +
  labs(
    title    = "Does inversion size predict whether a method is corroborated?",
    subtitle = sprintf("RO >= %.0f%%  |  each point = one raw inversion call",
                       RO_THRESHOLD * 100)
  ) +
  theme_cowplot(11) +
  theme(strip.background = element_rect(fill = "grey92"),
        plot.subtitle    = element_text(size = 9, color = "grey40"))
print(p_detect)
ggsave("fig_detection_by_size_RO.pdf", p_detect, width = 11, height = 5)

# Stats
kruskal.test(size ~ shared_methods, data = all_inv_RO) %>% print()


  # Proportion Covered ------------------------------------------------------
genome_size <- sum(scaffold_info$scaf_length) #1226560475

# For each method, sum inverted bases and find proportion of total genome
method_coverage_summary <- all_inv %>%
  group_by(method) %>%
  group_map(~ {
    gr <- GRanges(seqnames = .x$chr,
                  ranges   = IRanges(start = .x$start_pos, end = .x$end_pos))
    covered <- sum(width(reduce(gr)))
    tibble(method = .y$method,
           bases_covered = covered,
           prop_genome   = covered / genome_size)
  }) %>%
  bind_rows()

method_coverage_summary2 <- all_inv2 %>%
  group_by(method) %>%
  group_map(~ {
    gr <- GRanges(seqnames = .x$chr,
                  ranges   = IRanges(start = .x$start_pos, end = .x$end_pos))
    covered <- sum(width(reduce(gr)))
    tibble(method = .y$method,
           bases_covered = covered,
           prop_genome   = covered / genome_size)
  }) %>%
  bind_rows()

#test observed vs null for proportion inverted across genome
genome_length <- sum(scaffold_info$scaf_length)  # 1226560475
n1 <- 285363834      # bp called by comparative alignment
n2 <- 316242518      # bp called by local PCA
n3 <- 78499255       # bp called by pangenome
observed_overlap_3way <- 0.041285599 * genome_length    
cp_lp_overlap<- 0.111157468 * genome_length
pg_lp_overlap<- 0.010224070 * genome_length
cp_pg_overlap<- 0.006172773 * genome_length
  
# proportions covered by each method
p1 <- n1 / genome_length #cp
p2 <- n2 / genome_length #lp
p3 <- n3 / genome_length #pg

# null probability of 3-way overlap at any given bp, assuming independence
p_null_3way <- p1 * p2 * p3
expected_overlap_3way <- p_null_3way * genome_length

p_null_cp_lp <- p1 * p2
expected_overlap_cp_lp <- p_null_cp_lp * genome_length

p_null_pg_lp <- p2 * p3
expected_overlap_pg_lp <- p_null_pg_lp * genome_length

p_null_cp_pg <- p1 * p3
expected_overlap_cp_pg <- p_null_cp_pg * genome_length

# variance/SD under binomial null 3way
var_null <- genome_length * p_null_3way * (1 - p_null)_3way
sd_null <- sqrt(var_null_3way)

# z-score and p-value 3way
z <- (observed_overlap_3way - expected_overlap_3way) / sd_null
p_value_one_sided <- 1 - pnorm(z)      # is observed MORE than expected?
p_value_two_sided <- 2 * pnorm(-abs(z))

cat("Z-score:", z, "\n")
cat("One-sided p-value:", p_value_one_sided, "\n")

# variance/SD under binomial null 2way
var_null <- genome_length * p_null_cp_pg * (1 - p_null_cp_pg)
sd_null <- sqrt(var_null)

# z-score and p-value 2way
z <- (cp_pg_overlap - expected_overlap_cp_pg) / sd_null
p_value_one_sided <- 1 - pnorm(z)      # is observed MORE than expected?
p_value_two_sided <- 2 * pnorm(-abs(z))

cat("Z-score:", z, "\n")
cat("One-sided p-value:", p_value_one_sided, "\n")

#test observed vs null for number of inversions detected by any overlap
n1 <- 3300     # inversions called by comparative alignment
n2 <- 28      # inversions called by local PCA
n3 <- 41      # inversions called by pangenome
total_inversions<-3316
observed_overlap_3way <- 10    # observed 3-way overlap 
cp_lp_overlap<- 17
pg_lp_overlap<- 0
cp_pg_overlap<- 16

# proportions covered by each method
p1 <- n1 / total_inversions #cp
p2 <- n2 / total_inversions #lp
p3 <- n3 / total_inversions #pg

# null probability of 3-way overlap, assuming independence
p_null_3way <- p1 * p2 * p3
expected_overlap_3way <- p_null_3way * total_inversions

p_null_cp_lp <- p1 * p2
expected_overlap_cp_lp <- p_null_cp_lp * total_inversions

p_null_pg_lp <- p2 * p3
expected_overlap_pg_lp <- p_null_pg_lp * total_inversions

p_null_cp_pg <- p1 * p3
expected_overlap_cp_pg <- p_null_cp_pg * total_inversions

# variance/SD under binomial null 3way
var_null <- total_inversions * p_null_3way * (1 - p_null_3way)
sd_null <- sqrt(var_null)

# z-score and p-value 3way
z <- (observed_overlap_3way - expected_overlap_3way) / sd_null
p_value_one_sided <- 1 - pnorm(z)      # is observed MORE than expected?
p_value_two_sided <- 2 * pnorm(-abs(z))

cat("Z-score:", z, "\n")
cat("One-sided p-value:", p_value_one_sided, "\n")
cat("Two-sided p-value:", p_value_two_sided, "\n")

# variance/SD under binomial null 2way
var_null <- total_inversions * p_null_cp_pg * (1 - p_null_cp_pg)
sd_null <- sqrt(var_null)

# z-score and p-value 2way
z <- (cp_pg_overlap - expected_overlap_cp_pg) / sd_null
p_value_one_sided <- 1 - pnorm(z)      # is observed MORE than expected?
p_value_two_sided <- 2 * pnorm(-abs(z))

cat("Z-score:", z, "\n")
cat("One-sided p-value:", p_value_one_sided, "\n")
cat("Two-sided p-value:", p_value_two_sided, "\n")

#test observed vs null for number of inversions detected by 80% position 2/3 size overlap
n1 <- 4904     # inversions called by comparative alignment
n2 <- 32      # inversions called by local PCA
n3 <- 49      # inversions called by pangenome
total_inversions<-4967
observed_overlap_3way <- 0    # observed 3-way overlap 
cp_lp_overlap<- 5
pg_lp_overlap<- 1
cp_pg_overlap<- 12

# proportions covered by each method
p1 <- n1 / total_inversions #cp
p2 <- n2 / total_inversions #lp
p3 <- n3 / total_inversions #pg

# null probability of 3-way overlap, assuming independence
p_null_3way <- p1 * p2 * p3
expected_overlap_3way <- p_null_3way * total_inversions

p_null_cp_lp <- p1 * p2
expected_overlap_cp_lp <- p_null_cp_lp * total_inversions

p_null_pg_lp <- p2 * p3
expected_overlap_pg_lp <- p_null_pg_lp * total_inversions

p_null_cp_pg <- p1 * p3
expected_overlap_cp_pg <- p_null_cp_pg * total_inversions

# variance/SD under binomial null 3way
var_null <- total_inversions * p_null_3way * (1 - p_null_3way)
sd_null <- sqrt(var_null)

# z-score and p-value 3way
z <- (observed_overlap_3way - expected_overlap_3way) / sd_null
p_value_one_sided <- 1 - pnorm(z)      # is observed MORE than expected?
p_value_two_sided <- 2 * pnorm(-abs(z))

cat("Z-score:", z, "\n")
cat("One-sided p-value:", p_value_one_sided, "\n")
cat("Two-sided p-value:", p_value_two_sided, "\n")

# variance/SD under binomial null 2way
var_null <- total_inversions * p_null_cp_pg * (1 - p_null_cp_pg)
sd_null <- sqrt(var_null)

# z-score and p-value 2way
z <- (cp_pg_overlap - expected_overlap_cp_pg) / sd_null
p_value_one_sided <- 1 - pnorm(z)      # is observed MORE than expected?
p_value_two_sided <- 2 * pnorm(-abs(z))

cat("Z-score:", z, "\n")
cat("One-sided p-value:", p_value_one_sided, "\n")
cat("Two-sided p-value:", p_value_two_sided, "\n")


  # Plot shared proportion across genome ------------------------------------


# Build GRanges per method
methods_list <- unique(all_inv$method)

method_gr <- setNames(
  lapply(methods_list, function(m) {
    sub <- all_inv %>% filter(method == m)
    gr <- GRanges(seqnames = sub$chr,
                  ranges   = IRanges(start = sub$start_pos, end = sub$end_pos))
    reduce(gr)
  }),
  methods_list
)

method_coverage <- map_dbl(method_gr, ~ sum(width(.x)))

map_dfr(names(method_gr), ~ tibble(
     method = .x,
     bases_from_gr       = sum(width(method_gr[[.x]])),
     bases_from_coverage = method_coverage[.x]
   ))

# Build a long-format df with one row per method x category (unique vs shared with each other method)

shared_prop <- map_dfr(names(method_gr), function(m) {
  others <- names(method_gr)[names(method_gr) != m]
  o1 <- others[1]; o2 <- others[2]
  
  full  <- method_gr[[m]]
  
  # Shared with both
  shared_both  <- intersect(intersect(full, method_gr[[o1]]), method_gr[[o2]])
  
  # Shared with o1 only
  shared_o1    <- setdiff(intersect(full, method_gr[[o1]]), method_gr[[o2]])
  
  # Shared with o2 only
  shared_o2    <- setdiff(intersect(full, method_gr[[o2]]), method_gr[[o1]])
  
  # Unique to this method
  unique_only  <- setdiff(setdiff(full, method_gr[[o1]]), method_gr[[o2]])
  
  tibble(
    method   = m,
    category = c("unique",
                 paste0("shared with ", o1, " only"),
                 paste0("shared with ", o2, " only"),
                 "shared with both"),
    prop     = c(sum(width(unique_only)),
                 sum(width(shared_o1)),
                 sum(width(shared_o2)),
                 sum(width(shared_both))) / genome_size
  )
})

# plot
shared_prop <- shared_prop%>%
  mutate(category = factor(category, levels = c(
    "unique",
    paste0("shared with ", names(method_gr)[1], " only"),
    paste0("shared with ", names(method_gr)[2], " only"),
    paste0("shared with ", names(method_gr)[3], " only"),
    "shared with both"
  )))

shared_prop_plot<-ggplot(shared_prop, aes(x = method, y = prop, fill = category)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.5),
    name   = "Proportion of genome"
  ) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = "Method", title = "Inversion coverage by method") +
  theme_classic() +
  theme(legend.position = "right")
print(shared_prop_plot)
ggsave("fig_shared_prop_plot.svg", shared_prop_plot, width = 10, height = 5)


  # plot shared across chromosomes ------------------------------------------

gr_to_df <- function(gr, method, category) {
  if (length(gr) == 0) return(tibble())
  as.data.frame(gr) %>%
    transmute(
      chr      = as.character(seqnames),
      start_pos = start,
      end_pos   = end,
      method   = method,
      category = category
    )
}

SEG_Y    <- -0.4
SEG_YEND <-  0.4
method_height <- 0.25  
method_gap    <- 0.05 

# Assign y offsets for each method
methods_ordered <- unique(all_inv_anyoverlap$method)
method_offsets <- tibble(
  method   = methods_ordered,
  y_method = seq(0, by = method_height + method_gap, length.out = length(methods_ordered))
)

make_chr_plot <- function(chr_num) {
  chr_layout_chr <- chr_layout %>% filter(ref_chr == chr_num)
  
  chr_bg <- chr_layout_chr %>%
    cross_join(method_offsets) %>%
    mutate(
      ymin = y + y_method,
      ymax = ymin + method_height
    )
  
  chr_segments <- map_dfr(names(method_gr), function(m) {
    others <- names(method_gr)[names(method_gr) != m]
    o1 <- others[1]; o2 <- others[2]
    
    full <- method_gr[[m]]
    
    shared_both <- intersect(intersect(full, method_gr[[o1]]), method_gr[[o2]])
    shared_o1   <- setdiff(intersect(full, method_gr[[o1]]), method_gr[[o2]])
    shared_o2   <- setdiff(intersect(full, method_gr[[o2]]), method_gr[[o1]])
    unique_only <- setdiff(setdiff(full, method_gr[[o1]]), method_gr[[o2]])
    
    bind_rows(
      gr_to_df(unique_only, m, "unique"),
      gr_to_df(shared_o1,   m, paste0("shared with ", o1, " only")),
      gr_to_df(shared_o2,   m, paste0("shared with ", o2, " only")),
      gr_to_df(shared_both, m, "shared with both")
    )
  }) %>%
    filter(chr == as.character(chr_num)) %>%
    left_join(method_offsets, by = "method") %>%
    mutate(
      ymin = chr_layout_chr$y + y_method,
      ymax = ymin + method_height
    )
  
  ggplot() +
    geom_rect(
      data = chr_bg,
      aes(xmin = 0, xmax = chr_len,
          ymin = ymin, ymax = ymax),
      fill = "grey92", colour = "grey70", linewidth = 0.3
    ) +
    geom_rect(
      data = chr_segments,
      aes(xmin = start_pos, xmax = end_pos,
          ymin = ymin, ymax = ymax,
          fill = category),
      color = NA, alpha = 0.92
    ) +
    scale_x_continuous(
      name   = "Reference position (Mb)",
      labels = function(x) comma(x / 1e6, accuracy = 1)
    ) +
    scale_y_continuous(
      breaks = method_offsets$y_method + method_height / 2 + chr_layout_chr$y,
      labels = method_offsets$method
    ) +
    scale_fill_manual(
      values = c(
        "unique"            = "#66C2A5",
        "shared with both"  = "#A6D854FF",
        "shared with pangenome only" = "#FC8D62",
        "shared with comparativealignment only" = "#E78AC3",
        "shared with localPCA only" = "#8DA0CB"
      ),
      name = NULL
    )+
    labs(title = paste("Chromosome", chr_num, "inversions by overlap"), y = NULL) +
    theme_cowplot(11) +
    theme(
      legend.position    = "right",
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}


# Generate all plots
chr_plots <- map(1:13, make_chr_plot)
walk(1:13, function(chr_num) {
  ggsave(
    filename = paste0("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/analyses/SharedInversions_chromosomelevel/fig_shared_chr", chr_num, ".svg"),
    plot     = chr_plots[[chr_num]],
    width    = 10,
    height   = 5
  )
})


# Analyze density ---------------------------------------------------------

all_methods <- unique(all_inv_RO$method)

recall_df <- all_inv_RO %>%
  distinct(cluster_id, method) %>%
  mutate(detected = 1) %>%
  complete(cluster_id, method = all_methods, fill = list(detected = 0))

cluster_ranges <- all_inv_RO %>%
  group_by(cluster_id, chr) %>%
  summarise(start = min(start_pos), end = max(end_pos), .groups = "drop")

gr <- GRanges(seqnames = cluster_ranges$chr,
              ranges = IRanges(start = cluster_ranges$start, end = cluster_ranges$end),
              cluster_id = cluster_ranges$cluster_id)

# distance to nearest *other* inversion
nearest <- distanceToNearest(gr)
cluster_ranges$dist_nearest <- NA_real_
cluster_ranges$dist_nearest[queryHits(nearest)] <- mcols(nearest)$distance

# local density: number of other inversions within a window (e.g. 100kb)
window <- 100000
cluster_ranges$density_count <- countOverlaps(gr, gr, maxgap = window) - 1  # -1 to exclude self

#compute local density in inverted bp
cluster_ranges$width <- width(gr)  
overlaps <- findOverlaps(gr, gr, maxgap = window)
overlaps_df <- as.data.frame(overlaps) %>%
  filter(queryHits != subjectHits)  # exclude self-overlap

bases_nearby <- overlaps_df %>%
  mutate(w = cluster_ranges$width[subjectHits]) %>%
  group_by(queryHits) %>%
  summarise(inverted_bases_nearby = sum(w), .groups = "drop")

cluster_ranges$inverted_bases_nearby <- 0
cluster_ranges$inverted_bases_nearby[bases_nearby$queryHits] <- bases_nearby$inverted_bases_nearby

# join both into recall_df
recall_df <- recall_df %>%
  left_join(cluster_ranges %>% select(cluster_id, dist_nearest, density_count, width, inverted_bases_nearby),
            by = "cluster_id")


model <- glm(detected ~ dist_nearest * method, data = recall_df, family = binomial)
summary(model)
logLik(model) #-796.0966 (df=6)

model2 <- glm(detected ~ log1p(density_count) * method, data = recall_df, family = binomial)
summary(model2)
logLik(model2)

model3 <- glm(detected ~ density_count * method, data = recall_df, family = binomial)
summary(model3)
logLik(model3)#-736.5472 (df=6)

model4 <- glm(detected ~ inverted_bases_nearby * method, data = recall_df, family = binomial)
summary(model4)
logLik(model4)# -797.8677 (df=6)

#plot predictions
library(ggeffects)

densityplot<-plot(ggpredict(model3, terms = c("density_count", "method")))+
  labs(
    title = "Predicted detection probability by inversion density",
    x = "Local inversion density (neighboring inversions)",
    y = "Predicted probability of detection",
    color = "Method",
    fill = "Method"
  ) +
  scale_color_manual(values = c("comparativealignment" = "#CC79A7", "localPCA" = "#56B4E9", "pangenome" = "#E69F00")) +
  scale_fill_manual(values = c("comparativealignment" = "#CC79A7", "localPCA" = "#56B4E9", "pangenome" = "#E69F00")) +
  theme_minimal()
print(densityplot)
ggsave(
  filename = paste0("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/analyses/preddetectionbydensity.svg"),
  plot     = densityplot,
  width    = 10,
  height   = 8
)

# analyze gretl output ----------------------------------------------------


setwd("~/Desktop/Gompert_Lab_Research/TimemaSVmethods/Cactus Pangenome/gretl/")
#load datasheets for clipped graph
Clipped_graphstats<-read.delim("gretl_stats_HWY154_REF_4119Hap2.txt", sep="\t")
Clipped_pathstats<-read.delim("gretl_pathstats_HWY154_REF_4119Hap2.txt", sep="\t")

#combine all scaffolds for raw graphs
SCAFFS <- c(
  "Scaffold_1__1_contigs__length_160647932",
  "Scaffold_2__1_contigs__length_157594471",
  "Scaffold_3__2_contigs__length_137956696",
  "Scaffold_4__1_contigs__length_97222829",
  "Scaffold_5__1_contigs__length_83128659",
  "Scaffold_6__1_contigs__length_78844258",
  "Scaffold_7__1_contigs__length_75018798",
  "Scaffold_8__1_contigs__length_71271319",
  "Scaffold_9__2_contigs__length_79556474",
  "Scaffold_10__2_contigs__length_75648701",
  "Scaffold_11__2_contigs__length_80009992",
  "Scaffold_12__1_contigs__length_47609450",
  "Scaffold_13__3_contigs__length_82050896"
)

pathstats_list <- list()

for (scaff in SCAFFS) {
  
  file_path <- paste0("gretl_pathstats_", scaff, ".txt")
  
  pathstats <- read.delim(file_path, sep = "\t")
  
  pathstats_small <- subset(pathstats, grepl("t_crist", Path)) %>%
    mutate(Path = if_else(
      str_detect(Path, "Hap2_t_crist_hwy154_cen4119"),
      paste0(str_remove(Path, "^Hap2_"), "#0"),
      Path
    ))
  
  pathstats_clean <- separate_wider_delim(
    pathstats_small,
    cols = Path,
    delim = "#",
    names = c("Genome", "Hap", "Scaff", "Extra")
  )
  
  pathstats_clean$scaffold_id <- scaff  # track which scaffold this came from
  
  pathstats_list[[scaff]] <- pathstats_clean
}

# combine all scaffolds into one dataframe
pathstats_all <- bind_rows(pathstats_list)

pathstats_summed <- pathstats_all %>%
  group_by(Genome, Hap) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")

pathstats_summed_small<- pathstats_summed%>%
  dplyr::select(c("Genome", "Hap", "Sequence..bp.", "Nodes", "Edges","Inverted.nodes..bp."))%>%
  mutate(PropGenomeInverted=Inverted.nodes..bp./Sequence..bp.)
write.csv(pathstats_summed_small, "gretl_summary.csv", row.names = FALSE)




