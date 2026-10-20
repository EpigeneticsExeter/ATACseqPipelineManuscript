
### Figure 5: Threshold choice analysis
# 5a
# peak number versus Q value curve
### segmented regression to find elbow + 95% CI


library(kneedle)
library(segmented)
library(tidyr)
library(dplyr)
library(GenomicRanges)
library(ggplot2)

options(scipen = 5)
peak_files <- list(
  
  "NeuN_Dataset1" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/q0.05/dataset1/updated_0926/NEUN_5e-2.narrowPeak.filt",
  
  "SOX10_Dataset1" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/q0.05/dataset1/updated_0926/SOX_5e-2.narrowPeak.filt",
  
  "NeuN_Dataset2" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/q0.05/dataset2/neun_q5e-2.narrowPeak.filt",
  
  "SOX10_Dataset2" =
    "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/q0.05/dataset2/sox10_q5e-2.narrowPeak.filt"
)




peak_files

outdir <- "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/Figures/Figure5_Sept"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cell_colours <- c(
  "NeuN"  = "hotpink",
  "SOX10" = "#6B8E23"
)



## thresholds:


exponents <- -seq(2, 100, by = 2)

thresholds <- 10^exponents

############################################################
## STORE RESULTS
############################################################

results_list <- list()

summary_df <- data.frame()
plot_list <- list()



format_q_label <- function(q) {
  
  exponent <- round(log10(q))
  
  mantissa <- q / (10^exponent)
  
  nearest <- c(1, 5)[
    which.min(abs(c(1, 5) - mantissa))
  ]
  
  paste0(nearest, "e", exponent)
}


#plot_data <- list()
plot_df <- data.frame()
for(dataset_name in names(peak_files)) {
  
  cat("\n====================================\n")
  cat("Processing:", dataset_name, "\n")
  cat("====================================\n")
  
  
  ## read in peaks 
  
  
  peaks.gr <- makeGRangesFromDataFrame(
    read.table(peak_files[[dataset_name]]),
    keep.extra.columns = TRUE,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3",
    ignore.strand = TRUE
  )
  
  peaks.df <- as.data.frame(peaks.gr)
  
  colnames(peaks.df) <- c(
    "chr",
    "start",
    "end",
    "width",
    "strand",
    "peakID",
    "score",
    "strand2",
    "signal",
    "pvalue",
    "qvalue",
    "score2"
  )
  
  
  ## Remove sex chromosomes
  peaks.df <- peaks.df %>%
    filter(!chr %in% c("chrX", "chrY"))
  
  ## Convert -log10(q) back to q-values
  peaks.df$qvals <- 10^(-peaks.df$qvalue)
  
  
  
  counts_list <- sapply(
    thresholds,
    function(t) {
      sum(peaks.df$qvals < t, na.rm = TRUE)
    }
  )
  
  ## BUILD DATAFRAME
  
  
  df <- data.frame(
    threshold = thresholds,
    logQ = -log10(thresholds),
    count = counts_list
  )
  
  df$Dataset <- ifelse(grepl("Dataset1", dataset_name),
                       "Dataset 1", "Dataset 2")
  
  df$CellType <- ifelse(grepl("NeuN", dataset_name),
                        "NeuN", "SOX10")
  
  plot_df <- rbind(plot_df, df)
  ## SEGMENTED ANALYSIS
  seg_results <- data.frame()
  
  model <- lm(count ~ logQ, data = df)
  
  seg_model <- segmented(
    model,
    seg.Z = ~logQ
  )
  
  
  
  seg_breaks <- abs(seg_model$psi[, "Est."])
  
  
  
  ci_segmented <- confint(seg_model)
  
  ci_df <- as.data.frame(ci_segmented)
  
  ci_lower <- ci_df[, grep(
    "low",
    colnames(ci_df),
    ignore.case = TRUE
  )]
  
  ci_upper <- ci_df[, grep(
    "up",
    colnames(ci_df),
    ignore.case = TRUE
  )]
  
  
  
  
  
  seg_labels <- sapply(
    10^(-seg_breaks),
    format_q_label
  )
  
  ci_lower_labels <- sapply(
    10^(-ci_lower),
    format_q_label
  )
  
  ci_upper_labels <- sapply(
    10^(-ci_upper),
    format_q_label
  )
  
  seg_results <- rbind(
    seg_results,
    data.frame(
      Dataset = ifelse(grepl("Dataset1", dataset_name),
                       "Dataset 1", "Dataset 2"),
      CellType = ifelse(grepl("NeuN", dataset_name),
                        "NeuN", "SOX10"),
      Break = seg_breaks,
      Lower = ci_lower,
      Upper = ci_upper,
      BreakLabel = seg_labels,
      LowerLabel = ci_lower_labels,
      UpperLabel = ci_upper_labels
    )
  )
  
  
  results_list[[dataset_name]] <- list(
    segmented_breaks = seg_breaks,
    segmented_CI = ci_segmented
  )
  
  
  
  for(i in seq_len(length(seg_breaks))) {
    
    tmp <- data.frame(
      
      Dataset = dataset_name,
      
      Breakpoint_Number = i,
      
      Breakpoint_log10Q = seg_breaks[i],
      
      Breakpoint_Qvalue = signif(
        10^(-seg_breaks[i]),
        3
      ),
      
      Breakpoint_Label = seg_labels[i],
      
      CI95_Lower_log10Q = ci_lower[i],
      
      CI95_Upper_log10Q = ci_upper[i],
      
      CI95_Lower_Qvalue = signif(
        10^(-ci_upper[i]),
        3
      ),
      
      CI95_Upper_Qvalue = signif(
        10^(-ci_lower[i]),
        3
      ),
      
      CI95_Lower_Label = ci_upper_labels[i],
      
      CI95_Upper_Label = ci_lower_labels[i]
    )
    
    summary_df <- rbind(summary_df, tmp)
    
  }
  
  
  head(summary_df)
  
  
  ymax <- max(df$count)
  
  seg_y <- seq(
    ymax * 0.95,
    ymax * 0.80,
    length.out = length(seg_breaks)
  )
  
  ci_low_y <- seq(
    ymax * 0.75,
    ymax * 0.65,
    length.out = length(ci_lower)
  )
  
  ci_up_y <- seq(
    ymax * 0.60,
    ymax * 0.50,
    length.out = length(ci_upper)
  )
  
  
  
  clean_title <- dataset_name
  
  clean_title <- gsub("_", " ", clean_title)
  
  clean_title <- gsub("Dataset1", "Dataset 1", clean_title)
  
  clean_title <- gsub("Dataset2", "Dataset 2", clean_title)
  
  # Determine colour from dataset name
  plot_colour <- case_when(
    grepl("NeuN", dataset_name, ignore.case = TRUE) ~ cell_colours["NeuN"],
    grepl("SOX10", dataset_name, ignore.case = TRUE) ~ cell_colours["SOX10"]
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  p <- ggplot(
    df,
    aes(x = logQ, y = count)
  ) +
    
    geom_point(size = 1.3,
               color = plot_colour) +
    
    geom_line(linewidth = 1.2, color = plot_colour) +
    
    
    ########################################################
  ## Segmented breakpoints
  ########################################################
  
  
  geom_vline(
    xintercept = seg_breaks,
    color = "black",
    linewidth = 1
  ) +
    
    annotate(
      "text",
      x = seg_breaks + 0.15,
      y = ymax * 1.02,
      label = seg_labels,
      color = "black",
      vjust = 0,
      hjust = 0,
      size = 4.5,
      fontface = "bold"
    ) +
    
    ########################################################
  ## 95% CI LOWER
  ########################################################
  
  geom_vline(
    xintercept = ci_lower,
    color = "darkblue",
    linewidth = 1,
    linetype = "dashed",
    alpha = 0.7
  ) +
    
    annotate(
      "text",
      x = ci_lower - 1,
      y = ci_low_y,
      label = ci_lower_labels,
      color = "black",
      hjust = 1,
      vjust = -0.3,
      size = 4
    ) +
    
    ########################################################
  ## 95% CI UPPER
  ########################################################
  
  geom_vline(
    xintercept = ci_upper,
    color = "darkblue",
    linetype = "dashed",
    linewidth = 1,
    alpha = 0.7
  ) +
    
    annotate(
      "text",
      x = ci_upper + 1,
      y = ci_up_y,
      label = ci_upper_labels,
      color = "black",
      hjust = 0,
      vjust = -0.3,
      size = 4
    ) +
    
    
    
    ########################################################
  ## X AXIS
  ########################################################
  
  scale_x_continuous(
    
    name = expression(
      paste("Q value threshold (1e"^{-n}, ")")
    ),
    
    breaks = seq(2, 100, by = 2),
    
    labels = seq(2, 100, by = 2)
  ) +
    
    ########################################################
  ## LABELS
  ########################################################
  
  labs(
    y = "Number of peaks",
    title = clean_title
  ) +
    
    ########################################################
  ## THEME
  ########################################################
  
  theme_classic(base_size = 12) +
    
    theme(
      
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
      ),
      
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
      
    )
  
  plot_list[[dataset_name]] <- p
  
  

  ## SAVE FIGURE

  library(patchwork)
  outfile <- file.path(
    outdir,
    paste0(dataset_name, "_Figure5.pdf")
  )
  
  ggsave(
    filename = outfile,
    plot = p,
    width = 10,
    height = 6
  )
  
  cat("Saved figure:\n")
  cat(outfile, "\n")
}


#plot_data <- list()


for(dataset_name in names(peak_files)) {
  
  plot_data[[dataset_name]] <- df
  
  cat(dataset_name,
      unique(df$Dataset),
      head(df$count, 1), "\n")
  
}

## ==========================
## ==========================
neun_df <- subset(plot_df, CellType == "NeuN")
sox_df  <- subset(plot_df, CellType == "SOX10")
table(sox_df$Dataset)

head(summary_df)

seg_neun <- summary_df %>%
  filter(grepl("NeuN", Dataset))

seg_sox <- summary_df %>%
  filter(grepl("SOX10", Dataset))
seg_neun$DatasetType <- ifelse(
  grepl("Dataset1", seg_neun$Dataset),
  "Dataset 1",
  "Dataset 2"
)

seg_sox$DatasetType <- ifelse(
  grepl("Dataset1", seg_sox$Dataset),
  "Dataset 1",
  "Dataset 2"
)

seg_cols <- c(
  "Dataset 1" = "black",
  "Dataset 2" = scales::alpha("black", 0.7)
)

ci_cols <- c(
  "Dataset 1" = "darkblue",
  "Dataset 2" = scales::alpha("darkblue", 0.7)
)

seg_neun$Breakpoint_Expr <- sapply(
  seg_neun$Breakpoint_Label,
  function(x){
    parts <- strsplit(x, "e")[[1]]
    paste0(parts[1], '*e^{', parts[2], '}')
  }
)

seg_sox$Breakpoint_Expr <- sapply(
  seg_sox$Breakpoint_Label,
  function(x){
    parts <- strsplit(x, "e")[[1]]
    paste0(parts[1], '*e^{', parts[2], '}')
  }
)


p_neun <- ggplot(
  neun_df,
  aes(
    x = logQ,
    y = count,
    linetype = Dataset
  )
) +
  
  geom_line(
    aes(
      colour = Dataset,
      linetype = Dataset
    ),
    linewidth = 1.3,
    guide="none"
  ) +
  
  
  scale_colour_manual(
    name = "Dataset",
    values = c(
      "Dataset 1" = "hotpink",
      "Dataset 2" = scales::alpha("hotpink", 0.65)
    ),
    guide = guide_legend(
      override.aes = list(
        linetype = c("solid", "solid"),
        linewidth = 1.5
      )
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "Dataset 1" = "solid",
      "Dataset 2" = "dashed"
    ),
    guide = "none"
  )+
  
  scale_x_continuous(
    name = expression(paste("Q value threshold (1e"^{-n}, ")")),
    breaks = seq(2,100,2),
    labels = seq(2,100,2)
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  )+
  scale_alpha_manual(
    values = c(
      "Dataset 1" = 1,
      "Dataset 2" = 0.65
    ),
    guide = "none"
  )+
  
  geom_vline(
    data = seg_neun,
    aes(
      xintercept = Breakpoint_log10Q,
      alpha = DatasetType
    ),
    colour = "red",
    linewidth = 0.7,
    show.legend = FALSE
  )+
  
  geom_vline(
    data = seg_neun,
    aes(
      xintercept = CI95_Lower_log10Q,
      alpha = DatasetType
    ),
    colour = "darkslategrey",
    linewidth = 0.5,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  
  geom_vline(
    data = seg_neun,
    aes(
      xintercept = CI95_Upper_log10Q,
      alpha = DatasetType
    ),
    colour = "darkslategrey",
    linewidth = 0.5,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  
  labs(
    title = "NeuN",
    y = "Number of peaks"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave( "p_neunTest.pdf")

p_sox <- ggplot(
  sox_df,
  aes(
    x = logQ,
    y = count,
    linetype = Dataset
  )
) +
  
  geom_line(
    aes(
      colour = Dataset,
      linetype = Dataset
    ),
    linewidth = 1.3,
    guide = "none"
  ) +
  
  scale_colour_manual(
    name = "Dataset",
    values = c(
      "Dataset 1" = "#6B8E23",
      "Dataset 2" = scales::alpha("#6B8E23", 0.65)
    ),
    guide = guide_legend(
      override.aes = list(
        linetype = c("solid", "solid"),
        linewidth = 1.5
      )
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "Dataset 1" = "solid",
      "Dataset 2" = "dashed"
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    name = expression(paste("Q value threshold (1e"^{-n}, ")")),
    breaks = seq(2, 100, 2),
    labels = seq(2, 100, 2)
  ) +
  
  scale_alpha_manual(
    values = c(
      "Dataset 1" = 1,
      "Dataset 2" = 0.7
    ),
    guide = "none"
  ) +
  
  geom_vline(
    data = seg_sox,
    aes(
      xintercept = Breakpoint_log10Q,
      alpha = DatasetType
    ),
    colour = "red",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  geom_vline(
    data = seg_sox,
    aes(
      xintercept = CI95_Lower_log10Q,
      alpha = DatasetType
    ),
    colour = "darkslategrey",
    linewidth = 0.5,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  
  geom_vline(
    data = seg_sox,
    aes(
      xintercept = CI95_Upper_log10Q,
      alpha = DatasetType
    ),
    #colour = "darkgoldenrod4",
    colour = "darkslategrey",
    linewidth = 0.5,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  
  labs(
    title = "SOX10",
    y = "Number of peaks"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

combined_plot2 <- p_neun | p_sox

panelA <- combined_plot2


ggsave(
  filename = "combined2.pdf",
  plot = combined_plot2,
  width = 20,
  height = 6
)




## 2 by 3


## SAVE SUMMARY CSV


csv_outfile <- file.path(
  outdir,
  "Figure5_threshold_summary_with95CI.csv"
)

write.csv(
  summary_df,
  file = csv_outfile,
  row.names = FALSE
)