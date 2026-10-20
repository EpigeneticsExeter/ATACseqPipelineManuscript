# Fig 5 C: curated peak set widths
## Figure 5 C, peak widths

## lists of GRanges peak sets, 1 list per cell type, Using peak lists generated in script for 5A
peak_list_neun
peak_list_sox10



### peak width plots updated

library(dplyr)
library(ggplot2)


## FUNCTION: EXTRACT WIDTHS


extract_widths <- function(gr_list, celltype) {
  
  bind_rows(
    
    lapply(names(gr_list), function(nm) {
      
      gr <- gr_list[[nm]]
      
      widths <- width(gr)
      
      ## cap huge peaks
      widths[widths > 10000] <- 10000
      
      data.frame(
        
        PeakSet = nm,
        
        Width = widths,
        
        CellType = celltype,
        
        Dataset = ifelse(
          grepl("^Dat1", nm, ignore.case = TRUE),
          "Dat1",
          "Dat2"
        )
      )
    })
  )
}


## BUILD DATAFRAMES


widths_neun <- extract_widths(
  peak_list_neun,
  "NeuN"
)

widths_sox10 <- extract_widths(
  peak_list_sox10,
  "SOX10"
)




## CLEAN LABELS


clean_labels <- function(x) {
  
  x <- gsub("_", " ", x)
  x <- gsub("Dat1", "Dataset 1", x)
  x <- gsub("Dat2", "Dataset 2", x)
  x <- gsub("q", "", x)
  x <- gsub("q", "", x)
  x <- gsub("1e", "Q1e", x)
  
}

widths_neun$PeakLabel  <- clean_labels(widths_neun$PeakSet)
widths_sox10$PeakLabel <- clean_labels(widths_sox10$PeakSet)
widths_sox10$PeakLabel <- gsub(".narrowPeak.noSC.filt", "", widths_sox10$PeakLabel)
widths_sox10$PeakSet <- gsub(".narrowPeak.noSC.filt", "", widths_sox10$PeakSet)
breakpoint_df <- summary_df %>%
  mutate(
    PeakSet = case_when(
      Dataset == "NeuN_Dataset1"  ~ paste0("Dat1_NeuN_", Breakpoint_Label),
      Dataset == "NeuN_Dataset2"  ~ paste0("Dat2_NeuN_q", Breakpoint_Label),
      Dataset == "SOX10_Dataset1" ~ paste0("Dat1_SOX10_", Breakpoint_Label),
      Dataset == "SOX10_Dataset2" ~ paste0("Dat2_SOX10_q", Breakpoint_Label)
    )
  )

breakpoint_df[, c("Dataset", "Breakpoint_Label", "PeakSet")]

table(unique(widths_neun$PeakSet))
table(unique(widths_sox10$PeakSet))
#PeakSet Width CellType Dataset             PeakLabel
# Dat1_NeuN_1e-29   752     NeuN    Dat1 Dataset 1 NeuN Q1e-29
# Add breakpoint information to df
widths_neun <- widths_neun %>%
  mutate(
    is_breakpoint = PeakSet %in% breakpoint_df$PeakSet
  )

widths_sox10 <- widths_sox10 %>%
  mutate(
    is_breakpoint = PeakSet %in% breakpoint_df$PeakSet
  )

table(widths_neun$is_breakpoint)
table(widths_sox10$is_breakpoint)


## PLOT FUNCTION


plot_widths <- function(
    df,
    fill_colour,
    title_text
) {
  
  ## order datasets
  df <- df %>%
    arrange(
      Dataset,
      PeakLabel
    )
  
  df$PeakLabel <- factor(
    df$PeakLabel,
    levels = unique(df$PeakLabel)
  )
  
  
  ggplot(
    df,
    aes(
      x = Width,
      y = PeakLabel,
      fill = Dataset
    )
  ) +
    
    
    ## BOXPLOT
    
    geom_boxplot(
      width = 0.65,
      outlier.shape = NA,
      linewidth = 0.8,
      color = "black"
    ) +
    
    
    ## BREAKPOINT OUTLINE
    
    geom_boxplot(
      data = df %>% filter(is_breakpoint),
      width = 0.65,
      fill = NA,
      color = "black",
      linewidth = 3,
      outlier.shape = NA
    ) +
    
    
    ## MEDIAN POINT
    
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      color = "black"
    ) +
    
    
    ## COLORS
    
    scale_fill_manual(
      values = c(
        "Dat1" = fill_colour,
        "Dat2" = scales::alpha(fill_colour, 0.45)
      )
    ) +
    
    
    ## NO LEGEND
    
    guides(fill = "none") +
    
    
    ## AXES
    
    coord_cartesian(
      xlim = c(
        0,
        quantile(df$Width, 0.99)
      )
    ) +
    
    
    ## LABELS
    
    labs(
      x = "Peak width (bp)",
      y = NULL,
      title = title_text
    ) +
    
    
    ## THEME
    
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 10),
      legend.position = "none"
    )
}


## GENERATE PLOTS



widths_neun$Dataset <- factor(
  widths_neun$Dataset,
  levels = c("Dat2", "Dat1")
)
widths_sox10$Dataset <- factor(
  widths_sox10$Dataset,
  levels = c("Dat2", "Dat1")
)
p_width_neun <- plot_widths(
  widths_neun,
  fill_colour = "hotpink",
  title_text = "NeuN"
)

p_width_sox10 <- plot_widths(
  widths_sox10,
  fill_colour = "#6B8E23",
  title_text = "SOX10"
)




## COMBINE


library(patchwork)

panel_C <-
  
  p_width_neun +
  p_width_sox10  +
  
  plot_layout(nrow = 1)

pdf("widthtestplot.pdf", width=20, height=10)
panel_C
dev.off()

saveRDS(panel_C, "Fig.5.panel_C.rds")

## combine ABC



combined_plot.abc <-
  
  wrap_elements(full = panel_A) /
  wrap_elements(full = panel_B) /
  wrap_elements(full=panel_C) +
  
  plot_layout(
    heights = c(1.3, 1, 0.8)
  ) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    plot.tag = element_text(
      size = 22,
      face = "bold"
    )
  )



## SAVE

outdir <-  "~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/Figures/Figure5_Sept"
combined_outfile <- file.path(
  outdir,
  "Figure5_combined_ABC.pdf"
)
ggsave(
  filename = combined_outfile,
  plot = combined_plot.abc,
  width = 18,
  height = 12
)
