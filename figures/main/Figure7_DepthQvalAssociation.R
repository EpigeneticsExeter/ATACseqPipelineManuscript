## Figure 7
# Association between total aligned reads and peak set specific Q value
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/Figures/QC")
Meth_P_pass1 <- read_csv("~/Library/CloudStorage/OneDrive-UniversityofExeter/ATAC/methodsPaper/Figures/QC/Meth_P_pass1.csv")
Meth_P_pass1 <- as.data.frame(Meth_P_pass1)
head(Meth_P_pass1)
depth <- Meth_P_pass1 %>%
  group_by(CellType) %>%
  summarise(
    total_aligned_reads = sum(aligned_total, na.rm = TRUE),
    total_aligned_fragments = sum(aligned_total, na.rm = TRUE) / 2
  )

depth <- as.data.frame(depth)
depth$CellType
#"IRF8"  "NEUN"  "SOX10" "TN"

depth$Q <- c(1e-19, 1e-17, 1e-21, 1e-13 )
library(lme4)
# association between fragments and Q value
# statistical and visual
head(depth)
mod1 <- lm(Q ~total_aligned_fragments, data=depth )
summary(mod1)
cor.test()


####
depth <- depth %>%
  mutate(
    neglog10Q = -log10(Q)
  )

depth$dataset <- "D2"


depth
cor.test(
  depth$total_aligned_fragments,
  depth$neglog10Q,
  method = "pearson"
)

cor.test(
  depth$total_aligned_fragments,
  depth$neglog10Q,
  method = "spearman",
  exact = FALSE
)

library(ggplot2)

pdf("depth_Qvalue2.pdf")
ggplot(
  depth,
  aes(
    x = total_aligned_fragments,
    y = neglog10Q,
    label = CellType
  )
) +
  geom_point(size = 4) +
  geom_text(
    vjust = -0.8,
    size = 4
  ) +
  scale_x_continuous(
    labels = scales::comma
  ) +
  labs(
    x = "Total aligned fragments",
    y = expression(-log[10](Q)),
    title = "Sequencing depth and peak-calling threshold"
  ) +
  theme_classic()

ggplot(depth, aes(
  x = total_aligned_fragments,
  y = neglog10Q,
  label = CellType
)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 4) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1,
    show.legend = FALSE
  ) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Total aligned fragments",
    y = expression(-log[10](Q)),
    title = "Sequencing depth and peak-calling threshold"
  ) +
  theme_classic()
dev.off()


## add BDR
passAllQCstatsBDR <- read_csv("passAllQCstatsBDR.csv")
passAllQCstatsBDR <- as.data.frame(passAllQCstatsBDR)

depthBDR <- passAllQCstatsBDR %>%
  group_by(fraction) %>%
  summarise(
    total_aligned_reads = sum(TotalNumReadsStat, na.rm = TRUE),
    total_aligned_fragments = sum(TotalNumReadsStat, na.rm = TRUE) / 2
  )

## exact p value (not <0.001)
get_lm_pvalue <- function(df) {
  
  lm_fit <- lm(
    neglog10Q ~ total_aligned_fragments,
    data = df
  )
  
  p_value <- summary(lm_fit)$coefficients[
    "total_aligned_fragments",
    "Pr(>|t|)"
  ]
  
  p_label <- paste0(
    "p = ",
    format(p_value, scientific = TRUE, digits = 3)
  )
  
  return(p_label)
}

############################################################
# Calculate p-values for each panel
############################################################

pval_combined <- as.factor(get_lm_pvalue(depthD1D2))



# plot combined
pdf("depth_Qvalue_D1D2_smallLabel.pdf")

head(depthD1D2)
# colour poi
p1 <- ggplot(depthD1D2, aes(
  x = total_aligned_fragments,
  y = neglog10Q,
  label = data_label
)) +
  geom_point(aes(color = dataset), size = 4) +
  geom_text(vjust = -0.8, size = 2.5) +
  geom_smooth(
    aes(
      x = total_aligned_fragments,
      y = neglog10Q
    ),
    method = "lm",
    se = FALSE,
    linewidth = 1,
    colour = "black",
    inherit.aes = FALSE,
    show.legend = FALSE
  )  +
  scale_color_manual(
    values = c(
      "D1" = "goldenrod",
      "D2" = "brown"
    )
  ) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Total aligned fragments",
    y = expression(-log[10](Q)),
    color = "Dataset",
    title = "Datasets 1 & 2"
  ) +
  theme_classic() +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = pval_combined,
    hjust = 1.7,
    vjust = 1.5,
    size = 3.5
  ) 
p1
dev.off()


depthBDR$CellType
depthBDR$Q <- c(1e-9,1e-31,1e-23,1e-15)

depthBDR <- depthBDR %>%
  mutate(
    neglog10Q = -log10(Q)
  )

depthBDR$dataset <- "D1"
colnames(depthBDR) <- colnames(depth)
depthD1D2 <- rbind(depthBDR, depth)
depthD1D2$data_label <- paste0(depthD1D2$CellType,"_",depthD1D2$dataset)


## datasets separately
pval_D1 <- as.factor(get_lm_pvalue(BDR))
BDR <- depthD1D2[which(depthD1D2$dataset=="D1"),]
pdf("D1_depth_qval.pdf")
p2 <- ggplot(BDR, aes(
  x = total_aligned_fragments,
  y = neglog10Q,
  label = data_label
)) +
  geom_point(aes(color = dataset), size = 4) +
  geom_text(vjust = -0.8, size = 2.5) +
  geom_smooth(
    aes(
      x = total_aligned_fragments,
      y = neglog10Q
    ),
    method = "lm",
    se = FALSE,
    linewidth = 1,
    colour = "black",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "D1" = "goldenrod",
      "D2" = "brown"
    )
  ) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Total aligned fragments",
    y = expression(-log[10](Q)),
    color = "Dataset",
    title = "Dataset 1"
  ) +
  theme_classic()  +
  theme(
    legend.position = "none" 
  )+
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = pval_D1,
    hjust = 1.7,
    vjust = 1.5,
    size = 3.5
  ) 
p2
dev.off()

SFARI <- depthD1D2[which(depthD1D2$dataset=="D2"),]
pval_D2 <- as.factor(get_lm_pvalue(SFARI))
pdf("D2_depth_qval.pdf")

p3 <- ggplot(SFARI, aes(
  x = total_aligned_fragments,
  y = neglog10Q,
  label = data_label
)) +
  geom_point(aes(color = dataset), size = 4) +
  geom_text(vjust = -0.8, size = 2.5) +
  geom_smooth(
    aes(
      x = total_aligned_fragments,
      y = neglog10Q
    ),
    method = "lm",
    se = FALSE,
    linewidth = 1,
    colour = "black",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "D1" = "goldenrod",
      "D2" = "brown"
    )
  ) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Total aligned fragments",
    y = expression(-log[10](Q)),
    color = "Dataset",
    title = "Dataset 2"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = pval_D2,
    hjust = 1.7,
    vjust = 1.5,
    size = 3.5
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

print(p3)

dev.off()

# combine p1,p2,p3 in a single figure, 1 row of 3, with a single legend

# Combine plots into a single figure (1 row x 3 columns)


library(patchwork)

# Keep the legend only in p1
p1 <- p1 +
  theme(legend.position = "right")

p2 <- p2 +
  theme(legend.position = "none")

p3 <- p3 +
  theme(legend.position = "none")

# Combine plots and collect the single dataset legend
p_combined <- (p1 | p2 | p3) +
  plot_layout(
    nrow = 1,
    guides = "collect"
  ) +
  plot_annotation(
    title = "Relationship between sequencing depth and curated peak calling thresholds",
  ) 
#&
theme(
  legend.position = "bottom"
)


# Export combined figure


ggsave(
  filename = "depth_Qvalue_D1D2_combined.pdf",
  plot = p_combined,
  width = 18,
  height = 7,
  units = "in"
)

print(p_combined)