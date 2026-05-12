


# decoupleR + collecTRI
# In this script we are doing:
#            1-) DEG analysis among 21 unique idents
#            2-) Use LogFolChange values of the DEG list to feed into decoupleR and collecTRI
#            3-) results are not embeded into seurat as there are less genes
#            4-) we saved the matrix and dataframe of the results 

# https://saezlab.github.io/decoupleR/articles/tf_sc.html

library(Seurat)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(seecolor)
library(decoupleR)
library(viper)
library(colorRamp2)
library(ComplexHeatmap)
library(ggrepel)

#LOAD NC DATA all genes:
nc.data <-readRDS(file="ncObject.RDS")

Idents(nc.data) <- "combined"
my_levels <- c("HH5_undecided","1som_undecided",'4som_undecided',
               "HH5_nne",'1som_nne',"4som_nne","7som_nne","HH11_nne",
               "4som_npb","7som_nc","HH11_nc",
               "HH5_ne",'1som_ne',"4som_ne","7som_ne","HH11_ne",
               "HH5_Vnt",'1som_Vnt',"4som_Vnt","7som_Vnt","HH11_Vnt")
nc.data@active.ident <- factor(x = nc.data@active.ident, levels = my_levels)
nc.data@meta.data[["combined"]] <-Idents(nc.data)

# CollecTRI is a comprehensive resource containing a curated collection of TFs 
# and their transcriptional targets compiled from 12 different resources

net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)


 #------------ using DEG log fold changes ------------------------------------
#nc.data.markers <- FindAllMarkers(nc.data, only.pos = T,min.pct = 0.25, logfc.threshold = 0.2)
nc.data.markers <- read.table("USETHIS_deg_StagesCells.txt",header=TRUE,row.names=1) 
nc.data.markers <- subset(nc.data.markers, p_val_adj <0.05)

# Step 1: Pivot to wide matrix (genes × clusters)
# Make wide matrix: genes as rows, clusters as columns
deg_mat <- nc.data.markers %>%
  select(gene, cluster, avg_log2FC) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames('gene') %>%
  as.matrix()


deg_mat[is.na(deg_mat)] <- 0         # Replace NA with 0 (neutral fold change)
deg_mat[is.infinite(deg_mat)] <- 0   # Replace Inf/-Inf with 0 as w

# Step 2: Run decoupleR on this matrix
# Run decoupleR (ULM method) on log2FC matrix
acts_deg <- decoupleR::run_ulm(
  mat = deg_mat,
  net = net,
  .source = 'source',
  .target = 'target',
  .mor = 'mor'
)
#Now acts_deg will contain:  TF × cluster activities based on log2 fold changes.

tf_activity_df <- acts_deg %>%
  dplyr::rename(TF = source, Cluster = condition, Activity = score) %>%
  dplyr::arrange(desc(Activity))  # optional sorting

#write.csv(tf_activity_df, file = "tf_activity_deg.csv", row.names = FALSE)


#:::::::::::::::::::: test RNA exp vs TF activity scatter plot::::::::::::::::::::::::::::::::::::::::

# Average expression per cluster (ident). This gives log-normalized expression
tf_expr_avg <- AverageExpression(nc.data, features = unique(tf_activity_df$TF), group.by = "ident")$RNA

# Convert to long format to have ::: TF, Cluster, Expression
tf_expr_df <- tf_expr_avg %>%
  as.data.frame() %>%
  tibble::rownames_to_column("TF") %>%
  pivot_longer(cols = -TF, names_to = "Cluster", values_to = "Expression")

# for some reason here ident names start with number gets "g" letter in beginning lets remove that:
tf_expr_df$Cluster <- sub("^g(\\d)", "\\1", tf_expr_df$Cluster)  # remove 'g' from beginning
tf_expr_df$Cluster <- gsub("-", "_", tf_expr_df$Cluster) # format fix 


# Join TF activity with expression to have ::: TF, Cluster, Activity, Expression
tf_combined_df <- tf_activity_df %>%
  rename(Activity = Activity) %>%
  inner_join(tf_expr_df, by = c("TF", "Cluster"))

# Here we observe dataset distribution to label cutoff for "high" and "low" for activity and rna expression (if needed)
activity_high_cutoff <- quantile(tf_combined_df$Activity, 0.60, na.rm = TRUE)
activity_low_cutoff  <- quantile(tf_combined_df$Activity, 0.10, na.rm = TRUE)

expression_high_cutoff <- quantile(tf_combined_df$Expression, 0.85, na.rm = TRUE)
expression_low_cutoff  <- quantile(tf_combined_df$Expression, 0.40, na.rm = TRUE)

# Assign status (only 4 defined groups)
tf_combined_df <- tf_combined_df %>%
  mutate(Status = case_when(
    Expression > 0.5 & Activity > 0.5 ~ "High Expr + High Act",
    Expression < 0.2 & Activity > 0.5 ~ "Low Expr + High Act",
    Expression > 0.5 & Activity < 0.2 ~ "High Expr + Low Act",
    Expression < 0.2 & Activity < 0.2 ~ "Low Expr + Low Act"
  )) %>%
  filter(!is.na(Status))  # ✅ removes "Intermediate" rows

# order 
tf_combined_df$Cluster <- factor(tf_combined_df$Cluster, 
                                 levels = c("HH5_undecided","1som_undecided",'4som_undecided',
                                            "HH5_nne",'1som_nne',"4som_nne","7som_nne","HH11_nne",
                                            "4som_npb","7som_nc","HH11_nc",
                                            "HH5_ne",'1som_ne',"4som_ne","7som_ne","HH11_ne",
                                            "HH5_Vnt",'1som_Vnt',"4som_Vnt","7som_Vnt","HH11_Vnt"))

#And then plot with:
ggplot(tf_combined_df, aes(x = Expression, y = Activity, color = Status)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~Cluster, scales = "free") +
  theme_minimal()

# lets scale the x axis (expression)
ggplot(tf_combined_df, aes(x = Expression + 0.1, y = Activity, color = Status)) +  # add 0.1 to avoid log(0)
  geom_point(alpha = 0.7) +
  facet_wrap(~Cluster, scales = "free") +
  scale_x_log10() +  # 🔥 log-scale for expression
  geom_vline(xintercept = log10(0.5 + 0.1), linetype = "dashed", color = "gray60") +  # match your cutoffs
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c(
    "High Expr + High Act" = "#a83279",
    "Low Expr + High Act" = "#3b4cc0",
    "High Expr + Low Act" = "#208f8c",
    "Low Expr + Low Act" = "#dddddd"
  )) +
  theme_minimal(base_size = 12) +
  labs(
    title = "TF Activity vs Expression per Cluster (log-scaled)",
    x = "TF Expression (log10 scale)",
    y = "TF Activity Score",
    color = "TF Status"
  )

    # Make SMALL GROUPS for plots: 
group_colors <- c(
  "High Expr + High Act" = "#a83279",  # magenta-berry
  "Low Expr + High Act"  = "#3b4cc0",  # indigo blue
  "High Expr + Low Act"  = "#208f8c",  # teal
  "Low Expr + Low Act"   = "#cccccc"   # light grey
)



xxx_groups <- c("HH5_nne", "1som_nne", "4som_nne", "7som_nne", "HH11_nne")
xxx_groups <- c("4som_npb", "7som_nc", "HH11_nc")
xxx_groups <- c("HH5_ne",'1som_ne',"4som_ne","7som_ne","HH11_ne")
xxx_groups <- c("HH5_Vnt",'1som_Vnt',"4som_Vnt","7som_Vnt","HH11_Vnt")
xxx_groups <-c("HH5_undecided","1som_undecided",'4som_undecided')


plot_df <- tf_combined_df %>%
  filter(Cluster %in% xxx_groups)
plot_df$Cluster <- factor(plot_df$Cluster, levels = xxx_groups)


ggplot(plot_df, aes(x = Expression + 0.1, y = Activity, color = Status)) +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~Cluster, ncol = 5, scales = "free") +
  scale_x_log10() +
  geom_vline(xintercept = 0.5 + 0.1, linetype = "dashed", color = "gray60") +  
  geom_vline(xintercept = 0.2 + 0.1, linetype = "dashed", color = "gray60") +  
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "gray60") +
  scale_color_manual(values = group_colors) +  # defined earlier
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    legend.title = element_text(face = "bold")
  )+
  labs(
    title = "TF Activity vs Expression — NNE Clusters",
    x = "TF Expression (log10 scale)",
    y = "TF Activity Score",
    color = "TF Status"
  )




# save the groups from plot above::::::::::::::::::::::::::::::::::::::::::::::::::::::
# This gives you 30 total TF-cluster points (10 from each of the 3 groups), with the top scores.
# for the smaller groups!!! 
# we calculate each sep then merge:

top_high_expr_high_act <- plot_df %>%
  filter(Status == "High Expr + High Act") %>%
  group_by(Cluster) %>%
  arrange(desc(Activity), desc(Expression)) %>%
  slice_head(n = 5)

top_low_expr_high_act <- plot_df %>%
  filter(Status == "Low Expr + High Act") %>%
  group_by(Cluster) %>%
  arrange(desc(Activity), Expression) %>%  # activity high, expression low
  slice_head(n = 5)

top_high_expr_low_act <- plot_df %>%
  filter(Status == "High Expr + Low Act") %>%
  group_by(Cluster) %>%
  arrange(Activity, desc(Expression)) %>%  # activity low, expression high
  slice_head(n = 5)


top_labels <- bind_rows(
  top_high_expr_high_act,
  top_low_expr_high_act,
  top_high_expr_low_act
)

#add labels to scatter plot 




ggplot(plot_df, aes(x = Expression + 0.1, y = Activity, color = Status)) +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~Cluster, ncol = 5, scales = "free") +
  scale_x_log10() +
  geom_vline(xintercept = 0.5 + 0.1, linetype = "dashed", color = "gray60") +  
  geom_vline(xintercept = 0.2 + 0.1, linetype = "dashed", color = "gray60") +  
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "gray60") +
  geom_text_repel(
    data = top_labels,
    aes(label = TF),
    size = 3,       # Small font
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = "gray60"
  ) +
  scale_color_manual(values = group_colors) +  # defined earlier
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    legend.title = element_text(face = "bold")
  )+
  labs(
    title = "TF Activity vs Expression — NC Clusters",
    x = "TF Expression (log10 scale)",
    y = "TF Activity Score",
    color = "TF Status"
  )

#::::::::::::::::::::::::::: heatmaps for those small groups ::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Color scale settings
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- colorRampPalette(colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))
#........................


tfs_high <- unique(top_high_expr_high_act$TF)
tfs_low <- unique(top_high_expr_low_act$TF)
#


top_high_mat <- tf_combined_df %>%
  filter(TF %in% tfs_high, Cluster %in% xxx_groups) %>%
  pivot_wider(names_from = TF, values_from = Activity, id_cols = Cluster) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()
top_high_mat[is.na(top_high_mat)] <-0
top_high_mat <- top_high_mat[xxx_groups,]


top_low_mat <- tf_combined_df %>%
  filter(TF %in% tfs_low, Cluster %in% xxx_groups) %>%
  pivot_wider(names_from = TF, values_from = Activity, id_cols = Cluster) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()
top_low_mat[is.na(top_low_mat)] <-0
top_low_mat <- top_low_mat[xxx_groups,]


# Plot High Expr + High Act
pheatmap::pheatmap(top_high_mat,
                   color = colors.use,
                   breaks = my_breaks,
                   border_color = "white",
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   main = "High Expr + High Act TFs")

# Plot Low Expr + High Act
pheatmap::pheatmap(top_low_mat,
                   color = colors.use,
                   breaks = my_breaks,
                   border_color = "white",
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   main = "High Expr + Low Act TFs")


#dot plot rna exp and tf activity:::::::::::::::::::::::::::::::::::::::::::
x_groups <- c("HH5_nne", "1som_nne", "4som_nne", "7som_nne", "HH11_nne")

# Step 1: Create full TF × cluster grid
full_grid <- expand.grid(
  TF = tfs,
  Cluster = x_groups,
  stringsAsFactors = FALSE
)

# Step 2: Join actual filtered values (may be missing)
dot_df <- full_grid %>%
  left_join(
    tf_combined_df %>% filter(TF %in% tfs, Cluster %in% x_groups),
    by = c("TF", "Cluster")
  )

# Step 3: Fill in missing with original values from the full tf_expr_df 
dot_df <- dot_df %>%
  left_join(
    tf_expr_df  %>%
      select(TF, Cluster, Expression_full = Expression),
    by = c("TF", "Cluster")
  ) %>%
  mutate(
    Expression = coalesce(Expression, Expression_full)
  ) %>%
  select(-Expression_full)

# Step 4: Fill any remaining NA activity with 0 (assume neutral)
dot_df <- dot_df %>%
  mutate(Activity = replace_na(Activity, 0))

# Step 5: Bin expression into 3 levels for size
dot_df <- dot_df %>%
  mutate(ExpressionBin = cut(
    Expression,
    breaks = quantile(Expression, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE),
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE
  ))

# Step 6: Order TFs (Y-axis)
dot_df$TF <- factor(dot_df$TF, levels = rev(final_gene_names))
dot_df$Cluster <- factor(dot_df$Cluster, levels = x_groups)


# Step 7: Plot
ggplot(dot_df, aes(x = Cluster, y = TF)) +
  geom_point(
    aes(size = ExpressionBin, fill = Activity),
    shape = 21,
    color = "black",
    stroke = 0.3
  ) +
  scale_fill_gradient2(
    low = "#00ff00",   # greenish low
    mid = "white",
    high = "#bf8bff",  # purple high
    midpoint = 0,
    name = "TF Activity"
  ) +
  scale_size_manual(
    values = c("Low" = 5, "Medium" = 7, "High" = 10),
    name = "RNA Expression"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "TF Activity (color) and RNA Expression (size)",
    x = "Cluster",
    y = "Transcription Factor"
  )


#Plot heatmap for all clusters from list of TFs::::::::::::::::::::::::::::::::::::::::::::::::::
tfs <- common_all
tfs <- c("NANOG","POU5F3","KLF4")
tfs <- c("COL9A1","VCP","PMEL","FABP7")


orderGrp = c("HH5_undecided","1som_undecided",'4som_undecided',
           "HH5_nne",'1som_nne',"4som_nne","7som_nne","HH11_nne",
           "4som_npb","7som_nc","HH11_nc",
           "HH5_ne",'1som_ne',"4som_ne","7som_ne","HH11_ne",
           "HH5_Vnt",'1som_Vnt',"4som_Vnt","7som_Vnt","HH11_Vnt")

activity_mat <- tf_combined_df %>%
  filter(TF %in% tfs) %>%
  select(Cluster, TF, Activity) %>%
  tidyr::pivot_wider(names_from = TF, values_from = Activity, values_fill = 0) %>%
  tibble::column_to_rownames("Cluster") %>%
  as.matrix()

activity_mat <-activity_mat[orderGrp,]

# Color scale settings
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- colorRampPalette(colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))
#........................
pheatmap::pheatmap(
  mat = activity_mat,
  color = colors.use,       # your custom RdBu or viridis palette
  breaks = my_breaks,       # same breakpoints for consistency
  border_color = "white",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  cellwidth = 15,
  cellheight = 15,
  treeheight_row = 10,
  treeheight_col = 10,
  main = "TF Activity"
)


#



















# Creative plot: Activity + p-value comparison per group :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

tf_combined_df <- read.csv("tf_combined_deg.csv")
# Make SMALL GROUPS for plots: 

xxx_groups <- c("HH5_nne", "1som_nne", "4som_nne", "7som_nne", "HH11_nne")

plot_df <- tf_combined_df %>%
  filter(Cluster %in% xxx_groups)

# PLOT option: Horizontal bar plot (diverging bars, color shade is pval)__________
#filter to top & bottom TFs per cluster for clarity
# These are absolute 10 vals meaning + or - 
plot_df_top_bottom <- plot_df %>%
  group_by(Cluster) %>%
  slice_max(order_by = abs(Activity), n = 10) %>%
  ungroup()

ggplot(plot_df_top_bottom, aes(x = Activity, y = reorder(TF, Activity))) +
  geom_col(aes(fill = -log10(p_value)), width = 0.7) +
  facet_wrap(~Cluster, nrow = 1, scales = "free_y") +  # force 1-row layout
  scale_fill_viridis_c(option = "C", name = "-log10(p-value)") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top/Bottom TF Activity per Cluster with Significance",
    x = "TF Activity Score",
    y = "TF (ordered by activity)"
  )

# # PLOT option: Horizontal bubble plot (bubble size is pval)__________
# These are absolute 10 vals meaning + or - 

xxx_groups <- c("HH5_nne", "1som_nne", "4som_nne", "7som_nne", "HH11_nne")
xxx_groups <- c("4som_npb", "7som_nc", "HH11_nc")
xxx_groups <- c("HH5_ne",'1som_ne',"4som_ne","7som_ne","HH11_ne")


plot_df2 <- tf_combined_df %>%
  filter(Cluster %in% xxx_groups) %>%
  group_by(Cluster) %>%
  slice_max(order_by = abs(Activity), n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Direction = ifelse(Activity >= 0, "High", "Low"),
    Significance = -log10(p_value)
  )


# more edited version of the plot: 
plot_df2$Cluster <- factor(plot_df2$Cluster, levels = xxx_groups)


ggplot(plot_df2, aes(x = Activity, y = reorder(TF, Activity))) +
  geom_segment(aes(x = 0, xend = Activity, y = TF, yend = TF),color = "gray70", size = 0.4) + # Midline-to-bubble segment
  geom_point(aes(size = Significance, color = Direction), alpha = 0.9) +  # Bubble
  facet_wrap(~Cluster, nrow = 1, scales = "free_y") + # Facets in one row with spacing
  scale_color_manual(values = c("High" = "#8add42", "Low" = "#c60356")) + # Color scale for activity direction
  scale_size_continuous(range = c(2, 8), name = "-log10(p-value)") + # Size scale for significance
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") + # Midline
  theme_minimal(base_size = 12) +# Clean minimal theme with no background grid
  theme(
    panel.grid = element_blank(),              # ❌ remove grid
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),  # ✅ facet border
    strip.background = element_blank(),        # cleaner strip
    panel.spacing = unit(1, "lines")           # spacing between facets
  ) +
  labs(# Labels
    title = "TF Activity (Bubble Plot) – Nt Clusters",
    x = "TF Activity Score",
    y = "TF",
    color = "Activity Direction"
  )

# :::::::::::::::::::: below is for more balanced top and bottom 5 TF activity scores:::::

# Top 5 high-activity TFs per cluster
top_tfs <- tf_combined_df %>%
  filter(Cluster %in% xxx_groups) %>%
  group_by(Cluster) %>%
  slice_max(order_by = Activity, n = 5, with_ties = FALSE) %>%
  ungroup()

# Top 5 low-activity TFs per cluster
bottom_tfs <- tf_combined_df %>%
  filter(Cluster %in% xxx_groups) %>%
  group_by(Cluster) %>%
  slice_min(order_by = Activity, n = 5, with_ties = FALSE) %>%
  ungroup()

# Combine and prepare for plotting
plot_df_balanced <- bind_rows(top_tfs, bottom_tfs) %>%
  mutate(
    Direction = ifelse(Activity >= 0, "High", "Low"),
    Significance = -log10(p_value)
  )

# plot the same plot: :::
plot_df_balanced$Cluster <- factor(plot_df_balanced$Cluster, levels = xxx_groups)

ggplot(plot_df_balanced, aes(x = Activity, y = reorder(TF, Activity))) +
  geom_segment(aes(x = 0, xend = Activity, y = TF, yend = TF),color = "gray70", size = 0.4) + # Midline-to-bubble segment
  geom_point(aes(size = Significance, color = Direction), alpha = 0.9) +  # Bubble
  facet_wrap(~Cluster, nrow = 1, scales = "free_y") + # Facets in one row with spacing
  scale_color_manual(values = c("High" = "#8add42", "Low" = "#c60356")) + # Color scale for activity direction
  scale_size_continuous(range = c(2, 8), name = "-log10(p-value)") + # Size scale for significance
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") + # Midline
  theme_minimal(base_size = 12) +# Clean minimal theme with no background grid
  theme(
    panel.grid = element_blank(),              # ❌ remove grid
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),  # ✅ facet border
    strip.background = element_blank(),        # cleaner strip
    panel.spacing = unit(1, "lines")           # spacing between facets
  ) +
  labs(# Labels
    title = "TF Activity (Bubble Plot) – NC Clusters",
    x = "TF Activity Score",
    y = "TF",
    color = "Activity Direction"
  )



#:::::::::::::::::::::: lets do abs top 15 and balanced top and bottom n = 10 for globalllll::::::::::::::::::::::::::::::

tf_combined_df <- read.csv("tf_combined_deg.csv")
tf_combined_df$Cluster <- factor(tf_combined_df$Cluster, 
                                 levels = c("HH5_undecided","1som_undecided",'4som_undecided',
                                            "HH5_nne",'1som_nne',"4som_nne","7som_nne","HH11_nne",
                                            "4som_npb","7som_nc","HH11_nc",
                                            "HH5_ne",'1som_ne',"4som_ne","7som_ne","HH11_ne",
                                            "HH5_Vnt",'1som_Vnt',"4som_Vnt","7som_Vnt","HH11_Vnt"))
# Color scale settings
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- colorRampPalette(colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

#___________________Top 15 TFs by absolute mean activity across all clusters_________________________________________________
# Compute mean activity per TF across all clusters
top15_tfs <- tf_combined_df %>%
  group_by(TF) %>%
  summarise(mean_activity = mean(Activity, na.rm = TRUE)) %>%
  arrange(desc(abs(mean_activity))) %>%
  slice_head(n = 30) %>%
  pull(TF)

top15_mat <- tf_combined_df %>%
  filter(TF %in% top15_tfs) %>%
  select(Cluster, TF, Activity) %>%
  pivot_wider(names_from = TF, values_from = Activity, id_cols = Cluster) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

pheatmap(top15_mat,
         color = colors.use,
         breaks = my_breaks,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         cellwidth = 15,
         cellheight = 15,
         border_color = "white",
         main = "Top TFs by Absolute Mean Activity")

#___________________Balanced version — top 15 and bottom 15 TFs by mean activity across all clusters_________________________________________________
# Get top 15 highest mean activity TFs
top_up <- tf_combined_df %>%
  group_by(TF) %>%
  summarise(mean_activity = mean(Activity, na.rm = TRUE)) %>%
  arrange(desc(mean_activity)) %>%
  slice_head(n = 20) %>%
  pull(TF) %>%
  unique()


top_high_mat <- tf_combined_df %>%
  filter(TF %in% top_up) %>%
  pivot_wider(names_from = TF, values_from = Activity, id_cols = Cluster) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

# Plot High Expr + High Act
pheatmap::pheatmap(top_high_mat,
                   color = colors.use,
                   breaks = my_breaks,
                   border_color = "white",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   main = "Top TF activity (balanced)")

# Get bottom 15 lowest mean activity TFs
top_down <- tf_combined_df %>%
  group_by(TF) %>%
  summarise(mean_activity = mean(Activity, na.rm = TRUE)) %>%
  arrange(mean_activity) %>%
  slice_head(n = 20)%>%
  pull(TF) %>%
  unique()

top_low_mat <- tf_combined_df %>%
  filter(TF %in% top_down) %>%
  pivot_wider(names_from = TF, values_from = Activity, id_cols = Cluster) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

# Plot Low Expr + High Act
pheatmap::pheatmap(top_low_mat,
                   color = colors.use,
                   breaks = my_breaks,
                   border_color = "white",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   main = "Low TF activity (balanced)")



#_______________________________________________________-

# Median activity and expression across all TF-cluster pairs
median_activity <- median(tf_combined_df$Activity, na.rm = TRUE)
median_expression <- median(tf_combined_df$Expression, na.rm = TRUE)

cat("Median Activity:", median_activity, "\n")
cat("Median Expression:", median_expression, "\n")

quantile(tf_combined_df$Expression, probs = 0.75, na.rm = TRUE)



ggplot(tf_combined_df, aes(x = Activity)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Distribution of TF Activity Scores", x = "Activity")

ggplot(tf_combined_df, aes(x = Expression)) +
  geom_histogram(bins = 50, fill = "seagreen", color = "white") +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Distribution of TF Expression", x = "Expression")
