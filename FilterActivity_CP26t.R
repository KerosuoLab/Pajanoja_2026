
# new fig test 29sept
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)

# Read data
tf_combined_df <- read.csv("tf_combined_deg.csv")

# Choose a single cell type
cell_to_cluster <- "nne"

# Filter for selected cell type and only rows with expression > 0
df_cell <- tf_combined_df %>%
  mutate(
    StageOnly = sub("_(.*)", "", Cluster),
    CellType = str_extract(Cluster, "[^_]+$")
  ) %>%
  filter(CellType == cell_to_cluster, Expression > 0)

# Identify TFs with:
# - At least one positive activity
# - At least one negative activity
# - Non-NA activity in all 5 stages
# - Expression > 0 in all 5 stages

tf_stage_switchers <- df_cell %>%
  filter(!is.na(Activity)) %>%
  group_by(TF) %>%
  summarise(
    stages_with_activity = n_distinct(StageOnly),
    any_positive = any(Activity > 0),
    any_negative = any(Activity < 0),
    expr_positive_all = all(Expression > 0)
  ) %>%
  filter(
    stages_with_activity == 5,
    any_positive == TRUE,
    any_negative == TRUE,
    expr_positive_all == TRUE
  )

interesting_tfs <- tf_stage_switchers$TF

# Subset the original filtered data for only the interesting TFs
interesting_tfs_df <- df_cell %>%
  filter(TF %in% interesting_tfs)



# Plot
stage_colors <- c('HH5' = '#e6c029', '1som' = '#fa7148', '4som' = '#4cd2ff','7som' = '#2f7cff',"HH11"="#d04dfc")
df_cell$StageOnly <- factor(df_cell$StageOnly, levels = names(stage_colors))

ggplot(interesting_tfs_df, aes(x = StageOnly, y = Activity, fill = StageOnly)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ TF, scales = "free_y") +
  scale_fill_manual(values = stage_colors) +
  theme_minimal() +
  labs(
    title = paste("Stage-Switching TFs in", cell_to_cluster),
    y = "TF Activity Score"
  )



# Select top 20 TFs based on average RNA expression across all 5 stages
top_20_tfs <- df_cell %>%
  filter(TF %in% interesting_tfs, Expression > 0) %>%
  group_by(TF) %>%
  summarise(
    AvgExpression = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(AvgExpression)) %>%
  slice_head(n = 15) %>%
  pull(TF)

# Filter only top 20 for plotting
top_20_df <- df_cell %>%
  filter(TF %in% top_20_tfs)


ggplot(top_20_df, aes(x = StageOnly, y = Activity, fill = StageOnly)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ TF, scales = "free_y") +
  scale_fill_manual(values = stage_colors) +
  theme_minimal() +
  labs(
    title = paste("Stage-Switching TFs in", cell_to_cluster),
    y = "TF Activity Score"
  )



ggplot(top_20_df, aes(x = TF, Activity, y = Activity)) +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.4) +
  geom_segment(aes(xend = TF, y = 0, yend = Activity,),
               color = "gray40", linewidth = 1) +
  geom_point(aes(size = Expression, fill = StageOnly),
             shape = 21, color = "black", alpha = 0.9, stroke = 0.3) +
  scale_size_continuous(range = c(4, 8)) +
  scale_fill_manual(values = stage_colors) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    x = NULL,
    y = "TF Activity Score",
    size = "RNA Expression",
    fill = "Stage"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 13),
    legend.position = "bottom"
  )




# ---- USER INPUT ---- choose the stage and gene list to analyze ---------------------
stage_to_cluster <- "HH5"  
gene_list <- unique(readLines("HH5genes.txt"))
#
stage_to_cluster <- "1som"  
gene_list <- unique(readLines("1ssgenes.txt"))
#
stage_to_cluster <- "4som"  
gene_list <- unique(readLines("4ssgenes.txt"))
#
stage_to_cluster <- "7som"  
gene_list <- unique(readLines("7ssgenes.txt"))
#
stage_to_cluster <- "HH11"  
gene_list <- unique(readLines("HH11genes.txt"))



# ---- Filter for chosen stage ----
df_stage <- tf_combined_df %>%
  mutate(
    StageOnly = sub("_(.*)", "", Cluster),
    CellType = str_extract(Cluster, "[^_]+$")) %>%
  filter(StageOnly == stage_to_cluster, Expression > 0)

#--- run below---------
df_geneplot <- df_stage %>%
  filter(TF %in% gene_list) %>%
  mutate(
    TF = factor(TF, levels = rev(gene_list))  # maintain order
  )

ggplot(df_geneplot, aes(x = TF, y = Activity)) +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.4) +
  geom_segment(aes(xend = TF, y = 0, yend = Activity, color = CellType),
               linewidth = 1.1, show.legend = FALSE) +  # no legend for lines
  geom_point(aes(size = Expression, fill = CellType),
             shape = 21, color = "black", alpha = 0.9, stroke = 0.3,
             show.legend = c(fill = FALSE)) +           # drop fill legend
  scale_fill_manual(values = celltype_colors) +
  scale_color_manual(values = celltype_colors) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    x = NULL,
    y = "TF Activity Score",
    size = "RNA Expression"  # only this will remain
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom"
  )

#----------------------------------------------------------------------------
table(df_stage$CellType)
# Define
target <- "nne"
neighbors <- c("ne", "nc","Vnt")

# 1. TFs active in target
active_in_target <- df_stage %>%
  filter(CellType == target, Expression > 0, Activity > 0) %>%
  pull(TF) %>% unique()

# 2. For each TF, check if it is expressed in ALL neighbors and Activity ≤ 0 in ALL
inactive_in_all_neighbors <- df_stage %>%
  filter(CellType %in% neighbors, Expression > 0) %>%
  group_by(TF) %>%
  summarise(AllInactive = all(Activity <= 0), .groups = "drop") %>%
  filter(AllInactive) %>%
  pull(TF)

# 3. Final intersect
differential_tfs <- intersect(active_in_target, inactive_in_all_neighbors)

top_x2 <- df_stage %>%
  filter(TF %in% differential_tfs, CellType == target,Expression > 0) %>%
  group_by(TF) %>%
  summarise(AvgExpression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(AvgExpression)) %>%
  slice_head(n = 5) %>%
  pull(TF)


df_geneplot <- df_stage %>%
  filter(TF %in% top_x2) %>%
  mutate(
    TF = factor(TF, levels = rev(top_x2))  # maintain order
  )


