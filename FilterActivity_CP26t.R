
# new fig test 29sept
library(dplyr)
library(stringr)
library(ggplot2)

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


