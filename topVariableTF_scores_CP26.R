

# Most variable TFs without filter
library(stringr)
tf_combined_df <- read.csv("tf_combined_deg.csv") 

# extract stage and cell type
tf_combined_df <- tf_combined_df %>%
  mutate(
    Stage = sub("_(.*)", "", Cluster),
    CellType = sub(".*?_", "", Cluster)
  )

# Function: top varying TFs per group
get_top_varying_tfs <- function(df, group_var, top_n = 7) {
  df %>%
    group_by(!!sym(group_var), TF) %>%
    summarise(AvgAct = mean(Activity), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(group_var), values_from = AvgAct) %>%
    column_to_rownames("TF") %>%
    as.data.frame() -> mat
  
  # Replace NAs with row-wise minimum
  for (i in 1:nrow(mat)) {
    row_min <- min(mat[i, ], na.rm = TRUE)
    mat[i, is.na(mat[i, ])] <- row_min
  }
  
  mat$StdDev <- apply(mat, 1, sd)
  mat %>%
    arrange(desc(StdDev)) %>%
    slice_head(n = top_n * (ncol(mat) - 1)) %>%  # subtract 1 for StdDev col
    select(-StdDev)
}

# 1) AMONG EACH STAGE ------------------------------------------------------------------------------------------------------------
mat_stage <- as.matrix(get_top_varying_tfs(tf_combined_df, "Stage"))
pheatmap(mat_stage, scale = "row", main = "Top Varying TFs per Stage")


# 2) AMONG EACH CELL TYPE ---------------------------------------------------------------------------------------------------------
mat_cell <- as.matrix(get_top_varying_tfs(tf_combined_df, "CellType"))
pheatmap(mat_cell, scale = "row", main = "Top Varying TFs per Cell Type")




#:: Lollipop plots in one row, boxed panels ::::::::::::::::::::::::::::
df_cell <- tf_combined_df %>%
  mutate(
    StageOnly = sub("_(.*)", "", Cluster),
    CellType = str_extract(Cluster, "[^_]+$"))%>%
  filter(Expression > 0)

# Desired stage order
stage_levels <- c("HH5", "1som", "4som", "7som", "HH11")

# Average per StageOnly x TF (ignore cell type)
df_stage <- df_cell %>%
  mutate(StageOnly = factor(StageOnly, levels = stage_levels)) %>%
  group_by(StageOnly, TF) %>%
  summarise(
    avg_activity = mean(Activity, na.rm = TRUE),
    avg_expr     = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  )

# Select top 5 and bottom 5 by activity per stage
top5 <- df_stage %>%
  group_by(StageOnly) %>%
  slice_max(order_by = avg_activity, n = 30, with_ties = FALSE) %>%
  mutate(Group = "Top 5 (positive)")


for (st in stage_levels) {
  tf_list <- unique(top5$TF[top5$StageOnly == st])
  out_file <- paste0("Fig2f_Top50_", st, ".txt")
  write.table(tf_list, out_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE)}

for (st in stage_levels) {
  df_stage <- subset(top5, StageOnly == st)
  out_file <- paste0("Fig2F_Top30all_", st, ".csv")
  write.csv(df_stage, out_file, row.names = FALSE)}


bottom5 <- df_stage %>%
  group_by(StageOnly) %>%
  slice_min(order_by = avg_activity, n = 30, with_ties = FALSE) %>%
  mutate(Group = "Bottom 5 (negative)")


for (st in stage_levels) {
  tf_list <- unique(bottom5$TF[top5$StageOnly == st])
  out_file <- paste0("Fig2f_Bottom30_", st, ".txt")
  write.table(tf_list, out_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE)}
#csv file for supp
for (st in stage_levels) {
  df_stage <- subset(bottom5, StageOnly == st)
  out_file <- paste0("Fig2F_Bottom30all_", st, ".csv")
  write.csv(df_stage, out_file, row.names = FALSE)}

plot_df <- bind_rows(top5, bottom5) %>%
  ungroup() %>%
  mutate(
    Group = factor(Group, levels = c("Top 5 (positive)", "Bottom 5 (negative)")),
    # order TFs within each facet by their activity
    TF = fct_reorder(TF, avg_activity)
  )

ggplot(plot_df, aes(x = avg_activity, y = TF)) +
  geom_segment(aes(x = 0, xend = avg_activity, y = TF, yend = TF),
               color = "gray65", linewidth = 1) +
  geom_point(aes(size = avg_expr, fill = Group),
             shape = 21, color = "black", stroke = 0.2, alpha = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray50") +
  scale_fill_manual(values = c("Top 5 (positive)" = "#7A3B0C",
                               "Bottom 5 (negative)" = "#2D86C0"),
                    name = NULL) +
  scale_size_continuous(range = c(3, 9), name = "Avg RNA (stage)") +
  facet_wrap(~ StageOnly, nrow = 1, scales = "free_y") +
  labs(x = "Average TF Activity (per stage)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),                    # remove grid
    strip.background = element_blank(),              # remove facet label box
    strip.text = element_text(face = "bold"),        # keep stage names bold
    panel.spacing = unit(20, "pt"),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )
