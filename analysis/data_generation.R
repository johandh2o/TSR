# ------------------------------
# Libraries
source("used_libraries.R")

# ------------------------------
# Get the path to the directories
parent_dir = here::here()
sister_folder = file.path(parent_dir, "plots")
# Create the folder if it doesn't exist
if (!dir.exists(sister_folder)) {
  dir.create(sister_folder)
}

# ------------------------------
# Print graph
plot_1 = dagify(
  A ~ W,
  M ~ A,
  Z ~ A + M,
  Y ~ W + A + M + Z,
  R ~ M + Z
) %>% tidy_dagitty(layout = "kk") %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(color='white',size=0.5) +
  geom_dag_edges() +
  geom_dag_text(color='black') +
  theme_dag()

# Save the plot
ggsave(filename = file.path(sister_folder, "plot_1_graph.png"),
    plot = plot_1, width = 6, height = 4, dpi = 300)