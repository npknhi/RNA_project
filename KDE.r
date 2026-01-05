#!/usr/bin/env Rscript

# ------------------------
# Parse command-line args
# ------------------------
args <- commandArgs(trailingOnly = TRUE)

# Defaults
distance_file <- ifelse(length(args) >= 1, args[1], "./distances.csv")
header_flag   <- ifelse(length(args) >= 2, as.logical(args[2]), TRUE)
sep_char      <- ifelse(length(args) >= 3, args[3], "\t")

bw_value <- 0.5
plot_dir <- "data/plots"

# ------------------------
# Create output directory
# ------------------------
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# ------------------------
# Read input data
# ------------------------
data.df <- read.table(
  file = distance_file,
  header = header_flag,
  sep = sep_char,
  stringsAsFactors = FALSE
)

if (nrow(data.df) == 0) {
  stop("Input file is empty")
}

# Expected columns:
# structure_file | base_pair | distance

# ------------------------
# 1) KDE + histogram (all base pairs pooled)
# ------------------------
data <- data.df$distance
kde <- density(data, bw = bw_value)

png(
  filename = file.path(plot_dir, "kde_all_base_pairs.png"),
  width = 900,
  height = 600
)

plot(
  kde,
  main = "Kernel Density Estimate (all base pairs)",
  xlab = "Distance (Å)",
  ylab = "Estimated density",
  col = "blue",
  lwd = 2
)

polygon(kde, col = rgb(0, 0, 1, 0.25), border = NA)

hist(
  data,
  probability = TRUE,
  col = rgb(1, 0, 0, 0.3),
  add = TRUE,
  breaks = 30
)

dev.off()

# ------------------------
# 2) KDE + histogram per base pair
# ------------------------
base_pairs <- unique(data.df$base_pair)

for (bp in base_pairs) {
  
  data_bp <- data.df$distance[data.df$base_pair == bp]
  if (length(data_bp) == 0) next
  
  kde <- density(data_bp, bw = bw_value)
  
  png(
    filename = file.path(plot_dir, paste0("kde_", bp, ".png")),
    width = 900,
    height = 600
  )
  
  plot(
    kde,
    main = paste("KDE for base pair", bp),
    xlab = "Distance (Å)",
    ylab = "Estimated density",
    col = "blue",
    lwd = 2
  )
  
  polygon(kde, col = rgb(0, 0, 1, 0.25), border = NA)
  
  hist(
    data_bp,
    probability = TRUE,
    col = rgb(1, 0, 0, 0.3),
    add = TRUE,
    breaks = 30
  )
  
  dev.off()
}

# ------------------------
# 3) Combined KDE plot
# ------------------------
max_distance <- max(data.df$distance)
grid_step <- 0.1
grid <- seq(0, max_distance + grid_step, by = grid_step)

png(
  filename = file.path(plot_dir, "kde_all_base_pairs_combined.png"),
  width = 1000,
  height = 700
)

plot(
  NA,
  xlim = c(0, max_distance),
  ylim = c(0, 1),
  xlab = "Distance (Å)",
  ylab = "Density",
  main = paste("KDE of inter-residue distances\n", distance_file)
)

colors <- rainbow(length(base_pairs))
names(colors) <- base_pairs

for (bp in base_pairs) {
  
  distances <- data.df$distance[data.df$base_pair == bp]
  if (length(distances) == 0) next
  
  diff <- outer(grid, distances, "-")
  kernel <- exp(-0.5 * (diff / bw_value)^2)
  density <- rowSums(kernel)
  density <- density / (length(distances) * bw_value * sqrt(2 * pi))
  
  lines(grid, density, col = colors[bp], lwd = 2)
}

legend(
  "topright",
  legend = base_pairs,
  col = colors,
  lwd = 2,
  cex = 0.8
)

dev.off()

cat("Plots successfully written to", plot_dir, "\n")