plot_kde_from_file <- function(
    distance_file = "distances.csv",
    bandwidth = 0.5,
    sep = "\t"
) {
  # --- Read data ---
  df <- read.table(distance_file, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  
  if (nrow(df) == 0) {
    stop("Distance file is empty")
  }
  
  # Expected columns:
  # structure_file | base_pair | distance
  
  base_pairs <- unique(df$base_pair)
  
  max_distance <- max(df$distance)
  grid_step <- 0.1
  grid <- seq(0, max_distance + grid_step, by = grid_step)
  
  # --- Prepare plot ---
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
  
  # --- KDE per base pair ---
  for (bp in base_pairs) {
    
    distances <- df$distance[df$base_pair == bp]
    
    if (length(distances) == 0) {
      next
    }
    
    # Gaussian KDE (manual)
    diff <- outer(grid, distances, "-")
    kernel <- exp(-0.5 * (diff / bandwidth)^2)
    density <- rowSums(kernel)
    density <- density / (length(distances) * bandwidth * sqrt(2 * pi))
    
    lines(grid, density, col = colors[bp], lwd = 2)
  }
  
  legend(
    "topright",
    legend = base_pairs,
    col = colors,
    lwd = 2,
    cex = 0.8
  )
}

plot_kde_from_file()
