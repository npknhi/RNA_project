# install.packages("languageserver")

# R Script: Compute and Plot Kernel Density Estimate (KDE)

set.seed(123)
# Example data: 1000 random points from normal distribution
data <- rnorm(1000, mean = 0, sd = 1)

# Compute KDE
kde <- density(data)                            # density() computes kernel density estimate
kde <- density(data, kernel = "epanechnikov")   # alternative kernel : Guassian Kernel by default
kde <- density(data, bw = 0.5)                  # smaller bw = more sensitive, larger = smoother | Control bandwith smoothing

# Print KDE values (optional)
head(kde)

# Plot KDE
plot(kde,
     main = "Kernel Density Estimate",
     xlab = "Value",
     ylab = "Estimated density",
     col = "blue",
     lwd = 2)

# Fill under the curve
polygon(kde, col = rgb(0,0,1,0.2), border = NA)

# 5Overlay histogram
hist(data, probability = TRUE, col = rgb(1,0,0,0.3), add = TRUE)
