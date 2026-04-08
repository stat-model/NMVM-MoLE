# =========================================
# US County Poverty Data Analysis (NMVM–MoLE)
# =========================================

rm(list = ls())

#================= Load packages =================
library(ggplot2)
library(usmap)
library(dplyr)

#================= Load custom functions =================
WD.PATH <- paste0(getwd(), "/Functions")
source(paste0(WD.PATH, "/Additional.r"))
source(paste0(WD.PATH, "/GHST.r"))
source(paste0(WD.PATH, "/NIG.r"))
source(paste0(WD.PATH, "/NMVBS.r"))
source(paste0(WD.PATH, "/NMVL.r"))
source(paste0(WD.PATH, "/normal.r"))
source(paste0(WD.PATH, "/SMSN_fiting.r"))
source(paste0(WD.PATH, "/SL.r"))

if(!dir.exists("figures")) dir.create("figures")
if(!dir.exists("results")) dir.create("results")

#================= Load data =================
county_pop=read_csv("county_population_poverty.csv")

df <- data.frame(
  population = county_pop$pop_2022,
  log_pop = county_pop$log_population,
  poverty = county_pop$pct_pov_2021
)

df <- na.omit(df)

#================= Define variables =================

x <- cbind(1, df$log_pop)
r <- x
y <- df$poverty

n <- length(y)

#================= Exploratory plot =================
library(ggplot2)
library(patchwork)

# Scatter plot
p_scatter <- ggplot(df, aes(x = log_pop, y = poverty)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(
    x = "Log(Population)",
    y = "Poverty Rate (%)",
    title = "Scatter plot"
  )

# Histogram
p_hist <- ggplot(df, aes(x = poverty)) +
  geom_histogram(bins = 30, fill = "gray70", color = "black") +
  theme_minimal() +
  labs(
    x = "Poverty Rate (%)",
    y = "Frequency",
    title = "Distribution of poverty rate"
  )

p_combined <- p_scatter + p_hist + plot_layout(ncol = 2)

ggsave("figures/ex2_combined.png", p_combined, width = 10, height = 4, dpi = 300)
ggsave("figures/ex2_combined.pdf", p_combined, width = 10, height = 4)
#================= Model selection (G = 1,...,4) =================
G_list <- 1:4
CAIC_values <- matrix(NA, nrow = length(G_list), ncol = 7)
colnames(CAIC_values) <- c("Gaussian", "SN", "GHST", "NIG", "SL", "NMVBS","NMVL")

for(i in seq_along(G_list)){
  g <- G_list[i]
  cat("Fitting G =", g, "\n")
  
  fit_norm <- mix.reg.norm.EM(y, x, r, g, verbose = FALSE)
  fit_sn   <- mix.reg.SMSN.EM(y, x, r, g, family = "Skew.n", verbose = FALSE)
  fit_ghst <- mix.GHST.MoE.EM(y, x, r, g, verbose = FALSE)
  fit_nig  <- mix.NIG.MoE.EM(y, x, r, g, verbose = FALSE)
  fit_sl   <- mix.SL.MoE.EM(y, x, r, g, verbose = FALSE)
  fit_nmvbs<- mix.NMVBS.MoE.EM(y, x, r, g, verbose = FALSE)
  fit_nmvl<- mix.NMVL.MoE.EM(y, x, r, g, verbose = FALSE)
  
  
  CAIC_values[i, ] <- c(fit_norm$CAIC,
                        fit_sn$CAIC,
                        fit_ghst$CAIC,
                        fit_nig$CAIC,
                        fit_sl$CAIC,
                        fit_nmvbs$CAIC,
                        fit_nmvl$CAIC)
}
CAIC_table <- data.frame(G = G_list, CAIC_values)
write.csv(CAIC_table, "results/ex2_CAIC_table.csv", row.names = FALSE)
print(CAIC_table)
# =====================================
# CAIC Plot for Example 2 (US Counties)
# =====================================

library(ggplot2)
library(tidyr)

# CAIC data
CAIC_df <- CAIC_table
# Convert to long format for ggplot
CAIC_long <- pivot_longer(CAIC_df, cols = -G, names_to = "Model", values_to = "CAIC")

# Plot
p <- ggplot(CAIC_long, aes(x = G, y = CAIC, color = Model)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "CAIC Values by Number of Clusters (G)",
    x = "Number of Clusters (G)",
    y = "CAIC",
    color = "Model"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# Save
ggsave("figures/CAIC_plot_example2.pdf", p, width = 8, height = 5)
ggsave("figures/CAIC_plot_example2.png", p, dpi = 300)

# Show plot
print(p)
#================= Fit final model =================
best_G <- 3

fit_norm <- mix.reg.norm.EM(y, x, r, best_G, verbose = FALSE)
fit_sn   <- mix.reg.SMSN.EM(y, x, r, best_G, family = "Skew.n", verbose = FALSE)
fit_ghst <- mix.GHST.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_nig  <- mix.NIG.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_sl   <- mix.SL.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_nmvbs<- mix.NMVBS.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_nmvl<- mix.NMVL.MoE.EM(y, x, r, best_G, verbose = FALSE)

#================= Model comparison =================
model_results <- data.frame(
  Model = c("Gaussian", "SN", "GHST", "NIG", "SL", "NMVBS","NMVL"),
  LogLik = c(fit_norm$loglik, fit_sn$loglik, fit_ghst$loglik,
             fit_nig$loglik, fit_sl$loglik, fit_nmvbs$loglik,fit_nmvl$loglik),
  CAIC = c(fit_norm$CAIC, fit_sn$CAIC, fit_ghst$CAIC,
           fit_nig$CAIC, fit_sl$CAIC, fit_nmvbs$CAIC, fit_nmvl$CAIC)
)

write.csv(model_results, "results/ex2_model_comparison.csv", row.names = FALSE)
print(model_results)

#================= Choose best model =================
best_model <- fit_nmvbs
df$cluster <- as.factor(best_model$clusters)
print(best_model)
#================= Cluster plot =================
p_cluster <- ggplot(df, aes(x = log_pop, y = poverty, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(x = "Log(Population)", y = "Poverty Rate", color = "Cluster")

ggsave("figures/ex2_clusters.png", p_cluster, dpi = 300)
ggsave("figures/ex2_clusters.pdf", p_cluster, width = 8, height = 5)

#================= Regression lines =================
beta <- best_model$beta

p_reg <- ggplot(df, aes(x = log_pop, y = poverty, color = cluster)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = beta[1,1], slope = beta[2,1], color = "red", size = 1.2) +
  geom_abline(intercept = beta[1,2], slope = beta[2,2], color = "green", size = 1.2) +
  geom_abline(intercept = beta[1,3], slope = beta[2,3], color = "blue", size = 1.2) +
  
  theme_minimal()

ggsave("figures/ex2_regression.png", p_reg, dpi = 300)
ggsave("figures/ex2_regression.pdf", p_reg, width = 8, height = 5)

#================= Gating function =================
tau <- best_model$alpha

grid_x <- seq(min(df$log_pop), max(df$log_pop), length.out = 200)

pi2 <- 1 / (1 + exp(-(tau[1,1] + tau[2,1]*grid_x)))
pi3 <- 1 / (1 + exp(-(tau[1,2] + tau[2,2]*grid_x)))

pi1 <- 1 - pi2 - pi3

df_gate <- data.frame(
  log_pop = rep(grid_x, 3),
  prob = c(pi1, pi2, pi3),
  cluster = factor(rep(1:3, each = length(grid_x)))
)

library(ggplot2)
p_gate <- ggplot(df_gate, aes(x = log_pop, y = prob, color = cluster)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  labs(y = "Probability", color = "Cluster") +
  scale_color_manual(values = c("red", "green", "blue"))

ggsave("figures/ex2_gating.png", p_gate, dpi = 300)
ggsave("figures/ex2_gating.pdf", p_gate, width = 8, height = 5)

# ================= Uncertainty Analysis =================
posterior_prob <- best_model$tau  

uncertainty <- 1 - apply(posterior_prob, 1, max)
df$uncertainty <- uncertainty

df$dominant_cluster <- apply(posterior_prob, 1, which.max)

library(ggplot2)
p_uncertainty <- ggplot(df, aes(x = log_pop, y = uncertainty, color = factor(dominant_cluster))) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(y = "Classification Uncertainty", color = "Dominant Cluster") +
  scale_color_manual(values = c("red", "green", "blue"))

ggsave("figures/ex2_uncertainty.png", p_uncertainty, dpi = 300)
ggsave("figures/ex2_uncertainty.pdf", p_uncertainty, width = 8, height = 5)

#================= Residuals =================
y_pred <- rowSums(best_model$tau * (x %*% beta))
residuals <- y - y_pred
rmse <- sqrt(mean((y - y_pred)^2))
cat("RMSE =", rmse, "\n")

df_res <- data.frame(index = 1:n, res = residuals)

p_res <- ggplot(df_res, aes(x = index, y = res)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal()

ggsave("figures/ex2_residuals.png", p_res, dpi = 300)
ggsave("figures/ex2_residuals.pdf", p_res, width = 8, height = 5)

# ================= US County Clustering Map =================
library(ggplot2)
library(sf)
library(tigris)
library(dplyr)
library(stringi)
library(grid)

options(tigris_use_cache = TRUE)

# ================= Prepare map =================
counties_map <- counties(cb = TRUE, resolution = "20m", year = 2022, class = "sf")

# ================= Prepare df_map =================
df_map <- df
df_map$county <- county_pop$county

df_map$county <- iconv(df_map$county, from = "", to = "UTF-8", sub = "byte")
df_map$county <- stri_trans_general(df_map$county, "Latin-ASCII")
df_map$county <- gsub(" County", "", df_map$county)
df_map$county <- tolower(df_map$county)

df_map <- df_map %>% filter(county != "hawaii")

df_map <- df_map %>% group_by(county) %>% slice(1) %>% ungroup()

counties_map$NAME <- tolower(counties_map$NAME)

# ================= Join data with map =================
map_data <- counties_map %>%
  left_join(df_map, by = c("NAME" = "county"))

map_data <- map_data %>% filter(!is.na(cluster))

coords <- st_coordinates(st_centroid(map_data))
xlim <- quantile(coords[,1], probs = c(0.01, 0.99), na.rm = TRUE)
ylim <- quantile(coords[,2], probs = c(0.01, 0.99), na.rm = TRUE)

# ================= Plot map =================
p_map <- ggplot(map_data) +
  geom_sf(aes(fill = as.factor(cluster)), color = "white", size = 0.1) +
  scale_fill_manual(
    values = c("1" = "#8A0000", "2" = "#008A00", "3" = "#004D8A"),
    name = "Cluster",
    labels = c("1", "2", "3")
  ) +
  labs(title = "US County Clustering Map") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE
  )

# ================= Save plots =================
ggsave("figures/us_county_clusters.png", plot = p_map, width = 12, height = 6, dpi = 300)
cairo_pdf("figures/us_county_clusters.pdf", width = 12, height = 6)
grid.draw(p_map)
dev.off()
#================= End =================
cat("Example 2 completed.\n")
