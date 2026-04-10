# =========================================
# US County Poverty Data Analysis (NMVM–MoLE)
# =========================================

rm(list = ls())

# ================= Load packages =================
library(ggplot2)
library(usmap)
library(dplyr)
library(patchwork)
library(tidyr)
library(GGally)
library(sf)
library(tigris)
library(stringi)

# ================= Create folders =================
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")

# ================= Read and merge data =================
edu   <- read.csv("Education2023.csv")
pop   <- read.csv("PopulationEstimates.csv")
pov   <- read.csv("Poverty2023.csv")
unemp <- read.csv("Unemployment2023.csv")

edu   <- edu   %>% rename(FIPS = FIPS.Code)
pop   <- pop   %>% rename(FIPS = FIPStxt)
pov   <- pov   %>% rename(FIPS = FIPS_Code)
unemp <- unemp %>% rename(FIPS = FIPS_Code)

edu$FIPS   <- as.character(edu$FIPS)
pop$FIPS   <- as.character(pop$FIPS)
pov$FIPS   <- as.character(pov$FIPS)
unemp$FIPS <- as.character(unemp$FIPS)

edu <- edu %>%
  filter(Attribute == "Percent of adults with a bachelor's degree or higher, 2019-23") %>%
  select(FIPS, education = Value)

pop <- pop %>%
  filter(Attribute == "POP_ESTIMATE_2023") %>%
  select(FIPS, population = Value)

pov <- pov %>%
  filter(Attribute == "PCTPOVALL_2023") %>%
  select(FIPS, poverty = Value)

unemp <- unemp %>%
  filter(Attribute == "Unemployment_rate_2023") %>%
  select(FIPS, unemployment = Value)

df <- edu %>%
  inner_join(pop,  by = "FIPS") %>%
  inner_join(pov,  by = "FIPS") %>%
  inner_join(unemp, by = "FIPS") %>%
  mutate(
    population   = as.numeric(population),
    poverty      = as.numeric(poverty),
    unemployment = as.numeric(unemployment),
    education    = as.numeric(education),
    log_pop      = log(population)
  )

summary(df)
nrow(df)

write.csv(df, "merged_county_data.csv", row.names = FALSE)

# ================= Load custom functions =================
WD.PATH <- paste0(getwd(), "/Functions")
source(paste0(WD.PATH, "/Additional.r"))
source(paste0(WD.PATH, "/GHST.r"))
source(paste0(WD.PATH, "/NIG.r"))
source(paste0(WD.PATH, "/NMVBS.r"))
source(paste0(WD.PATH, "/NMVL.r"))
source(paste0(WD.PATH, "/normal.r"))
source(paste0(WD.PATH, "/SMSN_fiting.r"))
source(paste0(WD.PATH, "/SL.r"))

# ================= Reload merged data =================
df <- read.csv("merged_county_data.csv")

# ================= Define variables =================
y <- df$poverty
x <- cbind(1, df$log_pop, df$unemployment, df$education)
r <- cbind(1, df$log_pop)

# ================= Exploratory analysis =================

# ---------- Scatter plots (NO clustering) ----------
p1 <- ggplot(df, aes(x = log_pop, y = poverty)) +
  geom_point(alpha = 0.5, size = 1.2, color = "#3E5C76") +
  theme_minimal(base_size = 9) +
  labs(x = "Log(Population)", y = "Poverty Rate (%)") +
  theme(panel.grid.minor = element_blank())

p2 <- ggplot(df, aes(x = unemployment, y = poverty)) +
  geom_point(alpha = 0.5, size = 1.2, color = "#3E5C76") +
  theme_minimal(base_size = 9) +
  labs(x = "Unemployment (%)", y = "Poverty Rate (%)") +
  theme(panel.grid.minor = element_blank())

p3 <- ggplot(df, aes(x = education, y = poverty)) +
  geom_point(alpha = 0.5, size = 1.2, color = "#3E5C76") +
  theme_minimal(base_size = 9) +
  labs(x = "Education (%)", y = "Poverty Rate (%)") +
  theme(panel.grid.minor = element_blank())

p_exploratory <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave("figures/exploratory_scatterplots.pdf",
       plot = p_exploratory,
       width = 6.5,
       height = 2.5)

ggsave("figures/exploratory_scatterplots.png",
       plot = p_exploratory,
       width = 6.5,
       height = 2.5,
       dpi = 300)

# ---------- Histogram ----------
p_hist <- ggplot(df, aes(x = poverty)) +
  geom_histogram(
    bins = 30,
    fill = "#3E5C76",
    color = "white",
    alpha = 0.85
  ) +
  theme_minimal(base_size = 9) +
  labs(
    x = "Poverty Rate (%)",
    y = "Frequency"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave("figures/ex3_histogram.pdf",
       plot = p_hist,
       width = 3.5,
       height = 2.8)

ggsave("figures/ex3_histogram.png",
       plot = p_hist,
       width = 3.5,
       height = 2.8,
       dpi = 300)

# ================= Model selection (G = 1,...,4) =================
G_list <- 1:4
CAIC_values <- matrix(NA, nrow = length(G_list), ncol = 7)
colnames(CAIC_values) <- c("Gaussian", "SN", "GHST", "NIG", "SL", "NMVBS", "NMVL")

for (i in seq_along(G_list)) {
  g <- G_list[i]
  cat("Fitting G =", g, "\n")
  
  fit_norm  <- mix.reg.norm.EM(y, x, r, g, verbose = FALSE)
  fit_sn    <- mix.reg.SMSN.EM(y, x, r, g, family = "Skew.n", verbose = FALSE)
  fit_ghst  <- mix.GHST.MoE.EM(y, x, r, g, verbose = FALSE)
  fit_nig   <- mix.NIG.MoE.EM(y, x, r, g, verbose = FALSE)
  fit_sl    <- mix.SL.MoE.EM(y, x, r, g, verbose = FALSE)
  fit_nmvbs <- mix.NMVBS.MoE.EM(y, x, r, g, verbose = FALSE)
  fit_nmvl  <- mix.NMVL.MoE.EM(y, x, r, g, verbose = FALSE)
  
  CAIC_values[i, ] <- c(
    fit_norm$CAIC,
    fit_sn$CAIC,
    fit_ghst$CAIC,
    fit_nig$CAIC,
    fit_sl$CAIC,
    fit_nmvbs$CAIC,
    fit_nmvl$CAIC
  )
}

CAIC_table <- data.frame(G = G_list, CAIC_values)
write.csv(CAIC_table, "results/ex3_CAIC_table.csv", row.names = FALSE)
print(CAIC_table)

# ================= CAIC plot =================
CAIC_long <- pivot_longer(CAIC_table, cols = -G, names_to = "Model", values_to = "CAIC")

p_caic <- ggplot(CAIC_long, aes(x = G, y = CAIC, color = Model)) +
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

ggsave("figures/CAIC_plot_example3.pdf", p_caic, width = 8, height = 5)
ggsave("figures/CAIC_plot_example3.png", p_caic, dpi = 300)
print(p_caic)

# ================= Fit final model =================
best_G <- 3

fit_norm  <- mix.reg.norm.EM(y, x, r, best_G, verbose = FALSE)
fit_sn    <- mix.reg.SMSN.EM(y, x, r, best_G, family = "Skew.n", verbose = FALSE)
fit_ghst  <- mix.GHST.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_nig   <- mix.NIG.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_sl    <- mix.SL.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_nmvbs <- mix.NMVBS.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_nmvl  <- mix.NMVL.MoE.EM(y, x, r, best_G, verbose = FALSE)

# ================= Model comparison =================
model_results <- data.frame(
  Model = c("Gaussian", "SN", "GHST", "NIG", "SL", "NMVBS", "NMVL"),
  LogLik = c(
    fit_norm$loglik, fit_sn$loglik, fit_ghst$loglik,
    fit_nig$loglik, fit_sl$loglik, fit_nmvbs$loglik, fit_nmvl$loglik
  ),
  CAIC = c(
    fit_norm$CAIC, fit_sn$CAIC, fit_ghst$CAIC,
    fit_nig$CAIC, fit_sl$CAIC, fit_nmvbs$CAIC, fit_nmvl$CAIC
  )
)

write.csv(model_results, "results/ex3_model_comparison.csv", row.names = FALSE)
print(model_results)

# ================= Choose best model =================
best_model <- fit_nmvbs
df$cluster <- factor(best_model$clusters)
print(best_model)

# ================= Cluster plot =================
cluster_colors <- c("1" = "#4C6A67", "2" = "#C06C2B", "3" = "#3A4F63")

p1_cluster <- ggplot(df, aes(x = log_pop, y = poverty, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = cluster_colors) +
  theme_minimal(base_size = 9) +
  labs(x = "Log(Population)", y = "Poverty Rate", color = "Cluster") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p2_cluster <- ggplot(df, aes(x = unemployment, y = poverty, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = cluster_colors) +
  theme_minimal(base_size = 9) +
  labs(x = "Unemployment", y = "Poverty Rate", color = "Cluster") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p3_cluster <- ggplot(df, aes(x = education, y = poverty, color = cluster)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = cluster_colors) +
  theme_minimal(base_size = 9) +
  labs(x = "Education (%)", y = "Poverty Rate", color = "Cluster") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_combined_cluster <- (p1_cluster + p2_cluster + p3_cluster) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "right")

ggsave("figures/scatterplots_clusters.png",
       plot = p_combined_cluster,
       width = 6.5,
       height = 2.5,
       dpi = 300)

ggsave("figures/scatterplots_clusters.pdf",
       plot = p_combined_cluster,
       width = 6.5,
       height = 2.5)

# ================= Pair plot =================
p_pair_cluster <- ggpairs(
  df,
  columns = c("log_pop", "unemployment", "education", "poverty"),
  mapping = aes(color = cluster)
)

ggsave("figures/pair_clusters.png", p_pair_cluster, width = 10, height = 5, dpi = 300)
ggsave("figures/pair_clusters.pdf", p_pair_cluster, width = 10, height = 5)

# ================= Gating function =================
tau <- best_model$alpha
grid_x <- seq(min(df$log_pop), max(df$log_pop), length.out = 200)

eta1 <- tau[1, 1] + tau[2, 1] * grid_x
eta2 <- tau[1, 2] + tau[2, 2] * grid_x

den <- 1 + exp(eta1) + exp(eta2)

pi1 <- exp(eta1) / den
pi2 <- exp(eta2) / den
pi3 <- 1 / den

df_gate <- data.frame(
  log_pop = rep(grid_x, 3),
  prob = c(pi1, pi2, pi3),
  cluster = factor(rep(1:3, each = length(grid_x)))
)

p_gate <- ggplot(df_gate, aes(x = log_pop, y = prob, color = cluster)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = cluster_colors) +
  theme_minimal(base_size = 9) +
  labs(
    x = "Log(Population)",
    y = "Probability",
    color = "Cluster"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave("figures/ex3_gating.pdf",
       plot = p_gate,
       width = 3.5,
       height = 2.8)

ggsave("figures/ex3_gating.png",
       plot = p_gate,
       width = 3.5,
       height = 2.8,
       dpi = 300)

# ================= Uncertainty analysis =================
posterior_prob <- best_model$tau
uncertainty <- 1 - apply(posterior_prob, 1, max)

df$uncertainty <- uncertainty
df$dominant_cluster <- apply(posterior_prob, 1, which.max)

uncertainty_scale <- scale_color_gradient(
  low = "#F7F7F7",
  high = "#3E5C76"
)

p1_uncertainty_heatmap <- ggplot(df, aes(x = log_pop, y = education, color = uncertainty)) +
  geom_point(alpha = 0.7, size = 3) +
  uncertainty_scale +
  theme_minimal() +
  labs(
    x = "Log(Population)",
    y = "Education (%)",
    color = "Uncertainty"
  ) +
  theme(panel.grid.minor = element_blank())

p2_uncertainty_heatmap <- ggplot(df, aes(x = log_pop, y = unemployment, color = uncertainty)) +
  geom_point(alpha = 0.7, size = 3) +
  uncertainty_scale +
  theme_minimal() +
  labs(
    x = "Log(Population)",
    y = "Unemployment (%)",
    color = "Uncertainty"
  ) +
  theme(panel.grid.minor = element_blank())

p3_uncertainty_heatmap <- ggplot(df, aes(x = education, y = unemployment, color = uncertainty)) +
  geom_point(alpha = 0.7, size = 3) +
  uncertainty_scale +
  theme_minimal() +
  labs(
    x = "Education (%)",
    y = "Unemployment (%)",
    color = "Uncertainty"
  ) +
  theme(panel.grid.minor = element_blank())

p_uncertainty_heatmap <- p1_uncertainty_heatmap +
  p2_uncertainty_heatmap +
  p3_uncertainty_heatmap +
  plot_layout(ncol = 3)

ggsave("figures/ex3_uncertainty_heatmap.png",
       p_uncertainty_heatmap,
       width = 10,
       height = 5,
       dpi = 300)

ggsave("figures/ex3_uncertainty_heatmap.pdf",
       p_uncertainty_heatmap,
       width = 10,
       height = 5)

# ================= Prediction and residuals =================
beta <- best_model$beta
y_pred <- rowSums(best_model$tau * (x %*% beta))

residuals <- y - y_pred
df_res <- data.frame(index = 1:length(residuals), res = residuals)

p_residual <- ggplot(df_res, aes(x = index, y = res)) +
  geom_point(size = 1.2, alpha = 0.7, color = "#3E5C76") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#C06C2B", linewidth = 0.6) +
  theme_minimal(base_size = 9) +
  labs(
    x = "Observation Index",
    y = "Residuals"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave("figures/ex3__residual.pdf",
       plot = p_residual,
       width = 3.5,
       height = 2.8)

ggsave("figures/ex3__residual.png",
       plot = p_residual,
       width = 3.5,
       height = 2.8,
       dpi = 300)

# ================= US county clustering map =================
options(tigris_use_cache = TRUE)

counties_map <- counties(cb = TRUE, resolution = "20m", year = 2022, class = "sf")

df$county <- read.csv("merged_county_data.csv")$county
df_map <- df

df_map$county <- iconv(df_map$county, from = "", to = "UTF-8", sub = "byte")
df_map$county <- stri_trans_general(df_map$county, "Latin-ASCII")
df_map$county <- gsub(" County", "", df_map$county)
df_map$county <- tolower(df_map$county)

df_map <- df_map %>%
  filter(county != "hawaii") %>%
  group_by(county) %>%
  slice(1) %>%
  ungroup()

counties_map$NAME <- tolower(counties_map$NAME)

map_data <- counties_map %>%
  left_join(df_map, by = c("NAME" = "county")) %>%
  filter(!is.na(cluster))

coords <- st_coordinates(st_centroid(map_data))
xlim <- quantile(coords[, 1], probs = c(0.01, 0.99), na.rm = TRUE)
ylim <- quantile(coords[, 2], probs = c(0.01, 0.99), na.rm = TRUE)

p_map <- ggplot(map_data) +
  geom_sf(aes(fill = as.factor(cluster)), color = "white", linewidth = 0.1) +
  scale_fill_manual(
    values = c("1" = "#5E7D6A", "2" = "#BFA27A", "3" = "#7C90A6"),
    name = "Cluster",
    labels = c("1", "2", "3")
  ) +
  labs(title = "U.S. County Clustering Map") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    axis.text = element_text(size = 8, color = "gray40"),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE
  )

ggsave("figures/ex3_us_county_clusters.png",
       plot = p_map,
       width = 12,
       height = 6,
       dpi = 300)

cairo_pdf("figures/ex3_us_county_clusters.pdf", width = 12, height = 6)
print(p_map)
dev.off()

# ================= End =================
cat("Example 3 completed.\n")
