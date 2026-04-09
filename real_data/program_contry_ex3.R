# =========================================
# US County Poverty Data Analysis (NMVM–MoLE)
# =========================================

rm(list = ls())

#================= Load packages =================
library(ggplot2)
library(usmap)
library(dplyr)
#================= Merge data =================

#======== Read data =================
edu   <- read.csv("Education2023.csv")
pop   <- read.csv("PopulationEstimates.csv")
pov   <- read.csv("Poverty2023.csv")
unemp <- read.csv("Unemployment2023.csv")

#========= Rename FIPS =================
edu   <- edu   %>% rename(FIPS = FIPS.Code)
pop   <- pop   %>% rename(FIPS = FIPStxt)
pov   <- pov   %>% rename(FIPS = FIPS_Code)
unemp <- unemp %>% rename(FIPS = FIPS_Code)

edu$FIPS   <- as.character(edu$FIPS)
pop$FIPS   <- as.character(pop$FIPS)
pov$FIPS   <- as.character(pov$FIPS)
unemp$FIPS <- as.character(unemp$FIPS)

#====== FILTER  =================
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

#======= Merge =================
df <- edu %>%
  inner_join(pop,  by = "FIPS") %>%
  inner_join(pov,  by = "FIPS") %>%
  inner_join(unemp, by = "FIPS")

#======= Clean =================
df <- df %>%
  mutate(
    population = as.numeric(population),
    poverty = as.numeric(poverty),
    unemployment = as.numeric(unemployment),
    education = as.numeric(education),
    log_pop = log(population)
  )

#====== Check =================
summary(df)
nrow(df)

#===== Save =================
write.csv(df, "merged_county_data.csv", row.names = FALSE)

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
df <- read.csv("merged_county_data.csv")

#================= Define variables =================

y <- df$poverty

x <- cbind(1, df$log_pop, df$unemployment, df$education)
r <- cbind(1, df$log_pop)


#================= Exploratory plot =================

library(ggplot2)
library(patchwork)

p1 <- ggplot(df, aes(x = log_pop, y = poverty)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(x = "Log(Population)", y = "Poverty Rate (%)", title = "Poverty vs Log(Population)")

p2 <- ggplot(df, aes(x = unemployment, y = poverty)) +
  geom_point(alpha = 0.6, size = 2, color = "darkred") +
  theme_minimal() +
  labs(x = "Unemployment (%)", y = "Poverty Rate (%)", title = "Poverty vs Unemployment")

p3 <- ggplot(df, aes(x = education, y = poverty)) +
  geom_point(alpha = 0.6, size = 2, color = "darkgreen") +
  theme_minimal() +
  labs(x = "Education (%)", y = "Poverty Rate (%)", title = "Poverty vs Education")

p_combined <- p1 + p2 + p3 + plot_layout(ncol = 3)
ggsave("figures/exploratory_scatterplots.png", p_combined, width = 15, height = 5, dpi = 300)
ggsave("figures/exploratory_scatterplots.pdf", p_combined, width = 10, height = 4)

# Histogram
p_hist <- ggplot(df, aes(x = poverty)) +
  geom_histogram(bins = 30, fill = "gray70", color = "black") +
  theme_minimal() +
  labs(
    x = "Poverty Rate (%)",
    y = "Frequency",
    title = "Distribution of poverty rate"
  )

ggsave("figures/ex3_histogram.png", p_hist, width = 10, height = 4, dpi = 300)
ggsave("figures/ex3_histogram.pdf", p_hist, width = 10, height = 4)
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
write.csv(CAIC_table, "results/ex3_CAIC_table.csv", row.names = FALSE)
print(CAIC_table)
# =====================================
# CAIC Plot for Example 3 (US Counties)
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
ggsave("figures/CAIC_plot_example3.pdf", p, width = 8, height = 5)
ggsave("figures/CAIC_plot_example3.png", p, dpi = 300)

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

write.csv(model_results, "results/ex3_model_comparison.csv", row.names = FALSE)
print(model_results)

#================= Choose best model =================
best_model <- fit_nmvbs
df$cluster <- as.factor(best_model$clusters)
print(best_model)
#================= Cluster plot =================
p1_cluster <- ggplot(df, aes(x = log_pop, y = poverty, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(x = "Log(Population)", y = "Poverty Rate", color = "Cluster")

p2_cluster <- ggplot(df, aes(x = unemployment, y = poverty, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(x = "Unemployment", y = "Poverty Rate", color = "Cluster")

p3_cluster <- ggplot(df, aes(x = education, y = poverty, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(x = "Education", y = "Poverty Rate", color = "Cluster")

p_combined_cluster <- p1_cluster + p2_cluster + p3_cluster + plot_layout(ncol = 3)
ggsave("figures/scatterplots_clusters.png", p_combined_cluster, width = 15, height = 5, dpi = 300)
ggsave("figures/scatterplots_clusters.pdf", p_combined_cluster, width = 10, height = 4)

#########################
library(GGally)
p_pair_cluster=ggpairs(df, columns = c("log_pop", "unemployment", "education", "poverty"),
        mapping = aes(color = cluster))
ggsave("figures/pair_clusters.png", p_pair_cluster, width = 10, height = 5, dpi = 300)
ggsave("figures/pair_clusters.pdf", p_pair_cluster, width = 10, height = 5)

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

ggsave("figures/ex3_gating.png", p_gate, dpi = 300)
ggsave("figures/ex3_gating.pdf", p_gate, width = 8, height = 5)

# ================= Uncertainty Analysis =================
posterior_prob <- best_model$tau  

uncertainty <- 1 - apply(posterior_prob, 1, max)
df$uncertainty <- uncertainty

df$dominant_cluster <- apply(posterior_prob, 1, which.max)

library(ggplot2)

p1_uncertainty_heatmap <- ggplot(df, aes(x = log_pop, y = education, color = uncertainty)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(
    x = "Log(Population)",
    y = "Education (%)",
    color = "Uncertainty",
    title = ""
  )

p2_uncertainty_heatmap <- ggplot(df, aes(x = log_pop, y = unemployment, color = uncertainty)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(
    x = "Log(Population)",
    y = "Unemployment",
    color = "Uncertainty",
    title = ""
  )

p3_uncertainty_heatmap <- ggplot(df, aes(x = education, y = unemployment, color = uncertainty)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(
    x = "Log(Population)",
    y = "Unemployment",
    color = "Uncertainty",
    title = ""
  )

p_uncertainty_heatmap <- p1_uncertainty_heatmap + p2_uncertainty_heatmap + p3_uncertainty_heatmap + plot_layout(ncol = 3)

ggsave("figures/ex3_uncertainty_heatmap.png", p_uncertainty_heatmap, width = 10, height = 5, dpi = 300)
ggsave("figures/ex3_uncertainty_heatmap.pdf", p_uncertainty_heatmap, width = 10, height = 5)

#================= Residuals =================
y_pred <- rowSums(best_model$tau * (x %*% beta))
residuals <- y - y_pred
rmse <- sqrt(mean((y - y_pred)^2))
cat("RMSE =", rmse, "\n")

df_res <- data.frame(index = 1:nrow(df), res = residuals)

p_res <- ggplot(df_res, aes(x = index, y = res)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal()

ggsave("figures/ex3_residuals.png", p_res, dpi = 300)
ggsave("figures/ex3_residuals.pdf", p_res, width = 8, height = 5)

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
df$county <- read.csv("merged_county_data.csv")$county

df_map <- df

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
  geom_sf(aes(fill = as.factor(cluster)), color = "black", size = 0.1) +
  scale_fill_manual(
    values = c("1" = "#5E7D6A", "2" = "#BFA27A", "3" = "#7C90A6"),
    name = "Cluster",
    labels = c("1", "2", "3")
  ) +
  labs(title = "US County Clustering Map") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)


# ================= Save plots =================
ggsave("figures/ex3_us_county_clusters.png", plot = p_map, width = 12, height = 6, dpi = 300)
cairo_pdf("figures/ex3_us_county_clusters.pdf", width = 12, height = 6)
grid.draw(p_map)
dev.off()
#================= End =================
cat("Example 3 completed.\n")
