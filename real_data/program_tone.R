# =========================================
# Tone Perception Data Analysis (NMVM–MoLE)
# =========================================

rm(list = ls())

#================= Load packages =================
library(ggplot2)

#================= Load custom functions =================
WD.PATH <- paste0(getwd(), "/Functions")
source(paste0(WD.PATH, "/Additional.r"))
source(paste0(WD.PATH, "/GHST.r"))
source(paste0(WD.PATH, "/NIG.r"))
source(paste0(WD.PATH, "/NMVBS.r"))
source(paste0(WD.PATH, "/normal.r"))
source(paste0(WD.PATH, "/SMSN_fiting.r"))
source(paste0(WD.PATH, "/SL.r"))

if(!dir.exists("figures")) dir.create("figures")
if(!dir.exists("results")) dir.create("results")

#================= Load data =================
tone_dataset <- read.csv("tone_dataset.csv")
n <- nrow(tone_dataset)

x <- cbind(1, tone_dataset$stretchratio)
r <- x
y <- tone_dataset$tuned

#================= Exploratory analysis =================
p_scatter <- ggplot(tone_dataset, aes(x = stretchratio, y = tuned)) +
  geom_point(alpha = 0.7, size = 2.5) +
  theme_minimal() +
  labs(
    x = "Actual Tone Ratio (Stretch Ratio)",
    y = "Perceived Tone Ratio"
  )

ggsave("figures/01_scatter_plot.pdf", p_scatter, width = 6, height = 5)
ggsave("figures/01_scatter_plot.png", p_scatter, dpi = 300)

#================= Model selection (varying G) =================
G_list <- 1:5
CAIC_values <- matrix(NA, nrow = length(G_list), ncol = 6)
colnames(CAIC_values) <- c("Gaussian", "SN", "GHST", "NIG", "SL", "NMVBS")

for(i in seq_along(G_list)){
  g <- G_list[i]
  cat("Fitting models with G =", g, "\n")
  
  fit_norm <- mix.reg.norm.EM(y, x, r, g, verbose = FALSE)
  fit_sn   <- mix.reg.SMSN.EM(y, x, r, g, family = "Skew.n", verbose = FALSE)
  fit_ghst <- mix.GHST.MoE.EM(y, x, r, g, verbose = FALSE)
  fit_nig  <- mix.NIG.MoE.EM(y, x, r, g,  verbose = FALSE)
  fit_sl   <- mix.SL.MoE.EM(y, x, r, g,  verbose = FALSE)
  fit_nmvbs<- mix.NMVBS.MoE.EM(y, x, r, g, verbose = FALSE)
  
  CAIC_values[i, ] <- c(fit_norm$CAIC,
                        fit_sn$CAIC,
                        fit_ghst$CAIC,
                        fit_nig$CAIC,
                        fit_sl$CAIC,
                        fit_nmvbs$CAIC)
}

CAIC_table <- data.frame(G = G_list, CAIC_values)
write.csv(CAIC_table, "results/CAIC_table.csv", row.names = FALSE)
print(CAIC_table)

#================= Fit final model (best G) =================
best_G <- 2 

fit_norm <- mix.reg.norm.EM(y, x, r, best_G,  verbose = FALSE)
fit_sn   <- mix.reg.SMSN.EM(y, x, r, best_G, family = "Skew.n",  verbose = FALSE)
fit_ghst <- mix.GHST.MoE.EM(y, x, r, best_G,  verbose = FALSE)
fit_nig  <- mix.NIG.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_sl   <- mix.SL.MoE.EM(y, x, r, best_G, verbose = FALSE)
fit_nmvbs<- mix.NMVBS.MoE.EM(y, x, r, best_G,  verbose = FALSE)

#================= Model comparison =================
model_results <- data.frame(
  Model = c("Gaussian", "SN", "GHST", "NIG", "SL", "NMVBS"),
  LogLik = c(fit_norm$loglik, fit_sn$loglik, fit_ghst$loglik, fit_nig$loglik, fit_sl$loglik, fit_nmvbs$loglik),
  CAIC   = c(fit_norm$CAIC, fit_sn$CAIC, fit_ghst$CAIC, fit_nig$CAIC, fit_sl$CAIC, fit_nmvbs$CAIC)
)

write.csv(model_results, "results/model_comparison.csv", row.names = FALSE)
print(model_results)

#================= Choose best model =================
best_model <- fit_nig  
tone_dataset$cluster <- as.factor(best_model$clusters)

#================= Scatter plot with clusters =================
p_cluster <- ggplot(tone_dataset, aes(x = stretchratio, y = tuned, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  labs(x = "Actual Tone Ratio", y = "Perceived Tone Ratio", color = "Cluster")

ggsave("figures/02_cluster_plot.pdf", p_cluster, width = 6, height = 5)
ggsave("figures/02_cluster_plot.png", p_cluster, dpi = 300)

#================= Regression lines =================
beta <- best_model$beta
p_reg <- ggplot(tone_dataset, aes(x = stretchratio, y = tuned, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_abline(intercept = beta[1,1], slope = beta[2,1], color = "red", size = 1.2) +
  geom_abline(intercept = beta[1,2], slope = beta[2,2], color = "blue", size = 1.2) +
  theme_minimal() +
  labs(x = "Actual Tone Ratio", y = "Perceived Tone Ratio")

ggsave("figures/03_regression_lines.pdf", p_reg, width = 6, height = 5)
ggsave("figures/03_regression_lines.png", p_reg, dpi = 300)

#================= Posterior probabilities & uncertainty =================
posterior_prob <- best_model$tau
uncertainty <- 1 - apply(posterior_prob, 1, max)
tone_dataset$uncertainty <- uncertainty

p_uncertainty <- ggplot(tone_dataset, aes(x = stretchratio, y = uncertainty)) +
  geom_point(alpha = 0.7, size = 2.5) +
  theme_minimal() +
  labs(y = "Uncertainty")

ggsave("figures/04_uncertainty_plot.pdf", p_uncertainty, width = 6, height = 5)
ggsave("figures/04_uncertainty_plot.png", p_uncertainty, dpi = 300)

#================= Prediction & RMSE =================
predict_mix <- function(x, beta, group){
  yhat <- numeric(nrow(x))
  for(i in 1:nrow(x)){
    j <- group[i]
    yhat[i] <- x[i,] %*% beta[,j]
  }
  return(yhat)
}

y_pred <- predict_mix(x, beta, best_model$clusters)
rmse <- sqrt(mean((y - y_pred)^2))
cat("RMSE =", rmse, "\n")

#================= Residual plot =================
residuals <- y - y_pred
df_res <- data.frame(index = 1:length(residuals), res = residuals)

p_residual <- ggplot(df_res, aes(x = index, y = res)) +
  geom_point(size = 2.5, alpha = 0.7, color = "#2C3E50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#D62728") +
  theme_minimal() +
  labs(x = "Observation index", y = "Residuals")

ggsave("figures/05_residual_plot.pdf", p_residual, width = 6.5, height = 5)
ggsave("figures/05_residual_plot.png", p_residual, width = 6.5, height = 5, dpi = 300)

#================= End of analysis =================
cat("Analysis completed. Plots saved in /figures, tables in /results.\n")
cat("RMSE:", rmse, "\n")
