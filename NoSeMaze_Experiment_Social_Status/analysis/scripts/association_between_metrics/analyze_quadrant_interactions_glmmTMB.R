# ============================================================
# Script: analyze_quadrant_interactions_glmmTMB.R
# Purpose: Analyze quadrant-level interaction metrics 
#          using beta regression mixed models (glmmTMB)
# Author: [Your Name]
# Date: [Insert Date]
# ============================================================

# 0. Load required packages
library(glmmTMB)     # For GLMMs with flexible distributions
library(readr)       # For reading CSV files
library(DHARMa)      # For model diagnostics
library(emmeans)     # For estimated marginal means
library(ggplot2)     # For plotting

# 1. Load and prepare the data
data <- read_csv("myRootPath/NoSeMaze_Experiment/results/figures/cross_cohort/associations_between_metrics/tube/match_matrices/interaction_data.csv")
# --> CHANGE TO YOUR ROOT PATH

# Convert relevant variables to factors
data$quadrant <- as.factor(data$quadrant)
data$group_id <- as.factor(data$group_id)

# 2. Adjust values exactly at 0 or 1 (required for beta distribution)
epsilon <- 1e-4
data$ChasByComp <- pmin(pmax(data$ChasByComp, epsilon), 1 - epsilon)
data$ChasByChas <- pmin(pmax(data$ChasByChas, epsilon), 1 - epsilon)
data$CompByComp <- pmin(pmax(data$CompByComp, epsilon), 1 - epsilon)

# 3. Fit GLMM with beta distribution and logit link
model <- glmmTMB(
  ChasByComp ~ quadrant + (1 | group_id),
  data = data,
  family = beta_family(link = "logit")
)

# Summarize model results
summary(model)

# 4. Estimate marginal means on response scale (0–1)
marginal <- emmeans(model, ~ quadrant, type = "response")
summary(marginal)

# 5. Pairwise comparisons between quadrants
emmeans(model, pairwise ~ quadrant, type = "response")

# 6. Model diagnostics using DHARMa
sim_res <- simulateResiduals(model)
plot(sim_res)

# 7. Plot estimated marginal means with confidence intervals
df_emm <- as.data.frame(marginal)

ggplot(df_emm, aes(x = quadrant, y = response)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  ylab("Estimated Mean Interaction") +
  xlab("Quadrant") +
  ggtitle("Estimated Mean Interaction by Quadrant (Beta GLMM)") +
  theme_minimal()

# Save plot
ggsave("interaction_plot.png", width = 6, height = 4, dpi = 300)

# ============================================================
# OPTIONAL: Model Variations and Alternative Specifications
# ============================================================

# Alternative model using different link function (probit)
model_alt <- glmmTMB(
  interaction ~ quadrant + (1 | group_id),
  data = data,
  family = beta_family(link = "probit")
)

# Diagnostics for alternative model
sim_res_alt <- simulateResiduals(model_alt)
plot(sim_res_alt)

# Random slope model: allow quadrant effect to vary by group
model_slope <- glmmTMB(
  interaction ~ quadrant + (quadrant | group_id),
  data = data,
  family = beta_family(link = "logit")
)
summary(model_slope)

# Diagnostics for random slope model
sim_res_slope <- simulateResiduals(model_slope)
plot(sim_res_slope)

# Alternative link: log
model_log <- glmmTMB(
  interaction ~ quadrant + (1 | group_id),
  data = data,
  family = beta_family(link = "log")
)