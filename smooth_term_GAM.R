# Load required packages
library(mgcv)         # For GAM modeling
library(ggplot2)      # For visualizations
library(patchwork)    # For combining plots
library(gratia)       # For GAM diagnostics and visualization
library(MuMIn)        # For model selection metrics (e.g., AIC)
library(readr)        # For reading CSV files
library(dplyr)        # For data manipulation
library(sf)           # For spatial data handling

# -----------------------------
# Load and prepare the dataset
# -----------------------------
recovery <- read_csv("/recovery.csv")

# Filter: only use disturbances before 2013, create binary variable for recovery success (<=10 years)
recovery_filtered <- recovery %>%
  group_by(ID) %>%
  filter(yod < 2013) %>%
  mutate(recovered_10y = ifelse(recovery_rate <= 10, 1, 0)) %>%
  ungroup()

# Compute pre-disturbance forest composition and post-disturbance bare ground
recovery_filtered <- recovery_filtered %>%
  group_by(ID) %>% 
  mutate(
    pre_coniferous = ifelse(year < yod, mean(coniferous[year < yod], na.rm = TRUE), NA),
    pre_broadleaved = ifelse(year < yod, mean(broadleaved[year < yod], na.rm = TRUE), NA),
    post_bare_ground = ifelse(year < yod, mean(bare_ground[year > yod], na.rm = TRUE), NA)
  ) %>%
  ungroup() %>%
  mutate(
    pre_dist_tree_cover = pre_coniferous + pre_broadleaved
  )

# Create unique entry per disturbance event (per ID)
disturbance_summary <- recovery_filtered %>%
  distinct(ID, .keep_all = TRUE)

# Convert to spatial features (sf)
disturbance_sf <- st_as_sf(disturbance_summary, coords = c("x", "y"), crs = 3035)

# Load hexagonal grid and keep only GRID_ID
hex_grid <- st_read("hex_500.shp") %>% select(GRID_ID)

# Spatially join recovery data to hex grid
hex_joined <- st_join(disturbance_sf, hex_grid, join = st_intersects)

# -----------------------------
# Aggregate predictor variables
# -----------------------------
hex_summary <- hex_joined %>%
  group_by(GRID_ID) %>%
  summarise(
    elevation = mean(height, na.rm = TRUE),
    severity = mean(severity_relative, na.rm = TRUE),
    VPD_yod1 = mean(VPD_yod1, na.rm = TRUE),
    temp_10y = mean(mean_temp10, na.rm = TRUE),
    prec_10y = mean(mean_prec10, na.rm = TRUE),
    temp_total = mean(mean_temp_total, na.rm = TRUE),
    prec_total = mean(mean_prec_total, na.rm = TRUE),
    recovery_rate = mean(recovery_rate, na.rm = TRUE),
    recovery_success = mean(recovered_10y, na.rm = TRUE) * 100,
    broadleaved = mean(pre_broadleaved, na.rm = TRUE),
    coniferous = mean(pre_coniferous, na.rm = TRUE),
    bare_ground = mean(post_bare_ground, na.rm = TRUE),
    forest_type = names(sort(table(forest_type), decreasing = TRUE))[1],
    geolocation = names(sort(table(geoloc), decreasing = TRUE))[1],
    .groups = "drop"
  )

# Join back to spatial grid
hex_model_data <- st_join(hex_grid, hex_summary, join = st_intersects)

# Remove incomplete entries
hex_model_data <- hex_model_data %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  mutate(pre_tree_cover = coniferous + broadleaved)

# -----------------------------
# Prepare data for modeling
# -----------------------------
hex_model_data <- hex_model_data %>%
  mutate(geolocation = as.factor(geolocation))

# Compute centroid coordinates
hex_model_data$centroid <- st_centroid(hex_model_data$geometry)
hex_model_data$long <- st_coordinates(hex_model_data$centroid)[,1]
hex_model_data$lat <- st_coordinates(hex_model_data$centroid)[,2]

# -----------------------------
# Fit Generalized Additive Model
# -----------------------------
fit_gam <- gam(
  recovery_success ~ 
    s(long, lat, bs = "tp") +
    s(severity) +
    s(VPD_yod1, by = geolocation) +
    s(temp_total) +
    s(prec_total) +
    s(elevation) +
    s(pre_tree_cover) +
    s(bare_ground),
  data = hex_model_data, method = "REML"
)

# Model summary
summary(fit_gam)

# Model metrics
cat("Model Metrics:\n",
    "AIC:", AIC(fit_gam), "\n",
    "Deviance explained:", round(summary(fit_gam)$dev.expl, 3), "\n",
    "Adjusted R²:", round(summary(fit_gam)$r.sq, 3), "\n",
    "GCV:", round(fit_gam$gcv.ubre, 3), "\n",
    "N:", fit_gam$n, "\n")

# Extract partial effects
partial_effects <- predict(fit_gam, type = "terms", se.fit = TRUE)

# Add VPD effect to data for mapping
hex_model_data$VPD_effect <- partial_effects$fit[, "s(VPD_yod1)"]

# Residuals
hex_model_data$residuals <- residuals(fit_gam)
hex_model_data$predicted <- fitted(fit_gam)

# -----------------------------
# Diagnostic plots
# -----------------------------
appraise(fit_gam)

# -----------------------------
# Predictors effect plots
# -----------------------------
predictors <- c("elevation", "severity", "VPD_yod1", "prec_total", "temp_total", "pre_tree_cover", "bare_ground")
plot_data <- hex_model_data %>%
  select(predicted, all_of(predictors)) %>%
  pivot_longer(-predicted, names_to = "predictor", values_to = "value") %>%
  mutate(predictor = recode(predictor,
                            VPD_yod1 = "VPD anomalies",
                            temp_total = "Temperature",
                            prec_total = "Precipitation",
                            severity = "Severity",
                            pre_tree_cover = "Pre-disturbance\ntree cover",
                            bare_ground = "Post-disturbance\nbare ground share",
                            elevation = "Elevation"
  )) %>%
  mutate(predictor = factor(predictor, levels = c("VPD anomalies", "Elevation", "Severity", "Temperature",
                                                  "Precipitation", "Pre-disturbance\ntree cover",
                                                  "Post-disturbance\nbare ground share")))

ggplot(plot_data, aes(x = value, y = predicted)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp"), color = "#11828A") +
  facet_wrap(~ predictor, scales = "free_x", nrow = 2) +
  theme_bw(base_size = 18) +
  labs(y = "Predicted recovery success", x = "Predictor value")


# -----------------------------
# Map predictions and residuals
# -----------------------------
ggplot(hex_model_data) +
  geom_sf(aes(fill = predicted), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "magma", name = "Predicted Recovery", trans = "sqrt") +
  theme_bw(base_size = 22) +
  theme(legend.position = "right")


# Residual maps
ggplot(hex_model_data) +
  geom_sf(aes(fill = residuals)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_bw(base_size = 22)


# Effect map for VPD
ggplot(hex_model_data) +
  geom_sf(aes(fill = VPD_effect), color = "grey") +
  scale_fill_gradient2(low = "#126D0E", mid = "white", high = "#E69E03", midpoint = 0) +
  theme_bw(base_size = 22)
