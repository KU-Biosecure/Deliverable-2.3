rm(list = objects() )
# ---- Load packages ----
library(sf)
library(tidyverse)
library(ggplot2)
library(gganimate)
library(gifski)
library(tmap)
library(readr)
library(httr)
library(utils)
library(gstat)
library(lattice)
library(readxl)
library(sp)
library(terra)
library(dplyr)
library(eurostat)
library(raster)
library(giscoR)
library(stringr)
library(eurostat)

setwd("~/Desktop/Funding 2025/Deliverables/untitled folder/Apenteng")

# ---- 1. Download & unzip shapefile ----
nuts_sf <- st_read("NUTS_RG_01M_2021_4326.shp/NUTS_RG_01M_2021_4326.shp")
# ---- 2. Filter to NUTS2 ----
nuts2_sf <- nuts_sf %>% filter(LEVL_CODE == 2)
nuts2 = nuts2_sf

AI_nuts2_outbreaks<- readRDS("AI_nuts2_outbreaks_poultry.rds")

AI_nuts2_outbreaks<- AI_nuts2_outbreaks %>%
  mutate(
    Free_range_layers_Mean_value = ifelse(is.na(Free_range_layers_Mean_value),
                                          mean(Free_range_layers_Mean_value, na.rm = TRUE), Free_range_layers_Mean_value),
    Free_range_broilers_Mean_value = ifelse(is.na(Free_range_broilers_Mean_value),
                                            mean(Free_range_broilers_Mean_value, na.rm = TRUE), Free_range_broilers_Mean_value),
    Broilers_Mean_value = ifelse(is.na(Broilers_Mean_value),
                                 mean(Broilers_Mean_value, na.rm = TRUE), Broilers_Mean_value),  
    Laying_hens_Mean_value = ifelse(is.na(Laying_hens_Mean_value),
                                    mean(Laying_hens_Mean_value, na.rm = TRUE), Laying_hens_Mean_value) 
  )

AI_model_data <- AI_nuts2_outbreaks %>%
  filter(!is.na(outbreak_count))

AI_model_data_baseline <- AI_model_data %>%
  mutate(outbreak_binary = ifelse(outbreak_count > 0, 1, 0))

cases <- AI_model_data_baseline %>% filter(outbreak_binary == 1)
controls_all <- AI_model_data_baseline %>% filter(outbreak_binary == 0)

n_controls <- 82
controls <- if (nrow(controls_all) >= n_controls) {
  sample_n(controls_all, n_controls)
} else {
  sample_n(controls_all, n_controls, replace = TRUE)
}
AI_model_data_balanced_baseline<- bind_rows(cases, controls_all)

AI_model_data_balanced_baseline <- AI_model_data_balanced_baseline %>%
  mutate(grouped_land_cover = case_when(
    land_cover_name %in% c("Arable land", "Pastures", "Permanent crops") ~ "Agricultural",
    land_cover_name %in% c("Forest", "Scrub and/or herbaceous vegetation associations", "Open spaces with little or no vegetation") ~ "Natural",
    land_cover_name %in% c("Urban fabric", "Industrial units", "Green urban areas", "Sport and leisure facilities") ~ "Urban",
    land_cover_name %in% c("Wetlands", "Water bodies") ~ "Water-related",
    TRUE ~ NA_character_
  )) %>%
  mutate(grouped_land_cover = factor(grouped_land_cover))

# Step 6: Fit logistic regression model
logit_model_baseline <- glm(
  outbreak_binary ~ poultry_density + grouped_land_cover + Broilers_Mean_value  + Laying_hens_Mean_value + Free_range_layers_Mean_value + Free_range_broilers_Mean_value,
  data = AI_model_data_balanced_baseline,
  family = binomial(link = "logit")
)

# Step 7: Get odds ratios and 95% confidence intervals
odds_ratio_table_case_baseline <- broom.mixed::tidy(
  logit_model_baseline,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

# Step 8: View results
summary(logit_model_baseline)
print(odds_ratio_table_case_baseline)

AI_model_data_balanced_baseline$risk_score <- predict(logit_model_baseline, newdata = AI_model_data_balanced_baseline, type = "response")

# Step 10: Extract coordinates
AI_model_data_balanced_baseline <- AI_model_data_balanced_baseline %>%
  mutate(
    lon = st_coordinates(st_centroid(geometry))[, 1],
    lat = st_coordinates(st_centroid(geometry))[, 2]
  ) %>%
  st_drop_geometry()

# Step 11: Convert to sf and spatial join
AI_sf_baseline <- st_as_sf(AI_model_data_balanced_baseline, coords = c("lon", "lat"), crs = 4326)
AI_joined_baseline <- st_join(AI_sf_baseline, nuts2, join = st_intersects)

# Step 12: Aggregate risk scores by NUTS2 region
nuts2_risk_baseline <- AI_joined_baseline %>%
  st_drop_geometry() %>%
  group_by(NUTS_ID.y) %>%
  summarise(mean_risk = mean(risk_score, na.rm = TRUE))

nuts2_map_baseline <- left_join(nuts2, nuts2_risk_baseline, by = c("NUTS_ID" = "NUTS_ID.y"))

# Step 13: Identify regions with outbreak OR control AND risk score
regions_with_condition <- AI_joined_baseline %>%
  st_drop_geometry() %>%
  group_by(NUTS_ID.y) %>%
  summarise(has_outbreak = any(outbreak_binary == 1),
            has_control = any(outbreak_binary == 0),
            risk_available = any(!is.na(risk_score))) %>%
  filter((has_outbreak | has_control) & risk_available) %>%
  pull(NUTS_ID.y)
regions_with_condition <- regions_with_condition[regions_with_condition != "MT00"]

# Step 15: Filter map and AI points
nuts2_filtered_baseline <- nuts2_map_baseline %>%
  filter(NUTS_ID %in% regions_with_condition & is.finite(mean_risk))

AI_filtered_baseline <- AI_joined_baseline %>%
  filter(NUTS_ID.y %in% nuts2_filtered_baseline$NUTS_ID & !is.na(risk_score))  # ✅ Filter here

# Step 16: Plot (only points with valid risk scores)
b_equl_baseline <- ggplot() +
  geom_sf(data = nuts2_map_baseline, fill = "lightgrey", color = "white") +
  geom_sf(data = nuts2_filtered_baseline, aes(fill = mean_risk), color = "white") +
  # geom_sf(data = AI_filtered_baseline |> filter(outbreak_binary == 1), 
  #        aes(color = "Outbreak Cases"), shape = 16, size = 2) +
  # geom_sf(data = AI_filtered_baseline |> filter(outbreak_binary == 0), 
  #         aes(color = "Control Regions"), shape = 16, size = 2) +
  scale_fill_viridis_c(name = "Risk scores", option = "inferno") +
  # scale_color_manual(
  #   name = "Outbreak Classification",
  #   values = c("Outbreak Cases" = "red", "Control Regions" = "blue")
  # ) +
  coord_sf(xlim = c(-10, 55), ylim = c(35, 70)) +
  theme_minimal() +
  labs(
    title = "HPAI Risk Map (Default in Broilers & Laying Hens)", #replace here see below
    caption = "Grey = all NUTS2 regions; colored = regions with outbreak OR control and valid risk score."
  )
#"HPAI Risk Map (Regions with Outbreak OR Control & Risk)"
#"HPAI Risk Map (20% Increase in Broilers & Laying Hens)"
#"HPAI Risk Map (20% Decrease in Broilers & Laying Hens)"

print(b_equl_baseline)

# Step 17: Save the plot
ggsave(
  filename = "AI_Risk_Figures/ai_poultry_Default.png", #baseline increase decrease
  plot = b_equl_baseline,
  width = 10,
  height = 8,
  dpi = 300
)

#---------------
# create biosecurity risk maps
df<- readRDS("AI_nuts2_outbreaks_poultry.rds")
df <- df %>%
  rowwise() %>%
  mutate(
    Combine_Biosecurity = mean(c_across(
      c(
        Free_range_layers_Mean_value,
        Free_range_broilers_Mean_value,
        Broilers_Mean_value,
        Laying_hens_Mean_value
      )
    ), na.rm = TRUE),
    Combine_Biosecurity = ifelse(is.nan(Combine_Biosecurity), NA, Combine_Biosecurity)
  ) %>%
  ungroup()

#-------------------------------------------------
# 3. Create outbreak_binary
#-------------------------------------------------
df <- df %>%
  mutate(outbreak_binary = ifelse(outbreak_count > 0, 1, 0))

logit_model <- glm(
  outbreak_binary ~ Combine_Biosecurity,
  data = df,
  family = binomial(link = "logit")
)

df <- df %>%
  mutate(predicted_prob = predict(logit_model, newdata = ., type = "response"))
df_no_mt <- df %>% filter(NUTS_ID != "MT00")
p_default <- ggplot(df_no_mt ) +
  geom_sf(aes(fill = predicted_prob), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "inferno", name = "Outbreak Probability") +
  coord_sf(xlim = c(-25, 45), ylim = c(34, 72)) +  # EU extent
  theme_minimal()

library(sf)
library(ggplot2)
library(dplyr)

# Transform df to WGS84 (EPSG:4326) if not already

df <- st_transform(df, 4326)

# Define Europe bounding box in WGS84

bbox_europe <- st_sfc(st_polygon(list(rbind(
  c(-25, 34),   # bottom-left
  c(45, 34),    # bottom-right
  c(45, 72),    # top-right
  c(-25, 72),   # top-left
  c(-25, 34)    # close polygon
))), crs = st_crs(df))

# Filter only regions intersecting Europe

df_europe <- df[st_intersects(df, bbox_europe, sparse = FALSE)[,1], ]

# Check that df_europe is not empty

nrow(df_europe)

# Predict probabilities

df_europe <- df_europe %>%
  mutate(predicted_prob = predict(logit_model, newdata = ., type = "response"))

# Plot

ggplot(df_europe) +
  geom_sf(aes(fill = predicted_prob), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "inferno", name = "Outbreak Probability", limits = c(0, 0.5)) +
  coord_sf(xlim = c(-25, 45), ylim = c(34, 72)) +
  theme_minimal()

#-------------------------------------------------
# 5. +20% and -20% Biosecuirty prediction risk maps
#-------------------------------------------------
library(dplyr)
library(sf)
library(ggplot2)
library(gridExtra)
library(grid)

# Define Europe bounding box (approximate)

bbox_europe <- st_sfc(st_polygon(list(rbind(
  c(-25, 34), c(45, 34), c(45, 72), c(-25, 72), c(-25, 34)
))), crs = st_crs(df))  # use CRS of your spatial data

# Keep only regions intersecting Europe

df_europe <- df[st_intersects(df, bbox_europe, sparse = FALSE)[,1], ]

#-----------------------------------------------------------
# Calculate +20% and -20% biosecurity predicted probabilities

df_up <- df_europe %>%
  mutate(predicted_prob = predict(
    logit_model,
    newdata = data.frame(Combine_Biosecurity = Combine_Biosecurity * 1.2),
    type = "response"
  ))

df_down <- df_europe %>%
  mutate(predicted_prob = predict(
    logit_model,
    newdata = data.frame(Combine_Biosecurity = Combine_Biosecurity * 0.8),
    type = "response"
  ))

# Ensure they are sf objects

df_up <- st_as_sf(df_up)
df_down <- st_as_sf(df_down)

# Shared color scale

prob_limits <- c(0, 0.5)

# Create plots

p_minus20 <- ggplot(df_down) +
  geom_sf(aes(fill = predicted_prob), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "inferno", limits = prob_limits, name = "Outbreak Probability") +
  theme_minimal() +
  labs(title = "a) -20% Biosecurity")

p_plus20 <- ggplot(df_up) +
  geom_sf(aes(fill = predicted_prob), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "inferno", limits = prob_limits, name = "Outbreak Probability") +
  theme_minimal() +
  labs(title = "b) +20% Biosecurity")

# Optional: add labels above plots

label_a <- textGrob("a) -20% Biosecurity", gp = gpar(fontsize = 14, fontface = "bold"))
label_b <- textGrob("b) +20% Biosecurity", gp = gpar(fontsize = 14, fontface = "bold"))

# Arrange side by side

grid.arrange(
  arrangeGrob(label_a, p_minus20, ncol = 1, heights = c(0.1, 0.9)),
  arrangeGrob(label_b, p_plus20, ncol = 1, heights = c(0.1, 0.9)),
  ncol = 2
)



