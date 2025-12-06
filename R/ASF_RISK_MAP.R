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
library(writexl)
library(readxl)
library(broom.mixed)

setwd("~/Desktop/BIOSECURE_FINAL_RISK_SCORE")


# ---- 1. Download & unzip shapefile ----
nuts_sf <- st_read("NUTS_RG_01M_2021_4326.shp/NUTS_RG_01M_2021_4326.shp")
# ---- 2. Filter to NUTS2 ----
nuts2_sf <- nuts_sf %>% filter(LEVL_CODE == 2)
nuts2 = nuts2_sf

ASF_nuts2_outbreaks <- readRDS("ASFf_nuts2_outbreaks.rds")

dim(ASF_nuts2_outbreaks)
head(ASF_nuts2_outbreaks)

library(dplyr)


ASF_nuts2_outbreaks <- ASF_nuts2_outbreaks %>%
  mutate(
    Indoor_Mean_value = ifelse(is.na(Indoor_Mean_value), mean(Indoor_Mean_value, na.rm = TRUE), Indoor_Mean_value)) 
#Outdoor_Mean_value = ifelse(is.na(Outdoor_Mean_value), mean(Outdoor_Mean_value, na.rm = TRUE), Outdoor_Mean_value)
# )

#-------------------------------------------------
# Scenario Decrease Biosecurity measures but glm name Indoor_Mean_value_25up have to change as well using this part

#ASF_nuts2_outbreaks <- ASF_nuts2_outbreaks %>%
#  mutate(
#    Indoor_Mean_value_25up = ifelse(
#      outbreak_status == "Outbreak",
#     Indoor_Mean_value * 0.8,
#      Indoor_Mean_value
#    )
#  )
#-------------------------------------------------
ASF_nuts2_outbreaks %>%
  group_by(outbreak_status) %>%
  summarise(mean_new = mean(Indoor_Mean_value, na.rm = TRUE))

# ---- 4. Filter complete cases ----
ASF_model_data <- ASF_nuts2_outbreaks %>%
  filter(!is.na(outbreak_count),
         #!is.na(pig_density),
         #!is.na(land_cover_name),
         !is.na(Indoor_Mean_value))#,
#!is.na(Outdoor_Mean_value))

# Drop geometry for modeling
ASF_model_data_no_geom <- ASF_model_data %>% st_drop_geometry()

# ---- 5. Create binary outcome ----
ASF_model_data_no_geom <- ASF_model_data_no_geom %>%
  mutate(outbreak_binary = ifelse(outbreak_count > 0, 1, 0),
         outbreak_status = ifelse(outbreak_binary == 1, "Outbreak", "No Outbreak"))

# ---- 6. Balance cases and controls ----
cases <- ASF_model_data_no_geom %>% filter(outbreak_binary == 1)
controls_all <- ASF_model_data_no_geom %>% filter(outbreak_binary == 0)

#set.seed(123)
n_controls <- 16
controls <- if (nrow(controls_all) >= n_controls) {
  sample_n(controls_all, n_controls)
} else {
  sample_n(controls_all, n_controls, replace = TRUE)
}
ASF_model_data_balanced <- bind_rows(cases, controls_all)

# ---- 7. Group land cover ----
ASF_model_data_balanced <- ASF_model_data_balanced %>%
  mutate(grouped_land_cover = case_when(
    land_cover_name %in% c("Arable land", "Pastures", "Permanent crops") ~ "Agricultural",
    land_cover_name %in% c("Forest", "Scrub and/or herbaceous vegetation associations", "Open spaces with little or no vegetation") ~ "Natural",
    land_cover_name %in% c("Urban fabric", "Industrial units", "Green urban areas", "Sport and leisure facilities") ~ "Urban",
    land_cover_name %in% c("Wetlands", "Water bodies") ~ "Water-related",
    TRUE ~ NA_character_
  )) %>%
  mutate(grouped_land_cover = factor(grouped_land_cover))

# ---- 8. Fit logistic regression ----
logit_model <- glm(
  outbreak_binary ~  Indoor_Mean_value+grouped_land_cover+pig_density,
  data = ASF_model_data_balanced,
  family = binomial(link = "logit")
)

# Step 7: Get odds ratios and 95% confidence intervals
#install.packages("broom.mixed")
odds_ratio_table_case <- broom.mixed::tidy(
  logit_model,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

# Step 8: View results
summary(logit_model)
print(odds_ratio_table_case)

# ---- 9. Predict risk scores for full dataset ----
ASF_model_data <- ASF_model_data %>%
  mutate(grouped_land_cover = case_when(
    land_cover_name %in% c("Arable land", "Pastures", "Permanent crops") ~ "Agricultural",
    land_cover_name %in% c("Forest", "Scrub and/or herbaceous vegetation associations", "Open spaces with little or no vegetation") ~ "Natural",
    land_cover_name %in% c("Urban fabric", "Industrial units", "Green urban areas", "Sport and leisure facilities") ~ "Urban",
    land_cover_name %in% c("Wetlands", "Water bodies") ~ "Water-related",
    TRUE ~ NA_character_
  )) %>%
  mutate(grouped_land_cover = factor(grouped_land_cover))

set.seed(123)
ASF_model_data$risk_score <- predict(logit_model, newdata = ASF_model_data, type = "response")

# ---- 10. Extract coordinates for full dataset ----
ASF_model_data <- ASF_model_data %>%
  mutate(
    lon = st_coordinates(st_centroid(geometry))[, 1],
    lat = st_coordinates(st_centroid(geometry))[, 2]
  )

ASF_model_data <- st_drop_geometry(ASF_model_data)

# ---- 11. Add coordinates to balanced dataset ----
ASF_model_data_balanced <- ASF_model_data_balanced %>%
  mutate(
    lon = ASF_model_data$lon[match(NUTS_ID, ASF_model_data$NUTS_ID)],
    lat = ASF_model_data$lat[match(NUTS_ID, ASF_model_data$NUTS_ID)]
  )

# Convert both to sf
ASF_sf <- st_as_sf(ASF_model_data, coords = c("lon", "lat"), crs = 4326)
ASF_balanced_sf <- st_as_sf(ASF_model_data_balanced, coords = c("lon", "lat"), crs = 4326)

# ---- 12. Join full dataset with NUTS2 ----
ASF_joined <- st_join(ASF_sf, nuts2, join = st_intersects)

# ---- 13. Aggregate risk scores by NUTS2 ----
nuts2_risk <- ASF_joined |>
  st_drop_geometry() |>
  group_by(NUTS_ID.y) |>
  summarise(mean_risk = mean(risk_score, na.rm = TRUE))

# ---- 14. Visualization: Gray polygons for context, valid scores only ----

# Full NUTS2 map for context (gray for missing risk scores)
nuts2_map <- nuts2 %>%
  left_join(nuts2_risk, by = c("NUTS_ID" = "NUTS_ID.y"))

# Filter outbreak/control regions with valid risk scores
nuts2_selected <- nuts2_map %>%
  filter(NUTS_ID %in% ASF_model_data_balanced$NUTS_ID & !is.na(mean_risk))

# Filter points to match these regions
ASF_balanced_filtered <- ASF_balanced_sf %>%
  filter(NUTS_ID %in% nuts2_selected$NUTS_ID)

# Join outbreak status
nuts2_selected <- nuts2_selected %>%
  left_join(ASF_model_data_balanced %>%
              dplyr::select(NUTS_ID, outbreak_status),
            by = "NUTS_ID")

# Plot
p_final <- ggplot() +
  # Gray polygons for all regions without risk scores
  geom_sf(data = nuts2_map, fill = "gray90", color = "white") +
  # Colored polygons for outbreak/control regions with valid risk scores
  geom_sf(data = nuts2_selected, aes(fill = mean_risk), color = "white") +
  # Points for outbreak/control regions
  # geom_sf(data = ASF_balanced_filtered, size = 2) +
  scale_fill_viridis_c(name = " Risk scores", option = "inferno") +
  coord_sf(xlim = c(-25, 45), ylim = c(34, 72)) +
  theme_minimal() +
  labs(
    title = "ASF Risk Map (Default) ", #please remember to set to defaut
    #caption = "Grey = all NUTS2 regions; colored = regions with outbreak OR control and valid risk score."
  )

#ASF Risk Map (20% increase in biosecurity score) 
#ASF Risk Map (20% decrease in biosecurity score)
print(p_final)
