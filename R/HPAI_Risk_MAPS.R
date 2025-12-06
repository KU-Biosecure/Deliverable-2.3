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

setwd("~/Desktop/BIOSECURE_FINAL_RISK_SCORE")

# ---- 1. Download & unzip shapefile ----
nuts_sf <- st_read("NUTS_RG_01M_2021_4326.shp/NUTS_RG_01M_2021_4326.shp")
# ---- 2. Filter to NUTS2 ----
nuts2_sf <- nuts_sf %>% filter(LEVL_CODE == 2)
nuts2 = nuts2_sf

AI_nuts2_outbreaks <- readRDS("AI_nuts2_outbreaks_poultry.rds")

write_xlsx(AI_nuts2_outbreaks, "AI_nuts2_outbreaks_poultry.xlsx")

AI_nuts2_outbreaks <- AI_nuts2_outbreaks %>%
  mutate(outbreak_status = ifelse(outbreak_count == 0,
                                  "No Outbreak",
                                  "Outbreak"))

head(AI_nuts2_outbreaks)

#----  Aalysis start from here-------------------
AI_nuts2_outbreaks<- AI_nuts2_outbreaks %>%
  mutate(
    # Free_range_layers_Mean_value = ifelse(is.na(Free_range_layers_Mean_value),
    #                                       mean(Free_range_layers_Mean_value, na.rm = TRUE), Free_range_layers_Mean_value),
    # Free_range_broilers_Mean_value = ifelse(is.na(Free_range_broilers_Mean_value),
    #                                         mean(Free_range_broilers_Mean_value, na.rm = TRUE), Free_range_broilers_Mean_value),
    Broilers_Mean_value = ifelse(is.na(Broilers_Mean_value),
                                 mean(Broilers_Mean_value, na.rm = TRUE), Broilers_Mean_value),  
    Laying_hens_Mean_value = ifelse(is.na(Laying_hens_Mean_value),
                                    mean(Laying_hens_Mean_value, na.rm = TRUE), Laying_hens_Mean_value) 
  )

#-------------------------------------------------
# Scenario Decrease Biosecurity measures but glm names Indoor_Broilers_25up ND Indoor_Laying_25up haVE to change as well using this part

#AI_nuts2_outbreaks <- AI_nuts2_outbreaks %>%
#  mutate(
#    Indoor_Broilers_25up = ifelse(outbreak_status == "Outbreak",
#                                  Broilers_Mean_value * 0.8,
#                                  Broilers_Mean_value),
    
#    Indoor_Laying_25up = ifelse(outbreak_status == "Outbreak",
 #                               Laying_hens_Mean_value * 0.8,
 #                               Laying_hens_Mean_value)
 # )
#-------------------------------------------------

# Step 2: Filter complete cases
AI_model_data <- AI_nuts2_outbreaks%>%
  filter(!is.na(outbreak_count),
         #!is.na(poultry_density),
         #!is.na(Free_range_layers_Mean_value),
         #!is.na(land_cover_name),
         #!is.na(Free_range_broilers_Mean_value),
         !is.na(Broilers_Mean_value),
         !is.na(Laying_hens_Mean_value))

# Step 3: Create binary outcome
AI_model_data_baseline <- AI_model_data %>%
  mutate(outbreak_binary = ifelse(outbreak_count > 0, 1, 0))

# Step 4: Sample 50 cases and 50 controls
cases <- AI_model_data_baseline %>% filter(outbreak_binary == 1)
controls_all <- AI_model_data_baseline %>% filter(outbreak_binary == 0)

# Beate

n_controls <- 82
controls <- if (nrow(controls_all) >= n_controls) {
  sample_n(controls_all, n_controls)
} else {
  sample_n(controls_all, n_controls, replace = TRUE)
}
AI_model_data_balanced_baseline<- bind_rows(cases, controls_all)

#n_samples <- 82
#sampled_cases <- sample_n(cases, n_samples)
#sampled_controls <- sample_n(controls_all, n_samples, replace = TRUE)
#AI_model_data_balanced_baseline <- bind_rows(sampled_cases, sampled_controls)

# Step 5: Group land cover
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
  outbreak_binary ~ poultry_density + grouped_land_cover + Broilers_Mean_value  + Laying_hens_Mean_value ,
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

# Step 9: Predict risk scores
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

# Step 14: Remove Malta (MT00)
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
