rm(list = objects() )
# ---- Load packages ----
install.packages("gstat")
install.packages("eurostat")
install.packages("giscoR")

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
#install.packages("writexl")
library(writexl)
#install.packages("readxl")
library(readxl)
#install.packages("broom.mixed")
library(broom.mixed)

setwd("~/Desktop/Funding 2025/Deliverables/untitled folder/Apenteng")

# ---- 1. Download & unzip shapefile ----
nuts_sf <- st_read("NUTS_RG_01M_2021_4326.shp/NUTS_RG_01M_2021_4326.shp")
# ---- 2. Filter to NUTS2 ----
nuts2_sf <- nuts_sf %>% filter(LEVL_CODE == 2)
nuts2 = nuts2_sf

ASF_nuts2_outbreaks <- readRDS("ASFf_nuts2_outbreaks.rds")
dim(ASF_nuts2_outbreaks)
head(ASF_nuts2_outbreaks)

#-------------------------------------------
# pretest

ASF_nuts2_outbreaks <- ASF_nuts2_outbreaks %>%
  mutate(across(c(Indoor_Mean_value, Outdoor_Mean_value),
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# ---- 4. Filter complete cases ----
ASF_model_data <- ASF_nuts2_outbreaks %>%
  filter(!is.na(outbreak_count))

# Drop geometry for modeling
ASF_model_data_no_geom <- ASF_model_data %>% st_drop_geometry()

# ---- 5. Create binary outcome ----
ASF_model_data_no_geom <- ASF_model_data_no_geom %>%
  mutate(outbreak_binary = ifelse(outbreak_count > 0, 1, 0),
         outbreak_status = ifelse(outbreak_binary == 1, "Outbreak", "No Outbreak"))

# ---- 6. Balance cases and controls ----
cases <- ASF_model_data_no_geom %>% filter(outbreak_binary == 1)
controls_all <- ASF_model_data_no_geom %>% filter(outbreak_binary == 0)

n_controls <- 16
controls <- if (nrow(controls_all) >= n_controls) {
  sample_n(controls_all, n_controls)
} else {
  sample_n(controls_all, n_controls, replace = TRUE)
}
ASF_model_data_balanced <- bind_rows(cases, controls_all)

logit_model1 <- glm(
  outbreak_binary ~  Indoor_Mean_value + Outdoor_Mean_value+ pig_density + land_cover_code,
  data = ASF_model_data_balanced,
  family = binomial(link = "logit")
)

# Step 7: Get odds ratios and 95% confidence intervals
#install.packages("broom.mixed")
odds_ratio_table_case <- broom.mixed::tidy(
  logit_model1,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

summary(logit_model1)
print(odds_ratio_table_case)

#--------------------------------------------
# create biosecurity risk maps
ASF_nuts2_outbreaks1 <- readRDS("ASFf_nuts2_outbreaks.rds")
ASF_nuts2_outbreaks1 <- ASF_nuts2_outbreaks1 %>%
  rowwise() %>%
  mutate(
    Combine_Biosecurity = if (!is.na(Indoor_Mean_value) & !is.na(Outdoor_Mean_value)) {
      (Indoor_Mean_value + Outdoor_Mean_value) / 2
    } else {
      coalesce(Indoor_Mean_value, Outdoor_Mean_value)
    }
  ) %>%
  ungroup()

library(ggplot2)
library(dplyr)
library(sf)

# Ensure outbreak_binary exists

ASF_nuts2_outbreaks1 <- ASF_nuts2_outbreaks1 %>%
  mutate(outbreak_binary = ifelse(outbreak_status == "Outbreak", 1, 0))


logit_model <- glm(
  outbreak_binary ~ Combine_Biosecurity,
  data = ASF_nuts2_outbreaks1,
  family = binomial(link = "logit")
)

# Predict probabilities for each NUTS2 region

ASF_nuts2_outbreaks1 <- ASF_nuts2_outbreaks1 %>%
  mutate(predicted_prob = predict(logit_model, newdata = ., type = "response"))

# Plot EU map with probabilities
ggplot() +
  geom_sf(data = ASF_nuts2_outbreaks1, aes(fill = predicted_prob), color = "white") +
  scale_fill_viridis_c(option = "inferno", name = "Outbreak Probability") +
  coord_sf(xlim = c(-25, 45), ylim = c(34, 72)) +  # EU extent
  theme_minimal()

#------------------------------------------------
# Increase and decrease by 20% in biosecurity mean value 
library(ggplot2)
library(dplyr)
library(sf)
library(gridExtra)

# Ensure outbreak_binary exists

ASF_nuts2_outbreaks1 <- ASF_nuts2_outbreaks1 %>%
  mutate(outbreak_binary = ifelse(outbreak_status == "Outbreak", 1, 0))

# Fit logistic regression

logit_model <- glm(
  outbreak_binary ~ Combine_Biosecurity,
  data = ASF_nuts2_outbreaks1,
  family = binomial(link = "logit")
)

# Predict probabilities with +20% and -20% adjustments without losing geometry

ASF_up <- ASF_nuts2_outbreaks1 %>%
  mutate(Combine_Biosecurity_adj = Combine_Biosecurity * 1.2,
         predicted_prob = predict(logit_model, newdata = data.frame(Combine_Biosecurity = Combine_Biosecurity * 1.2), type = "response"))

ASF_down <- ASF_nuts2_outbreaks1 %>%
  mutate(Combine_Biosecurity_adj = Combine_Biosecurity * 0.8,
         predicted_prob = predict(logit_model, newdata = data.frame(Combine_Biosecurity = Combine_Biosecurity * 0.8), type = "response"))

# Make sure the objects are still sf

ASF_up <- st_as_sf(ASF_up)
ASF_down <- st_as_sf(ASF_down)

# Set common color scale

prob_limits <- range(c(ASF_up$predicted_prob, ASF_down$predicted_prob), na.rm = TRUE)

# Plot maps

p_up <- ggplot(ASF_up) +
  geom_sf(aes(fill = predicted_prob), color = "white") +
  scale_fill_viridis_c(option = "inferno", limits = prob_limits, name = "Outbreak Probability") +
  coord_sf(xlim = c(-25, 45), ylim = c(34, 72)) +
  theme_minimal()

p_down <- ggplot(ASF_down) +
  geom_sf(aes(fill = predicted_prob), color = "white") +
  scale_fill_viridis_c(option = "inferno", limits = prob_limits, name = "Outbreak Probability") +
  coord_sf(xlim = c(-25, 45), ylim = c(34, 72)) +
  theme_minimal()

grid.arrange(p_up, p_down, ncol = 2)

