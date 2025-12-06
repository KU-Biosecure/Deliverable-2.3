
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

#HPAI_Poultry data
HPAI_Poultry_Birds <- read_excel("HPAI_Poultry_Birds.xlsx")%>%  
  filter(Epi_unit == "Farm")

AI_Poultry_Birds <- HPAI_Poultry_Birds 
# Convert AI outbreak data to sf object
AI_Birds_sf <- AI_Poultry_Birds %>%
  filter(!is.na(Longitude), !is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(nuts2))

# Spatial join: assign each outbreak to a NUTS2 region
AI_Bird_with_NUTS2 <- st_join(AI_Birds_sf, nuts2, join = st_intersects)

# Aggregate outbreak counts by NUTS2 region
outbreak_counts <- AI_Bird_with_NUTS2 %>%
  st_drop_geometry() %>%
  group_by(NUTS_ID) %>%
  #group_by(NAME_LATN.x) %>%
  summarise(outbreak_count = n())
summary(outbreak_counts)


# Now join outbreak counts to nuts2 geometries
nuts2_outbreaks <- nuts2 %>%
  left_join(outbreak_counts, by = "NUTS_ID") %>%
  mutate(outbreak_count = replace_na(outbreak_count, 0))

# Create a new column to classify regions
nuts2_outbreaks_ai <- nuts2_outbreaks %>%
  mutate(outbreak_status = ifelse(outbreak_count > 0, "Outbreak", "No Outbreak"))


# Plot the map
b_cases <- ggplot(nuts2_outbreaks_ai) +
  geom_sf(aes(fill = outbreak_status), color = "white", size = 0.1) +
  scale_fill_manual(values = c("Outbreak" = "red", "No Outbreak" = "lightgray")) +
  labs(title = "HPAI Bird Outbreak Status by NUTS2 Region",
       fill = "Status") +
  coord_sf(xlim = c(-10, 55), ylim = c(35, 70)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())


# Print to viewer
print(b_cases)

# Create folder if it doesn't exist
if (!dir.exists("AI_Risk_Figures")) {
  dir.create("AI_Risk_Figures")
}

# Save the plot to the folder
ggsave(
  filename = "AI_Risk_Figures/asf_bird_cases_map.png",
  plot = b_cases,
  width = 10,
  height = 8,
  dpi = 300
)


#---------------------------------------------------------

# Step 1: Ensure land_cover is a SpatRaster
land_cover <- rast("corine_rat_LEVEL2.tif")

# Step 2: Convert nuts2_outbreaks to sf and then to SpatVector
nuts2_outbreaks <- st_as_sf(nuts2_outbreaks)
nuts2_vect <- vect(nuts2_outbreaks)

# Step 3: Extract modal land cover class per polygon
land_cover_extract <- extract(land_cover, nuts2_vect, fun = modal, na.rm = TRUE)

# Step 4: Assign numeric land cover code to nuts2_outbreaks
nuts2_outbreaks$land_cover_code <- land_cover_extract[, 2]  # Adjust index if needed

# Step 5: Create lookup table for LEVEL2 land cover codes
land_cover_lookup <- data.frame(
  code = c(11, 12, 13, 14, 21, 22, 23, 24, 25, 31, 32, 33, 41, 42, 51, 52),
  name = c("Urban fabric", "Industrial units", "Mine sites", "Artificial vegetation",
           "Arable land", "Permanent crops", "Pastures", "Heterogeneous agriculture",
           "Forests", "Scrub/herbaceous", "Open spaces", "Wetlands", "Inland wetlands",
           "Coastal wetlands", "Inland waters", "Marine waters")
)

# Step 6: Join lookup table to get readable names
nuts2_outbreaks <- nuts2_outbreaks %>%
  left_join(land_cover_lookup, by = c("land_cover_code" = "code"))

# Optional: Rename column for clarity
nuts2_outbreaks <- nuts2_outbreaks %>%
  rename(land_cover_name = name)


# Ensure pig_density is a SpatRaster
chicken_density <- rast("Chicken_density.tif")
# If nuts2_outbreaks is not already sf, convert it
nuts2_outbreaks <- st_as_sf(nuts2_outbreaks)
# Then convert to SpatVector
nuts2_vect <- vect(nuts2_outbreaks)
# Now extract pig density
chicken_density_extract <- terra::extract(chicken_density, nuts2_vect, fun = mean, na.rm = TRUE)
# Add pig density to original sf object
nuts2_outbreaks$chicken_density <- chicken_density_extract[, 2]  # Adjust column index if needed

library(readxl)
poultry_density <- read_excel("poultry_density.xlsx")%>%
  #View(poultry_density)
  mutate(
    poultry_count = as.numeric(gsub("[^0-9.]", "", value))  # remove non-numeric characters
  ) %>%
  group_by(NAME_LATN) %>%
  summarise(poultry_density = sum(poultry_count, na.rm = TRUE))

nuts2_outbreaks_poultry <- nuts2_outbreaks

nuts2_outbreaks <- nuts2_outbreaks_poultry %>%
  left_join(poultry_density,   by = "NAME_LATN") 


#----Biosecurity measurers-------
#--All birds based on countries-- 2023/2024 --------------
X2023_reg_Biocheck_EU_BIRDS_Free_range_layers <- read_excel("2023_reg_Biocheck_EU_BIRDS.xlsx", sheet = 'Free_range_layers')
X2023_reg_Biocheck_EU_BIRDS_Free_range_broilers <- read_excel("2023_reg_Biocheck_EU_BIRDS.xlsx", sheet = 'Free_range_broilers')
X2023_reg_Biocheck_EU_BIRDS_Broilers <- read_excel("2023_reg_Biocheck_EU_BIRDS.xlsx", sheet = 'Broilers')
X2023_reg_Biocheck_EU_BIRDS_Laying_hens <- read_excel("2023_reg_Biocheck_EU_BIRDS.xlsx", sheet = 'Laying_hens')
unique(X2023_reg_Biocheck_EU_BIRDS_Free_range_layers$Country)
unique(X2023_reg_Biocheck_EU_BIRDS_Free_range_broilers$Country)
unique(X2023_reg_Biocheck_EU_BIRDS_Broilers$Country)
unique(X2023_reg_Biocheck_EU_BIRDS_Laying_hens$Country)

Biocheck_Region_Mean_Free_range_layers_2023 <- X2023_reg_Biocheck_EU_BIRDS_Free_range_layers %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J","K","L","M")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  distinct() 

Biocheck_Region_Mean_Free_range_broilers_2023 <- X2023_reg_Biocheck_EU_BIRDS_Free_range_broilers %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  distinct() 

Biocheck_Region_Mean_Broilers_2023 <- X2023_reg_Biocheck_EU_BIRDS_Broilers %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J","K")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  distinct() 

Biocheck_Region_Mean_Laying_hens_2023 <- X2023_reg_Biocheck_EU_BIRDS_Laying_hens %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J","K","L","M","N")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  distinct() 

#---------------2024-----------------
X2024_reg_Biocheck_EU_BIRDS_Free_range_layers <- read_excel("2024_reg_Biocheck_EU_BIRDS.xlsx", sheet = 'Free_range_layers')
X2024_reg_Biocheck_EU_BIRDS_Free_range_broilers <- read_excel("2024_reg_Biocheck_EU_BIRDS.xlsx", sheet = 'Free_range_broilers')
X2024_reg_Biocheck_EU_BIRDS_Broilers <- read_excel("2024_reg_Biocheck_EU_BIRDS.xlsx", sheet = 'Broilers')
X2024_reg_Biocheck_EU_BIRDS_Laying_hens <- read_excel("2024_reg_Biocheck_EU_BIRDS.xlsx", sheet = 'Laying_hens')
unique(X2024_reg_Biocheck_EU_BIRDS_Free_range_layers$Country)
unique(X2024_reg_Biocheck_EU_BIRDS_Free_range_broilers$Country)
unique(X2024_reg_Biocheck_EU_BIRDS_Broilers$Country)
unique(X2024_reg_Biocheck_EU_BIRDS_Laying_hens$Country)

Biocheck_Region_Mean_Free_range_layers_2024 <- X2024_reg_Biocheck_EU_BIRDS_Free_range_layers %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J","K","L","M")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  distinct() 

Biocheck_Region_Mean_Free_range_broilers_2024 <- X2024_reg_Biocheck_EU_BIRDS_Free_range_broilers %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  distinct() 

Biocheck_Region_Mean_Broilers_2024 <- X2024_reg_Biocheck_EU_BIRDS_Broilers %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J","K")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  distinct() 

Biocheck_Region_Mean_Laying_hens_2024 <- X2024_reg_Biocheck_EU_BIRDS_Laying_hens %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J","K","L","M","N")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  #dplyr::rename(ADMIN = Country) %>% 
  distinct() 

df1 = Biocheck_Region_Mean_Free_range_layers_2023
df2 = Biocheck_Region_Mean_Free_range_layers_2024

#Sum their value columns, treating missing values as 0.
# Full outer join using full_join
full_merged_Free_range_layers <- full_join(df1, df2, by = c("ID","NAME_LATN"), suffix = c("_df1", "_df2"))

# Replace NAs with 0 and sum the values
Free_range_layers_Mean_Value <- full_merged_Free_range_layers %>%
  mutate(
    Mean_score_Region_df1 = replace_na(Mean_score_Region_df1, 0),
    Mean_score_Region_df2 = replace_na(Mean_score_Region_df2, 0),
    Free_range_layers_Mean_value = Mean_score_Region_df1 + Mean_score_Region_df2
  ) %>%
  dplyr::select(ID, NAME_LATN, Free_range_layers_Mean_value)%>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

print(Free_range_layers_Mean_Value )
saveRDS(Free_range_layers_Mean_Value, "Free_range_layers_Mean_Value.rds")

#------Free_range_broilers-------------
df3 = Biocheck_Region_Mean_Free_range_broilers_2023
df4 = Biocheck_Region_Mean_Free_range_broilers_2024

#Sum their value columns, treating missing values as 0.
# Full outer join using full_join
full_merged_Free_range_broilers <- full_join(df3, df4, by = c("ID","NAME_LATN"), suffix = c("_df3", "_df4"))

# Replace NAs with 0 and sum the values
Free_range_broilers_Mean_Value <- full_merged_Free_range_broilers %>%
  mutate(
    Mean_score_Region_df3 = replace_na(Mean_score_Region_df3, 0),
    Mean_score_Region_df4 = replace_na(Mean_score_Region_df4, 0),
    Free_range_broilers_Mean_value = Mean_score_Region_df3 + Mean_score_Region_df4
  ) %>%
  dplyr::select(ID, NAME_LATN, Free_range_broilers_Mean_value)%>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

print(Free_range_broilers_Mean_Value )
saveRDS(Free_range_broilers_Mean_Value, "Free_range_broilers_Mean_Value.rds")

#------Broilers-------------
df5 = Biocheck_Region_Mean_Broilers_2023
df6 = Biocheck_Region_Mean_Broilers_2024

#Sum their value columns, treating missing values as 0.
# Full outer join using full_join
full_merged_Broilers <- full_join(df5, df6, by = c("ID","NAME_LATN"), suffix = c("_df5", "_df6"))

# Replace NAs with 0 and sum the values
Broilers_Mean_Value <- full_merged_Broilers %>%
  mutate(
    Mean_score_Region_df5 = replace_na(Mean_score_Region_df5, 0),
    Mean_score_Region_df6 = replace_na(Mean_score_Region_df6, 0),
    Broilers_Mean_value = Mean_score_Region_df5 + Mean_score_Region_df6
  ) %>%
  dplyr::select(ID, NAME_LATN, Broilers_Mean_value)%>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

print(Broilers_Mean_Value)
saveRDS(Broilers_Mean_Value, "Broilers_Mean_Value.rds")

#------Laying_hens-------------
df7 = Biocheck_Region_Mean_Laying_hens_2023
df8 = Biocheck_Region_Mean_Laying_hens_2024

#Sum their value columns, treating missing values as 0.
# Full outer join using full_join
full_merged_Laying_hens <- full_join(df7, df8, by = c("ID","NAME_LATN"), suffix = c("_df7", "_df8"))

# Replace NAs with 0 and sum the values
Laying_hens_Mean_Value <- full_merged_Laying_hens %>%
  mutate(
    Mean_score_Region_df3 = replace_na(Mean_score_Region_df7, 0),
    Mean_score_Region_df4 = replace_na(Mean_score_Region_df8, 0),
    Laying_hens_Mean_value = Mean_score_Region_df7 + Mean_score_Region_df8
  ) %>%
  dplyr::select(ID, NAME_LATN, Laying_hens_Mean_value)%>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

print(Laying_hens_Mean_Value)
saveRDS(Laying_hens_Mean_Value, "Laying_hens_Mean_Value.rds")


#-------Down the BioCheck for the poulty-------------
#This is for the Biocheck data which have been cleaned 
# 1. Clean helper: summarise numeric columns only
clean_metric_table <- function(df) {
  df %>%
    group_by(NAME_LATN) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
}

Free_range_layers_Mean_Value <- readRDS("Free_range_layers_Mean_Value.rds")
Free_range_broilers_Mean_Value <- readRDS("Free_range_broilers_Mean_Value.rds")
Broilers_Mean_Value <- readRDS("Broilers_Mean_Value.rds")
Laying_hens_Mean_Value <- readRDS("Laying_hens_Mean_Value.rds")


# 2. Clean each metric table
Free_range_layer_Mean_Value       <- clean_metric_table(Free_range_layers_Mean_Value)
Free_range_broilers_Mean_Value      <- clean_metric_table(Free_range_broilers_Mean_Value)
Broilers_Mean_Value       <- clean_metric_table(Broilers_Mean_Value)
Laying_hens_Mean_Value      <- clean_metric_table(Laying_hens_Mean_Value)

# 3. Join all to nuts2_drivers

nuts2_outbreaks_ai <- nuts2_outbreaks
AI_nuts2_outbreaks <- nuts2_outbreaks_ai %>%
  left_join(Free_range_layer_Mean_Value,       by = "NAME_LATN") %>%
  left_join(Free_range_broilers_Mean_Value,      by = "NAME_LATN") %>%
  left_join(Broilers_Mean_Value,       by = "NAME_LATN") %>%
  left_join(Laying_hens_Mean_Value,      by = "NAME_LATN") 


saveRDS(AI_nuts2_outbreaks, file = "AI_nuts2_outbreaks_poultry.rds")

AI_nuts2_outbreaks <- readRDS("AI_nuts2_outbreaks_poultry.rds")

#------------Plotting the biosecurity values--------------

# Poultry density
p_poultry <- ggplot(AI_nuts2_outbreaks) +
  geom_sf(aes(fill = poultry_density), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "inferno") +
  labs(
    title = "Poultry density Level by NUTS2 Region",
    fill = bquote("Head / 10 km"^2)
  ) +
  coord_sf(
    xlim = c(-10, 55),
    ylim = c(35, 70)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

print(p_poultry)


# Step 17: Save the plot
if (!dir.exists("AI_Risk_Figures")) {
  dir.create("AI_Risk_Figures")
}

ggsave(
  filename = "AI_Risk_Figures/poultry density.png",
  plot = p_poultry,
  width = 10,
  height = 8,
  dpi = 300
)



