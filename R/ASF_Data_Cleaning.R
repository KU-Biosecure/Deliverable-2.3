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

#ASF-Pig data
ASF_Swine_Pigs <- read_excel("ASF_Swine_Pigs.xlsx")

# Convert ASF outbreak data to sf object
ASF_Swine_Pigs_sf <- ASF_Swine_Pigs %>%
  filter(!is.na(Longitude), !is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(nuts2))

# Spatial join: assign each outbreak to a NUTS2 region
ASF_with_NUTS2 <- st_join(ASF_Swine_Pigs_sf, nuts2, join = st_intersects)

# Aggregate outbreak counts by NUTS2 region
outbreak_counts <- ASF_with_NUTS2 %>%
  st_drop_geometry() %>%
  group_by(NUTS_ID) %>%
  #group_by(NAME_LATN.x) %>%
  summarise(outbreak_count = n())


# Now join outbreak counts to nuts2 geometries
nuts2_outbreaks <- nuts2 %>%
  left_join(outbreak_counts, by = "NUTS_ID") %>%
  mutate(outbreak_count = replace_na(outbreak_count, 0))

# Create a new column to classify regions
nuts2_outbreaks <- nuts2_outbreaks %>%
  mutate(outbreak_status = ifelse(outbreak_count > 0, "Outbreak", "No Outbreak"))

# Plot the map
p_cases <- ggplot(nuts2_outbreaks) +
  geom_sf(aes(fill = outbreak_status), color = "white", size = 0.1) +
  scale_fill_manual(values = c("Outbreak" = "red", "No Outbreak" = "lightgray")) +
  labs(title = "ASF Risk Map (Regions with Outbreak OR Control & Risk Score)",
       fill = "Status") +
  coord_sf(xlim = c(-10, 55), ylim = c(35, 70)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())


# Print to viewer
print(p_cases)

# Create folder if it doesn't exist
if (!dir.exists("ASF_Risk_Figures")) {
  dir.create("ASF_Risk_Figures")
}

# Save the plot to the folder
ggsave(
  filename = "ASF_Risk_Figures/asf_both_map.png",
  plot = p_cases,
  width = 10,
  height = 8,
  dpi = 300
)


#------------------------------------------------------

#Read in the bio climate and land cover data
pig_density <-rast("Pig_density.tif")
#wild_boar <-rast("wildboar_density.tif")
land_cover <- rast("corine_rat_LEVEL2.tif")

#extract raster values to point locations where the diseases were found
#using the terra package
#WOAH
#to extract the livestock for wild boar density

library(terra)
library(sf)

# Ensure pig_density is a SpatRaster
pig_density <- rast("Pig_density.tif")
# If nuts2_outbreaks is not already sf, convert it
nuts2_outbreaks <- st_as_sf(nuts2_outbreaks)
# Then convert to SpatVector
nuts2_vect <- vect(nuts2_outbreaks)
# Now extract pig density
pig_density_extract <- terra::extract(pig_density, nuts2_vect, fun = mean, na.rm = TRUE)
# Add pig density to original sf object
nuts2_outbreaks$pig_density <- pig_density_extract[, 2]  # Adjust column index if needed


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


#----Biosecurity measurers-------
#--Indoor Pigs based on countries-- 2023/2024 --------------
X2023_reg_Biocheck_EU <- read_excel("2023_reg_Biocheck_EU.xlsx")
colnames(X2023_reg_Biocheck_EU)
unique(X2023_reg_Biocheck_EU$Country)

Biocheck_Region_InPigMean_2023 <- X2023_reg_Biocheck_EU %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J","K","L")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  #dplyr::rename(ADMIN = Country) %>% 
  distinct() 

X2024_reg_Biocheck_EU <- read_excel("2024_reg_Biocheck_EU.xlsx", sheet= "Pigs_indoor")
colnames(X2024_reg_Biocheck_EU)
unique(X2024_reg_Biocheck_EU$Country)

Biocheck_Region_InPigMean_2024 <- X2024_reg_Biocheck_EU %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  #dplyr::group_by(Country,question ) %>% #dplyr::group_by(Country) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J","K","L")) %>% # Select multiple multiple questions
  dplyr::mutate(Mean_score_Region= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Mean_score_Region ) %>% 
  #dplyr::rename(ADMIN = Country) %>% 
  distinct() 

df1 = Biocheck_Region_InPigMean_2023
df2 = Biocheck_Region_InPigMean_2024

#Sum their value columns, treating missing values as 0.
# Full outer join using full_join
full_merged <- full_join(df1, df2, by = c("ID","NAME_LATN"), suffix = c("_df1", "_df2"))

# Replace NAs with 0 and sum the values
Indoor_Mean_Value <- full_merged %>%
  mutate(
    Mean_score_Region_df1 = replace_na(Mean_score_Region_df1, 0),
    Mean_score_Region_df2 = replace_na(Mean_score_Region_df2, 0),
    Indoor_Mean_value = Mean_score_Region_df1 + Mean_score_Region_df2
  ) %>%
  dplyr::select(ID, NAME_LATN, Indoor_Mean_value)%>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

#print(Indoor_Mean_Value )

saveRDS(Indoor_Mean_Value, file = "Indoor_Mean_Value.rds")


##-------Outdoor-------------------------
V2024_reg_Biocheck_EU <- read_excel("2024_reg_Biocheck_EU.xlsx", sheet= "Pigs_outdoor")
colnames(V2024_reg_Biocheck_EU)
unique(V2024_reg_Biocheck_EU$Country)

Biocheck_Region_OutPigMean_2024 <- V2024_reg_Biocheck_EU %>% 
  mutate(ID = sub("\\..*", "", question)) %>% 
  dplyr::group_by(NAME_LATN,ID) %>% #dplyr::group_by(Regional Level) %>% 
  dplyr::filter(ID %in% c("A","B","C","D","E", "F","G","H","I","J","K","L")) %>% # Select multiple multiple questions
  dplyr::mutate(Outdoor_Mean_value= mean(`Mean value`, na.rm = T)) %>% 
  dplyr::select(ID,NAME_LATN,Outdoor_Mean_value ) %>% 
  #dplyr::rename(ADMIN = Country) %>% 
  distinct() 

Outdoor_Mean_Value <- Biocheck_Region_OutPigMean_2024 %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

saveRDS(Outdoor_Mean_Value, file = "Outdoor_Mean_Value.rds")


#This is for the Biocheck data which have been cleaned 
# 1. Clean helper: summarise numeric columns only
clean_metric_table <- function(df) {
  df %>%
    group_by(NAME_LATN) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
}

Indoor_Mean_Value <- readRDS("Indoor_Mean_Value.rds")
Outdoor_Mean_Value <- readRDS("Outdoor_Mean_Value.rds")

# 2. Clean each metric table
Indoor_Mean_Value       <- clean_metric_table(Indoor_Mean_Value)
Outdoor_Mean_Value      <- clean_metric_table(Outdoor_Mean_Value)


# 3. Join all to nuts2_drivers
ASF_nuts2_outbreaks <- nuts2_outbreaks %>%
  left_join(Indoor_Mean_Value,       by = "NAME_LATN") %>%
  left_join(Outdoor_Mean_Value,      by = "NAME_LATN") 

saveRDS(ASF_nuts2_outbreaks, file = "ASFf_nuts2_outbreaks.rds")

ASF_nuts2_outbreaks <- readRDS("ASFf_nuts2_outbreaks.rds")

#------------Plotting the biosecurity values--------------

# Pig density
p_pig <- ggplot(ASF_nuts2_outbreaks) +
  geom_sf(aes(fill = pig_density), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "inferno") +
  labs(
    title = "Pig density Level by NUTS2 Region",
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

print(p_pig)


# Step 17: Save the plot
if (!dir.exists("ASF_Risk_Figures")) {
  dir.create("ASF_Risk_Figures")
}

ggsave(
  filename = "ASF_Risk_Figures/pig density.png",
  plot = p_pig,
  width = 10,
  height = 8,
  dpi = 300
)


# Outdoor mean value
p_Outdoor <- ggplot(ASF_nuts2_outbreaks) +
  geom_sf(aes(fill = Outdoor_Mean_value), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "inferno", limits = c(0, 100)) + #"plasma"
  labs(title = "Outdoor mean value for pig by NUTS2 Region",
       caption = "Grey = all NUTS2 regions without outdoor mean value.",
       fill = "Outdoor Mean Value") +
  coord_sf(
    xlim = c(-10, 55),   # longitude: west to east
    ylim = c(35, 70)     # latitude: south to north
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
print(p_Outdoor)


# Step 17: Save the plot
if (!dir.exists("ASF_Risk_Figures")) {
  dir.create("ASF_Risk_Figures")
}

ggsave(
  filename = "ASF_Risk_Figures/Outdoor mean value.png",
  plot = p_Outdoor,
  width = 10,
  height = 8,
  dpi = 300
)


# Indoor mean value
p_Indoor <- ggplot(ASF_nuts2_outbreaks) +
  geom_sf(aes(fill = Indoor_Mean_value), color = "white", size = 0.1) +
  scale_fill_viridis_c(option = "inferno", limits = c(0, 100)) + #"plasma"
  labs(title = "Indoor mean value for pig by NUTS2 Region", 
       caption = "Grey = all NUTS2 regions without indoor mean value.",
       fill = "Indoor Mean Value") +
  coord_sf(
    xlim = c(-10, 55),   # longitude: west to east
    ylim = c(35, 70)     # latitude: south to north
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
print(p_Indoor)

# Step 17: Save the plot
if (!dir.exists("ASF_Risk_Figures")) {
  dir.create("ASF_Risk_Figures")
}

ggsave(
  filename = "ASF_Risk_Figures/indoor mean value.png",
  plot = p_Indoor,
  width = 10,
  height = 8,
  dpi = 300
)

