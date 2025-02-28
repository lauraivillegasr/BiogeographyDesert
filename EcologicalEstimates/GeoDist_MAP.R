#Estimate Jaccard disimilarity and perform Mantel test to check IBD patterns#

library(tidyverse)
library(geodist)
library(sf)
library(vegan)
library(pheatmap)
library(fuzzySim)

#Load input data

Location_forcentroids <- read.csv2("/data/Location_forcentroids.csv")
ecovar_data_GLS_no_categories <- read.csv2("/data/EcoVarData_transect_testing_fixed_feb.csv")
# Obtain centroids of the different sampling spots

Location_forcentroids <- st_as_sf(Location_forcentroids, coords = c('longitude', 'latitude'))

centroids_sf <- Location_forcentroids %>%
  group_by(Location) %>% 
  summarize(geometry = st_union(geometry)) %>% 
  st_centroid

centroids_sf

# Extract coordinates and convert to a data frame
coordinates_df <- centroids_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(longitude = X, latitude = Y) %>%
  bind_cols(centroids_sf %>% st_drop_geometry() %>% select(Location))
coordinates_df

#calculate distances between the centroid of each of the transects
# Calculate the distance matrix - dividing so it is in km
distance_matrix <- geodist(coordinates_df[, c("longitude", "latitude")], measure = "geodesic")/1000

# Convert the distance matrix to a data frame for better readability
distance_matrix_df <- as.data.frame(distance_matrix)
rownames(distance_matrix_df) <- centroids_sf$Location
colnames(distance_matrix_df) <- centroids_sf$Location

# Print the distance matrix
print(distance_matrix_df)

# Obtain average MAP for each transect
# Compute mean MAP_chelsa per Transect (location)
average_precip <- ecovar_data_GLS_no_categories %>%
  group_by(Transect) %>%
  summarise(Mean_Precip = mean(MAP_chelsa, na.rm = TRUE))

# Print results
print(average_precip)

#Match names of both files (codes vs whole transect name)

# Create a lookup table to map short codes to full names
transect_mapping <- data.frame(
  ShortCode = c("ALT", "ARO", "EPT", "SLH", "LAG", "TDT", "PAP"),
  FullName  = c("Altiplano", "Aroma", "EaglePoint", "Salars", "Salars", "TotoralDunes", "Paposo")
)

# Print the mapping table
print(transect_mapping)

#Now get average precipitation with correct names
average_precip <- merge(average_precip, transect_mapping, by.x = "Transect", by.y = "ShortCode")

# Keep only the Full Name and Mean Precipitation
average_precip <- average_precip %>%
  select(FullName, Mean_Precip) %>%
  rename(Location = FullName)

# Print results
print(average_precip)
# Aggregate to merge duplicate locations (Salars)
average_precip <- average_precip %>%
  group_by(Location) %>%
  summarise(Mean_Precip = mean(Mean_Precip, na.rm = TRUE))

# Print to check the result
print(average_precip)

# Convert to a matrix of pairwise differences
precip_matrix <- as.matrix(dist(average_precip$Mean_Precip))

# Set row and column names to match distance matrix
rownames(precip_matrix) <- average_precip$Location
colnames(precip_matrix) <- average_precip$Location

# Print the precipitation difference matrix
print(precip_matrix)

# Ensure row and column names match between matrices
rownames(distance_matrix_df) <- colnames(distance_matrix_df) <- rownames(precip_matrix)
colnames(precip_matrix) <- rownames(precip_matrix)

# Get common sites
common_sites <- intersect(rownames(distance_matrix_df), rownames(precip_matrix))

# Subset both matrices
geo_matrix <- distance_matrix_df[common_sites, common_sites]
precip_matrix <- precip_matrix[common_sites, common_sites]

# Convert to distance objects
geo_dist <- as.dist(geo_matrix)
precip_dist <- as.dist(precip_matrix)

# Run Mantel test
mantel_result <- mantel(geo_dist, precip_dist, method = "spearman")

# Print results
print(mantel_result)
