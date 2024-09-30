Estimate Jaccard disimilarity and perform Mantel test to check IBD patterns#

library(tidyverse)
library(geodist)
library(sf)
library(vegan)
library(pheatmap)
library(fuzzySim)

#Load input data

Location_forcentroids <- read.csv2("/data/Location_forcentroids.csv")

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

##Estmiate Jaccard

Jaccard_list_corrected <- read.csv("/data/Jaccard_list_corrected.csv", sep=";")
