# Compare nucleotide diversity to geographic location#
# Test correlation of Pi against geographic location with Mantel's test#

library(tidyverse) 
library(geodist) 
library(sf) 
library(vegan) 
library(pheatmap)
library(fuzzySim) 


# CALCULATE CENTROIDS

# Read input data
Location_forcentroids <- read.csv2("/data/Location_forcentroids.csv")

# Convert the data frame to an sf object with spatial coordinates
Location_forcentroids <- st_as_sf(Location_forcentroids, coords = c('longitude', 'latitude'), crs = 4326)

# Obtain centroids of the different sampling spots
centroids_sf <- Location_forcentroids %>%
  group_by(Location) %>% 
  summarize(geometry = st_union(geometry)) %>% 
  st_centroid()

# Check the centroids_sf object
print(centroids_sf)

# Extract coordinates and convert to a data frame
coordinates_df <- centroids_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(longitude = X, latitude = Y) %>%
  bind_cols(centroids_sf %>% st_drop_geometry() %>% select(Location))

print(coordinates_df)

# Calculate distances between the centroid of each of the transects
# Calculate the distance matrix - dividing so it is in km
distance_matrix <- geodist(coordinates_df[, c("longitude", "latitude")], measure = "geodesic") / 1000

# Convert the distance matrix to a data frame for better readability
distance_matrix_df <- as.data.frame(distance_matrix)
rownames(distance_matrix_df) <- coordinates_df$Location
colnames(distance_matrix_df) <- coordinates_df$Location

# Print the distance matrix
print(distance_matrix_df)


##################################################################
# Nucleotide diversity compared to geographic location
# 1. Plectidae
# 2. Alaimidae - no test, only 1 haplotype
# 3. Acrobeloides
# 4. Dorylamida
# 5. Acrobeles
# 6. Panagrolaimus
# 7. Aphelenchoididae - no test, only 1 value in distance matrix
##################################################################
#------------------------------
# 1. PLECTIDAE
#------------------------------
# Read nucleotide diversity data
nucleotide_diversity <- read.csv("/data/nucleotide_diversity_salarstogether.csv", sep = ";")

# Print nucleotide diversity data
print(nucleotide_diversity)

# Remove reference
# Remove LagunaGrande & SalarDeHuasco --> Salars is the mean of both
nucleotide_diversity <- nucleotide_diversity %>%
  filter(Location != "Reference") 

# Filter for Acrobeles
nucleotide_diversity <- nucleotide_diversity %>%
  mutate(Location = Location,
         Pi = Plectidae)

# Remove rows with NA values in the Pi column
nucleotide_diversity <- nucleotide_diversity %>% drop_na(Pi)

# Create a pairwise distance matrix from Pi values using Euclidean distance
pi_dist_matrix <- dist(nucleotide_diversity$Pi, method = "euclidean")

print(pi_dist_matrix)

# Convert the Pi distance matrix to a data frame
pi_distance_matrix_df <- as.matrix(pi_dist_matrix)
rownames(pi_distance_matrix_df) <- nucleotide_diversity$Location
colnames(pi_distance_matrix_df) <- nucleotide_diversity$Location

# Filter geo_distance_matrix to only include locations in pi_distance_matrix
common_locations <- intersect(rownames(pi_distance_matrix_df), rownames(distance_matrix_df))
filtered_geo_distance_matrix <- distance_matrix_df[common_locations, common_locations]

# Convert the filtered geo_distance_matrix to a dist object
geo_distance_matrix <- as.dist(filtered_geo_distance_matrix)

# Convert the Pi distance matrix to a dist object
pi_distance_matrix <- as.dist(pi_distance_matrix_df)

# Perform the Mantel test
# Go with Spearman
mantel_test_result_plectidae <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "pearson", permutations = 9999)
mantel_test_result_plectidae <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "spearman", permutations = 9999)

# Print the Mantel test result
print(mantel_test_result_plectidae)

mantel_test_result_plectidae[["signif"]] < 0.05




#------------------------------
# 2. ALAIMIDAE
# Alaimidae just one haplotype
# No Mantel test
#------------------------------




#------------------------------
# 3. ACROBELOIDES
#------------------------------
# Read nucleotide diversity data
nucleotide_diversity <- read.csv("/data/nucleotide_diversity_salarstogether.csv", sep = ";")

# Print nucleotide diversity data
print(nucleotide_diversity)

# Remove reference
# Remove LagunaGrande & SalarDeHuasco --> Salars is the mean of both
nucleotide_diversity <- nucleotide_diversity %>%
  filter(Location != "Reference") 
# Filter for Acrobeles
nucleotide_diversity <- nucleotide_diversity %>%
  mutate(Location = Location,
         Pi = Acrobeloides)

# Remove rows with NA values in the Pi column
nucleotide_diversity <- nucleotide_diversity %>% drop_na(Pi)

# Create a pairwise distance matrix from Pi values using Euclidean distance
pi_dist_matrix <- dist(nucleotide_diversity$Pi, method = "euclidean")

print(pi_dist_matrix)

# Convert the Pi distance matrix to a data frame
pi_distance_matrix_df <- as.matrix(pi_dist_matrix)
rownames(pi_distance_matrix_df) <- nucleotide_diversity$Location
colnames(pi_distance_matrix_df) <- nucleotide_diversity$Location

# Filter geo_distance_matrix to only include locations in pi_distance_matrix
common_locations <- intersect(rownames(pi_distance_matrix_df), rownames(distance_matrix_df))
filtered_geo_distance_matrix <- distance_matrix_df[common_locations, common_locations]

# Convert the filtered geo_distance_matrix to a dist object
geo_distance_matrix <- as.dist(filtered_geo_distance_matrix)

# Convert the Pi distance matrix to a dist object
pi_distance_matrix <- as.dist(pi_distance_matrix_df)

# Perform the Mantel test
mantel_test_result_acrobeloides <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "pearson", permutations = 9999)
mantel_test_result_acrobeloides <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "spearman", permutations = 9999)

# Print the Mantel test result
print(mantel_test_result_acrobeloides)

mantel_test_result_acrobeloides[["signif"]] < 0.05




#------------------------------
# 4. DORYLAMIDA
#------------------------------
# Read nucleotide diversity data
nucleotide_diversity <- read.csv("/data/nucleotide_diversity_salarstogether.csv", sep = ";")

# Print nucleotide diversity data
print(nucleotide_diversity)

# Remove reference
# Remove LagunaGrande & SalarDeHuasco --> Salars is the mean of both
nucleotide_diversity <- nucleotide_diversity %>%
  filter(Location != "Reference") 

# Filter for Acrobeles
nucleotide_diversity <- nucleotide_diversity %>%
  mutate(Location = Location,
         Pi = Dorylamida)

# Remove rows with NA values in the Pi column
nucleotide_diversity <- nucleotide_diversity %>% drop_na(Pi)

# Create a pairwise distance matrix from Pi values using Euclidean distance
pi_dist_matrix <- dist(nucleotide_diversity$Pi, method = "euclidean")

print(pi_dist_matrix)

# Convert the Pi distance matrix to a data frame
pi_distance_matrix_df <- as.matrix(pi_dist_matrix)
rownames(pi_distance_matrix_df) <- nucleotide_diversity$Location
colnames(pi_distance_matrix_df) <- nucleotide_diversity$Location

# Filter geo_distance_matrix to only include locations in pi_distance_matrix
common_locations <- intersect(rownames(pi_distance_matrix_df), rownames(distance_matrix_df))
filtered_geo_distance_matrix <- distance_matrix_df[common_locations, common_locations]

# Convert the filtered geo_distance_matrix to a dist object
geo_distance_matrix <- as.dist(filtered_geo_distance_matrix)

# Convert the Pi distance matrix to a dist object
pi_distance_matrix <- as.dist(pi_distance_matrix_df)

# Perform the Mantel test
# Go with Spearman
mantel_test_result_dorylaimida <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "pearson", permutations = 9999)
mantel_test_result_dorylaimida <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "spearman", permutations = 9999)

# Print the Mantel test result
print(mantel_test_result_dorylaimida)

mantel_test_result_dorylaimida[["signif"]] < 0.05




#------------------------------
#  5. ACROBELES
#------------------------------
# Read nucleotide diversity data
nucleotide_diversity <- read.csv("/data/nucleotide_diversity_salarstogether.csv", sep = ";")

# Print nucleotide diversity data
print(nucleotide_diversity)

# Remove reference
# Remove LagunaGrande & SalarDeHuasco --> Salars is the mean of both
nucleotide_diversity <- nucleotide_diversity %>%
  filter(Location != "Reference") 

# Filter for Acrobeles
nucleotide_diversity <- nucleotide_diversity %>%
  mutate(Location = Location,
         Pi = Acrobeles)

# Remove rows with NA values in the Pi column
nucleotide_diversity <- nucleotide_diversity %>% drop_na(Pi)

# Create a pairwise distance matrix from Pi values using Euclidean distance
pi_dist_matrix <- dist(nucleotide_diversity$Pi, method = "euclidean")

print(pi_dist_matrix)

# Convert the Pi distance matrix to a data frame
pi_distance_matrix_df <- as.matrix(pi_dist_matrix)
rownames(pi_distance_matrix_df) <- nucleotide_diversity$Location
colnames(pi_distance_matrix_df) <- nucleotide_diversity$Location

# Filter geo_distance_matrix to only include locations in pi_distance_matrix
common_locations <- intersect(rownames(pi_distance_matrix_df), rownames(distance_matrix_df))
filtered_geo_distance_matrix <- distance_matrix_df[common_locations, common_locations]

# Convert the filtered geo_distance_matrix to a dist object
geo_distance_matrix <- as.dist(filtered_geo_distance_matrix)

# Convert the Pi distance matrix to a dist object
pi_distance_matrix <- as.dist(pi_distance_matrix_df)

# Perform the Mantel test
mantel_test_result_acrobeles <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "pearson", permutations = 9999)
mantel_test_result_acrobeles <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "spearman", permutations = 9999)

# Print the Mantel test result
print(mantel_test_result_acrobeles)

mantel_test_result_acrobeles[["signif"]] < 0.05


#------------------------------
# 6. PANAGROLAIMUS
#------------------------------
# Read nucleotide diversity data
nucleotide_diversity <- read.csv("/data/nucleotide_diversity_salarstogether.csv", sep = ";")

# Print nucleotide diversity data
print(nucleotide_diversity)

# Remove reference
# Remove LagunaGrande & SalarDeHuasco --> Salars is the mean of both
nucleotide_diversity <- nucleotide_diversity %>%
  filter(Location != "Reference") 

# Filter for Acrobeles
nucleotide_diversity <- nucleotide_diversity %>%
  mutate(Location = Location,
         Pi = Panagrolaimidae)

# Remove rows with NA values in the Pi column
nucleotide_diversity <- nucleotide_diversity %>% drop_na(Pi)

# Create a pairwise distance matrix from Pi values using Euclidean distance
pi_dist_matrix <- dist(nucleotide_diversity$Pi, method = "euclidean")

print(pi_dist_matrix)

# Convert the Pi distance matrix to a data frame
pi_distance_matrix_df <- as.matrix(pi_dist_matrix)
rownames(pi_distance_matrix_df) <- nucleotide_diversity$Location
colnames(pi_distance_matrix_df) <- nucleotide_diversity$Location


# Filter geo_distance_matrix to only include locations in pi_distance_matrix
common_locations <- intersect(rownames(pi_distance_matrix_df), rownames(distance_matrix_df))
filtered_geo_distance_matrix <- distance_matrix_df[common_locations, common_locations]

# Convert the filtered geo_distance_matrix to a dist object
geo_distance_matrix <- as.dist(filtered_geo_distance_matrix)

# Convert the Pi distance matrix to a dist object
pi_distance_matrix <- as.dist(pi_distance_matrix_df)

# Perform the Mantel test
mantel_test_result_panagrolaimus <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "pearson", permutations = 9999)
mantel_test_result_panagrolaimus <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "spearman", permutations = 9999)
?mantel()

# Print the Mantel test result
print(mantel_test_result_panagrolaimus)

mantel_test_result_panagrolaimus[["signif"]] < 0.05


#------------------------------
# 7. APHELENCHOIDIDAE 
# --> only 1 distance, no test possible
#------------------------------
# Read nucleotide diversity data
nucleotide_diversity <- read.csv("/data/nucleotide_diversity_salarstogether.csv", sep = ";")

# Print nucleotide diversity data
print(nucleotide_diversity)

# Remove reference
# Remove LagunaGrande & SalarDeHuasco --> Salars is the mean of both
nucleotide_diversity <- nucleotide_diversity %>%
  filter(Location != "Reference") 

# Filter for Acrobeles
nucleotide_diversity <- nucleotide_diversity %>%
  mutate(Location = Location,
         Pi = Aphelenchoididae)

# Remove rows with NA values in the Pi column
nucleotide_diversity <- nucleotide_diversity %>% drop_na(Pi)

# Create a pairwise distance matrix from Pi values using Euclidean distance
pi_dist_matrix <- dist(nucleotide_diversity$Pi, method = "euclidean")

print(pi_dist_matrix)

# Convert the Pi distance matrix to a data frame
pi_distance_matrix_df <- as.matrix(pi_dist_matrix)
rownames(pi_distance_matrix_df) <- nucleotide_diversity$Location
colnames(pi_distance_matrix_df) <- nucleotide_diversity$Location

# Filter geo_distance_matrix to only include locations in pi_distance_matrix
common_locations <- intersect(rownames(pi_distance_matrix_df), rownames(distance_matrix_df))
filtered_geo_distance_matrix <- distance_matrix_df[common_locations, common_locations]

# Convert the filtered geo_distance_matrix to a dist object
geo_distance_matrix <- as.dist(filtered_geo_distance_matrix)

# Convert the Pi distance matrix to a dist object
pi_distance_matrix <- as.dist(pi_distance_matrix_df)

# Perform the Mantel test
# Go with Spearman
mantel_test_result_aphelenchoididae <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "pearson", permutations = 9999)
mantel_test_result_aphelenchoididae <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "spearman", permutations = 9999)

# Print the Mantel test result
print(mantel_test_result_aphelenchoididae)

mantel_test_result_aphelenchoididae[["signif"]] < 0.05
