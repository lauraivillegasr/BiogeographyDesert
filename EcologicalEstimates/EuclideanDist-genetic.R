#------------------------------
# PANAGROLAIMUS
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
mantel_test_result <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "pearson", permutations = 9999)
mantel_test_result <- mantel(xdis = pi_distance_matrix, ydis = geo_distance_matrix, method = "spearman", permutations = 9999)
?mantel()

# Print the Mantel test result
print(mantel_test_result)

mantel_test_result[["signif"]] < 0.05

out_dir <- "../results/"
Panagrolaimusplot <- "Panagrolaimusplot.png" 
Panagrolaimusplot <- paste0(out_dir, Panagrolaimusplot)
png(Panagrolaimusplot, width = 800, height = 600)

# Define theme
mytheme <- theme(axis.text.x = element_text(size = 14, color = "black"),
                 axis.text.y = element_text(size = 14, color = "black"),
                 axis.title.y = element_text(size = 14,face = "bold", color = "black"),
                 axis.title.x = element_text(size = 14,face = "bold", color = "black"),
                 title = element_text(size = 12, color = "black"),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line("lightgrey"),
                 panel.grid.minor = element_line("white"),
                 panel.border = element_rect("black", fill = NA),
                 plot.background = element_rect(fill="white"),
                 legend.background = element_rect(fill="white"),
                 legend.position="bottom")

# Plot with ggplot
Panagrolaimusplot<-ggplot() + geom_point(aes(y = pi_distance_matrix, x = geo_distance_matrix), size = 3) +
  labs(x = "Geographic distance (km)", y = "Euclidean distance of \u03c0") +
  mytheme
print(Panagrolaimusplot)
dev.off()
