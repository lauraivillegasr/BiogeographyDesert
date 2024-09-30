#Variable correlation plot#

library(GGally)
library(factoextra)
library(ggpubr)
library(corrplot)
library(Hmisc)
library(tidyverse)

# Read input data

ecovar_data_GLS_no_categories <- read.csv2("/data/EcoVarData_generacurated_nocategories.csv")


#data rearrangement so I can have a nice plot of variable correlation where geomorphological is grouped and climatic is grouped
data_cleaned <- ecovar_data_GLS_no_categories %>% select(-contains("composition"))
cols_to_move <- names(data_cleaned)[grepl("tmp|pre|MAT|MAP", names(data_cleaned))]
data_rearranged <- data_cleaned %>% select(-one_of(cols_to_move), everything(), one_of(cols_to_move))

# Plot the correlation between variables
corvar<-cor(data_rearranged)
corrplot(corvar)
corr <- data.frame(lapply(data_rearranged, as.integer))
p_values <- corr$P 

# Define the output directory and file name
out_dir <- "../results/"
correlation_plot_file <- "correlation_plot.png"
correlation_plot_path <- paste0(out_dir, correlation_plot_file)  # Fixed missing parenthesis

# Open a PNG device
png(correlation_plot_path, width = 800, height = 600)

# Plot the correlation matrix 
corrplot(corvar, method = "circle", 
         p.mat = p_values, sig.level = 0.05, insig = "blank", 
         tl.cex = 0.8, tl.col = "black", cl.pos = "r", 
         addCoefasPercent = FALSE, number.cex = 0.7, 
         type="lower",
         addCoef.col = NULL,
         is.corr = FALSE,
         diag = FALSE,
         order="original")

# Close the device
dev.off()
