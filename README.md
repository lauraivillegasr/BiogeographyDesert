# Hierarchical Patterns of Soil Biodiversity in the Atacama Desert: Insights Across Biological Scales

Here you will find all the code implemented in the analysis of the manuscriot titled "Hierarchical Patterns of Soil Biodiversity in the Atacama Desert: Insights Across Biological Scales". The analysis is divided in two blocks: Ecological estimates and Modeling. 

1. Inside the Ecological Estimates folder you will find the scripts for obtaining jaccard dissimalirity index and perform a mantel test, as well as information on how to generate the haplotype networks to further estimate eucledian distance and nucleotide diversity.

* To run ecological estimates and the different models in R, the following packages are quired on R:

  <center>
  
| Package       | Version  |
|---------------|----------|
| FactoMineR    | 2.11     |
| GGally        | 2.2.1    |
| Hmisc         | 5.1-3    |
| Metrics       | 0.1.4    |
| arrow         | 16.1.0   |
| caret         | 6.0-94   |
| factoextra    | 1.0.7    |
| fuzzySim      | 4.10.7   |
| geodist       | 0.1.0    |
| ggplot2       | 3.5.1    |
| ggpubr        | 0.6.0    |
| nlme          | 3.1-165  |
| pheatmap      | 1.0.12   |
| plotrix       | 3.8-4    |
| randomForest  | 4.7-1.1  |
| ranger        | 0.16.0   |
| rpart         | 4.1.23   |
| rsample       | 1.2.1    |
| rstatix       | 0.7.2    |
| sf            | 1.0-16   |
| tidyverse     | 2.0.0    |
| vegan         | 2.6-6.1  |

</center>

2. In the modeling folder, you will find scripts for modeling using linear models and random forest for both genera richness and reproductive mode estimation. 

* To gather ecological factors for the set of coordinates of each of the sampling spots using python 3, the following are reuiqred:

 <center>

| Package             | Version   |
|---------------------|-----------|
| pandas              | 2.1.1     |
| xarray              | 2023.9.0  |
| matplotlib          | 3.8.0     |
| numpy               | 1.26.0    |
| richdem             | 0.3.4     |
| rioxarray           | 0.18.1    |
| cmcrameri           | 1.7.1     |
| geopandas           | 0.14.0    |
| shapely             | 2.0.1     |
| matplotlib.colors   | 3.8.0     |
| mpl_toolkits.axes_grid1 | 3.8.0 |

</center>

**The necessary input files for for Ecological estimates and modelling can be found in Zenodo ([10.5281/zenodo.13880073](https://zenodo.org/records/14131623)). File names on the zenodo folder correspond to the file names in the different scirpts.**




