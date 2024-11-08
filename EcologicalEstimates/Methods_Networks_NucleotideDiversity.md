# Haplotype networks and nucleotide diversity of desert nematodes

## Methods

•       Genetic marker: 18S

•       Sequences with at least 20% high quality bases (HQ%)\

-   trimmed at ends (**\>1%** chance of an error per base)

•       Blasted against curated 18S-NemaBase from Gattoni et al. (2023)

•       Best hits were assigned according to lowest E-value and a high percentage of identity + manual curation of sequences based on morphology

•       Alignments: Nematode families represented by five or more 18s rRNA sequences from different samples

•       Sequences per genera were aligned (multiple alignment) using ClustalOmega and trimmed according to sequence overlap to similar or same length

•       Networks:

•       Input: Trimmed and curated alignments with reference sequence based on blast result in .nex format

 

•       DnaSP6: Construct locations as traits

•       Open Data File -- Load nexus

•       Format - Haploid

•       Assign genetic code -- Nuclear Universal

•       Define Sequence Sets -- Add locations per group

•       Calculate nucletotide diversity and theta -- for some sort by location

•       Generate Haplotype Data File -- invariable sites included = Input for PopArt

 

If nexus file cannot be read, change missing to N

 

EPT = Eagle_Point

ARO = Aroma

TDT = Totoral_Dunes

ALT = Altiplano

PAP = Paposo

LAG = Laguna_Grande

SLH = Salar_de_Huasco

From Blast = Reference

 

•       Inspect haplotypes in text editor

•       Copy #Hap Freq. Sequences

•       Create Presence/Absence Matrix of each hapoytpes for the locations

                        LocA    LocB

                        Hap_1 0          10

                        Hap_2 1          0

•       Save as trait file

 

•       PopArt: Read in alignment file and trait-file

•       Create TCS network

•       Export summary

•       Save network as image

•       Test correlation of Pi against geographic location with Mantel's test

•       library(tidyverse) \# A collection of R packages for data manipulation and visualization.

•       library(geodist) \# For calculating geographical distances.

•       library(sf) \# For handling spatial data.

•       library(vegan) \# For ecological data analysis, including diversity analysis.

•       library(pheatmap) \# For creating heatmaps.

•       library(fuzzySim) \# For creating presence-absence matrices.

•       Obtain centroids of the different sampling spots

•       Calculate distances between the centroid of each of the transects

•       Calculate the distance matrix - dividing so it is in km

•       Remove reference

•       Remove LagunaGrande & SalarDeHuasco \--\> Salars is the mean of both

•       Filter for ASV

•       Remove rows with NA values in the Pi column

•       Create a pairwise distance matrix from Pi values using Euclidean distance

•       Perform the Mantel test

•       Make heatmap and scatterplot
