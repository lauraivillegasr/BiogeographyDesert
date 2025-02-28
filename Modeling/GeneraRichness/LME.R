##lme genera richness##


library(nlme)
library(Hmisc)
library(arrow)

# Read input data
ecovar_data_GLS_no_categories <- read.csv2("/data/EcoVarData_transect_testing_fixed_feb.csv")

library(ggplot2)
library(dplyr)

#Transform data to fit normal distribution 
ecovar_data_GLS_no_categories$soil_thick_log<-(log((ecovar_data_GLS_no_categories$soil_thick)+1))



# Null model, excluding highly correlated variables and variables which estimation was dependent on others (1)

M_1_0=lme(composition_genera ~ Latitude+dem_ele+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa,random= ~1 | Transect, method = "REML", data=ecovar_data_GLS_no_categories)
summary(M_1_0)
plot(M_1_0)

M_1_1=lme(composition_genera ~ Latitude+dem_ele+rock_type_int+soil_thick_log+MAP_chelsa,random= ~1 | Transect, data=ecovar_data_GLS_no_categories)
summary(M_1_1)
plot(M_1_1)

M_1_2=lme(composition_genera ~ Latitude+dem_ele+rock_type_int+soil_thick_log+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_1_2)
plot(M_1_2)

M_1_3=lme(composition_genera ~ Latitude+dem_ele+rock_type_int+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_1_3)
plot(M_1_3)

M_1_4=lme(composition_genera ~ Latitude+dem_ele+soil_thick_log+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_1_4)
plot(M_1_4)

M_1_5=lme(composition_genera ~ Latitude+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_1_5)
plot(M_1_5)

M_1_6=lme(composition_genera ~ dem_ele+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_1_6)
plot(M_1_6)

###Using 4 variables at a time###

M_2_1=lme(composition_genera ~ Latitude+dem_ele+rock_type_int+soil_thick_log, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_1)
plot(M_2_1)

M_2_2=lme(composition_genera ~ Latitude+dem_ele+rock_type_int+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_2)
plot(M_2_2)

M_2_3=lme(composition_genera ~ Latitude+dem_ele+rock_type_int+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_3)
plot(M_2_3)

M_2_4=lme(composition_genera ~ Latitude+dem_ele+soil_thick_log+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_4)
plot(M_2_4)

M_2_5=lme(composition_genera ~ Latitude+dem_ele+soil_thick_log+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_5)
plot(M_2_5)

M_2_6=lme(composition_genera ~ Latitude+dem_ele+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_6)
plot(M_2_6)

M_2_7=lme(composition_genera ~ Latitude+rock_type_int+soil_thick_log+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_7)
plot(M_2_7)

M_2_8=lme(composition_genera ~ Latitude+rock_type_int+soil_thick_log+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_8)
plot(M_2_8)

M_2_9=lme(composition_genera ~ Latitude+rock_type_int+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_9)
plot(M_2_9)

M_2_10=lme(composition_genera ~ Latitude+soil_thick_log+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_10)
plot(M_2_10)


M_2_11=lme(composition_genera ~ dem_ele+rock_type_int+soil_thick_log+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_11)
plot(M_2_11)


M_2_12=lme(composition_genera ~ dem_ele+rock_type_int+soil_thick_log+range_tmp_chelsa, random= ~1 | Transect, method = "REML", data=ecovar_data_GLS_no_categories)
summary(M_2_12)
plot(M_2_12)


M_2_13=lme(composition_genera ~ dem_ele+rock_type_int+MAP_chelsa+range_tmp_chelsa,random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_13)
plot(M_2_13)

M_2_14=lme(composition_genera ~ dem_ele+soil_thick_log+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_14)
plot(M_2_14)


M_2_15=lme(composition_genera ~ rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_2_15)
plot(M_2_15)

###Using 3 variables at a time###

M_3_1=lme(composition_genera ~ Latitude+dem_ele+rock_type_int, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_1)
plot(M_3_1)

M_3_2=lme(composition_genera ~ Latitude+dem_ele+soil_thick_log, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_2)
plot(M_3_2)

M_3_3=lme(composition_genera ~ Latitude+dem_ele+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_3)
plot(M_3_3)

M_3_4=lme(composition_genera ~ Latitude+dem_ele+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_4)
plot(M_3_4)

M_3_5=lme(composition_genera ~ Latitude+rock_type_int+soil_thick_log, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_5)
plot(M_3_5)

M_3_6=lme(composition_genera ~ Latitude+rock_type_int+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_6)
plot(M_3_6)

M_3_7=lme(composition_genera ~ Latitude+rock_type_int+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_7)
plot(M_3_7)

M_3_8=lme(composition_genera ~ Latitude+soil_thick_log+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_8)
plot(M_3_8)

M_3_9=lme(composition_genera ~ Latitude+soil_thick_log+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_9)
plot(M_3_9)


M_3_10=lme(composition_genera ~ Latitude+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_10)
plot(M_3_10)

M_3_11=lme(composition_genera ~ dem_ele+rock_type_int+soil_thick_log, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_11)
plot(M_3_11)

M_3_12=lme(composition_genera ~ dem_ele+rock_type_int+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_12)
plot(M_3_12)

M_3_13=lme(composition_genera ~ dem_ele+rock_type_int+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_13)
plot(M_3_13)

M_3_14=lme(composition_genera ~ dem_ele+soil_thick_log+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_14)
plot(M_3_14)

M_3_15=lme(composition_genera ~ dem_ele+soil_thick_log+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_15)
plot(M_3_15)

M_3_16=lme(composition_genera ~ dem_ele+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_16)
plot(M_3_16)

M_3_17=lme(composition_genera ~ rock_type_int+soil_thick_log+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_17)
plot(M_3_17)

M_3_18=lme(composition_genera ~ rock_type_int+soil_thick_log+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_18)
plot(M_3_18)

M_3_19=lme(composition_genera ~ rock_type_int+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_19)
plot(M_3_19)

M_3_20=lme(composition_genera ~ soil_thick_log+MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_3_20)
plot(M_3_20)

###Using 2 variables at a time###

M_4_1=lme(composition_genera ~ Latitude+dem_ele, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_1)
plot(M_4_1)

M_4_2=lme(composition_genera ~ Latitude+rock_type_int, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_2)
plot(M_4_2)

M_4_3=lme(composition_genera ~ Latitude+soil_thick_log, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_3)
plot(M_4_3)

M_4_4=lme(composition_genera ~ Latitude+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_4)
plot(M_4_4)

M_4_5=lme(composition_genera ~ Latitude+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_5)
plot(M_4_5)

M_4_6=lme(composition_genera ~ dem_ele+rock_type_int, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_6)
plot(M_4_6)

M_4_7=lme(composition_genera ~ dem_ele+soil_thick_log, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_7)
plot(M_4_7)

M_4_8=lme(composition_genera ~ dem_ele+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_8)
plot(M_4_8)

M_4_9=lme(composition_genera ~ dem_ele+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_9)
plot(M_4_9)

M_4_10=lme(composition_genera ~ rock_type_int+soil_thick_log, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_10)
plot(M_4_10)


M_4_11=lme(composition_genera ~ rock_type_int+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_11)
plot(M_4_11)

M_4_12=lme(composition_genera ~ rock_type_int+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_12)
plot(M_4_12)


M_4_13=lme(composition_genera ~ soil_thick_log+MAP_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_13)
plot(M_4_13)


M_4_14=lme(composition_genera ~ soil_thick_log+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_14)
plot(M_4_14)

M_4_15=lme(composition_genera ~ MAP_chelsa+range_tmp_chelsa, random= ~1 | Transect, method = "REML",data=ecovar_data_GLS_no_categories)
summary(M_4_15)
plot(M_4_15)


objects <- mget(ls(pattern = "M_*"), .GlobalEnv)
compatible_objects <- Filter(function(obj) {
  inherits(try(logLik(obj), silent = TRUE), "logLik")
}, objects)


BICs<-sapply(mget(ls(pattern="M_*"),.GlobalEnv),BIC)
Dbic<-c()
Wm<-c()

for(i in 1:length(BICs)){
  tBIC=BICs[i]-min(BICs)
  Dbic<-c(Dbic,tBIC)
}

for(p in 1:length(BICs)){
  tWm=exp(-Dbic[p]/2)/sum(exp(-Dbic/2))
  Wm<-c(Wm,tWm)
}

BICs
Wm

max_Wm_value <- max(Wm)
max_Wm_value
Wm_model_name <- names(Wm)[which.max(Wm)]
Wm_model_name
