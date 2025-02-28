##GLS genera richness##


library(nlme)
library(Hmisc)
library(arrow)

# Read input data

ecovar_data_GLS_no_categories <- read.csv2("/data/EcoVarData_transect_testing_fixed_feb.csv")

#Transform data to fit normal distribution 
ecovar_data_GLS_no_categories$soil_thick_log<-(log((ecovar_data_GLS_no_categories$soil_thick)+1))


#total number of models to be tested: 57 including null model#

# Null model, excluding highly correlated variables and variables which estimation was dependent on others (1)

vP0_1 <- varPower(form=~Latitude+dem_ele+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa)
M_1_0=gls(composition_genera ~ Latitude+dem_ele+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
             weights=vP0_1)
summary(M_1_0)
plot(M_1_0)

###Using 5 variables at a time###

vP1_1 <- varPower(form=~Latitude+dem_ele+rock_type_int+soil_thick_log+MAP_chelsa)
M_1_1=gls(composition_genera ~ Latitude+dem_ele+rock_type_int+soil_thick_log+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP1_1)
summary(M_1_1)
plot(M_1_1)

vP1_2 <- varPower(form=~Latitude+dem_ele+rock_type_int+soil_thick_log+range_tmp_chelsa)
M_1_2=gls(composition_genera ~ Latitude+dem_ele+rock_type_int+soil_thick_log+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP1_2)
summary(M_1_2)
plot(M_1_2)

vP1_3 <- varExp(form=~Latitude+dem_ele+rock_type_int+MAP_chelsa+range_tmp_chelsa)
M_1_3=gls(composition_genera ~ Latitude+dem_ele+rock_type_int+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP1_3)
summary(M_1_3)
plot(M_1_3)

vP1_4 <- varIdent(form=~Latitude+dem_ele+soil_thick_log+MAP_chelsa+range_tmp_chelsa)
M_1_4=gls(composition_genera ~ Latitude+dem_ele+soil_thick_log+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP1_4)
summary(M_1_4)
plot(M_1_4)

vP1_5 <- varExp(form=~Latitude+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa)
M_1_5=gls(composition_genera ~ Latitude+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP1_5)
summary(M_1_5)
plot(M_1_5)

vP1_6 <- varPower(form=~dem_ele+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa)
M_1_6=gls(composition_genera ~ dem_ele+rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP1_6)
summary(M_1_6)
plot(M_1_6)


###Using 4 variables at a time###

vP2_1 <- varIdent(form=~Latitude+dem_ele+rock_type_int+soil_thick_log)
M_2_1=gls(composition_genera ~ Latitude+dem_ele+rock_type_int+soil_thick_log, data=ecovar_data_GLS_no_categories, 
          weights=vP2_1)
summary(M_2_1)
plot(M_2_1)

vP2_2 <- varPower(form=~Latitude+dem_ele+rock_type_int+MAP_chelsa)
M_2_2=gls(composition_genera ~ Latitude+dem_ele+rock_type_int+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP2_2)
summary(M_2_2)
plot(M_2_2)

vP2_3 <- varPower(form=~Latitude+dem_ele+rock_type_int+range_tmp_chelsa)
M_2_3=gls(composition_genera ~ Latitude+dem_ele+rock_type_int+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP2_3)
summary(M_2_3)
plot(M_2_3)

vP2_4 <- varPower(form=~Latitude+dem_ele+soil_thick_log+MAP_chelsa)
M_2_4=gls(composition_genera ~ Latitude+dem_ele+soil_thick_log+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP2_4)
summary(M_2_4)
plot(M_2_4)

vP2_5 <- varPower(form=~Latitude+dem_ele+soil_thick_log+range_tmp_chelsa)
M_2_5=gls(composition_genera ~ Latitude+dem_ele+soil_thick_log+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP2_5)
summary(M_2_5)
plot(M_2_5)

vP2_6 <- varPower(form=~Latitude+dem_ele+MAP_chelsa+range_tmp_chelsa)
M_2_6=gls(composition_genera ~ Latitude+dem_ele+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP2_6)
summary(M_2_6)
plot(M_2_6)

vP2_7 <- varIdent(form=~Latitude+rock_type_int+soil_thick_log+MAP_chelsa)
M_2_7=gls(composition_genera ~ Latitude+rock_type_int+soil_thick_log+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP2_7)
summary(M_2_7)
plot(M_2_7)

vP2_8 <- varPower(form=~Latitude+rock_type_int+soil_thick_log+range_tmp_chelsa)
M_2_8=gls(composition_genera ~ Latitude+rock_type_int+soil_thick_log+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP2_8)
summary(M_2_8)
plot(M_2_8)

vP2_9 <- varIdent(form=~Latitude+rock_type_int+MAP_chelsa+range_tmp_chelsa)
M_2_9=gls(composition_genera ~ Latitude+rock_type_int+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP2_9)
summary(M_2_9)
plot(M_2_9)

vP2_10 <- varIdent(form=~Latitude+soil_thick_log+MAP_chelsa+range_tmp_chelsa)
M_2_10=gls(composition_genera ~ Latitude+soil_thick_log+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP2_10)
summary(M_2_10)
plot(M_2_10)

vP2_11 <- varPower(form=~dem_ele+rock_type_int+soil_thick_log+MAP_chelsa)
M_2_11=gls(composition_genera ~ dem_ele+rock_type_int+soil_thick_log+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP2_11)
summary(M_2_11)
plot(M_2_11)

vP2_12 <- varPower(form=~dem_ele+rock_type_int+soil_thick_log+range_tmp_chelsa)
M_2_12=gls(composition_genera ~ dem_ele+rock_type_int+soil_thick_log+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP2_12)
summary(M_2_12)
plot(M_2_12)

vP2_13 <- varPower(form=~dem_ele+rock_type_int+MAP_chelsa+range_tmp_chelsa)
M_2_13=gls(composition_genera ~ dem_ele+rock_type_int+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP2_13)
summary(M_2_13)
plot(M_2_13)

vP2_14 <- varPower(form=~dem_ele+soil_thick_log+MAP_chelsa+range_tmp_chelsa)
M_2_14=gls(composition_genera ~ dem_ele+soil_thick_log+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP2_14)
summary(M_2_14)
plot(M_2_14)

vP2_15 <- varPower(form=~rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa)
M_2_15=gls(composition_genera ~ rock_type_int+soil_thick_log+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP2_15)
summary(M_2_15)
plot(M_2_15)


###Using 3 variables at a time###

vP3_1 <- varPower(form=~Latitude+dem_ele+rock_type_int)
M_3_1=gls(composition_genera ~ Latitude+dem_ele+rock_type_int, data=ecovar_data_GLS_no_categories, 
          weights=vP3_1)
summary(M_3_1)
plot(M_3_1)

vP3_2 <- varIdent(form=~Latitude+dem_ele+soil_thick_log)
M_3_2=gls(composition_genera ~ Latitude+dem_ele+soil_thick_log, data=ecovar_data_GLS_no_categories, 
          weights=vP3_2)
summary(M_3_2)
plot(M_3_2)

vP3_3 <- varPower(form=~Latitude+dem_ele+MAP_chelsa)
M_3_3=gls(composition_genera ~ Latitude+dem_ele+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP3_3)
summary(M_3_3)
plot(M_3_3)

vP3_4 <- varPower(form=~Latitude+dem_ele+range_tmp_chelsa)
M_3_4=gls(composition_genera ~ Latitude+dem_ele+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP3_4)
summary(M_3_4)
plot(M_3_4)

vP3_5 <- varPower(form=~Latitude+rock_type_int+soil_thick_log)
M_3_5=gls(composition_genera ~ Latitude+rock_type_int+soil_thick_log, data=ecovar_data_GLS_no_categories, 
          weights=vP3_5)
summary(M_3_5)
plot(M_3_5)

vP3_6 <- varIdent(form=~Latitude+rock_type_int+MAP_chelsa)
M_3_6=gls(composition_genera ~ Latitude+rock_type_int+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP3_6)
summary(M_3_6)
plot(M_3_6)

vP3_7 <- varPower(form=~Latitude+rock_type_int+range_tmp_chelsa)
M_3_7=gls(composition_genera ~ Latitude+rock_type_int+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP3_7)
summary(M_3_7)
plot(M_3_7)

vP3_8 <- varIdent(form=~Latitude+soil_thick_log+MAP_chelsa)
M_3_8=gls(composition_genera ~ Latitude+soil_thick_log+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP3_8)
summary(M_3_8)
plot(M_3_8)

vP3_9 <- varIdent(form=~Latitude+soil_thick_log+range_tmp_chelsa)
M_3_9=gls(composition_genera ~ Latitude+soil_thick_log+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP3_9)
summary(M_3_9)
plot(M_3_9)

vP3_10 <- varIdent(form=~Latitude+MAP_chelsa+range_tmp_chelsa)
M_3_10=gls(composition_genera ~ Latitude+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_10)
summary(M_3_10)
plot(M_3_10)

vP3_11 <- varPower(form=~dem_ele+rock_type_int+soil_thick_log)
M_3_11=gls(composition_genera ~ dem_ele+rock_type_int+soil_thick_log, data=ecovar_data_GLS_no_categories, 
           weights=vP3_11)
summary(M_3_11)
plot(M_3_11)

vP3_12 <- varPower(form=~dem_ele+rock_type_int+MAP_chelsa)
M_3_12=gls(composition_genera ~ dem_ele+rock_type_int+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_12)
summary(M_3_12)
plot(M_3_12)

vP3_13 <- varPower(form=~dem_ele+rock_type_int+range_tmp_chelsa)
M_3_13=gls(composition_genera ~ dem_ele+rock_type_int+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_13)
summary(M_3_13)
plot(M_3_13)

vP3_14 <- varPower(form=~dem_ele+soil_thick_log+MAP_chelsa)
M_3_14=gls(composition_genera ~ dem_ele+soil_thick_log+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_14)
summary(M_3_14)
plot(M_3_14)

vP3_15 <- varIdent(form=~dem_ele+soil_thick_log+range_tmp_chelsa)
M_3_15=gls(composition_genera ~ dem_ele+soil_thick_log+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_15)
summary(M_3_15)
plot(M_3_15)

vP3_16 <- varPower(form=~dem_ele+MAP_chelsa+range_tmp_chelsa)
M_3_16=gls(composition_genera ~ dem_ele+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_16)
summary(M_3_16)
plot(M_3_16)

vP3_17 <- varIdent(form=~rock_type_int+soil_thick_log+MAP_chelsa)
M_3_17=gls(composition_genera ~ rock_type_int+soil_thick_log+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_17)
summary(M_3_17)
plot(M_3_17)

vP3_18 <- varIdent(form=~rock_type_int+soil_thick_log+range_tmp_chelsa)
M_3_18=gls(composition_genera ~ rock_type_int+soil_thick_log+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_18)
summary(M_3_18)
plot(M_3_18)

vP3_19 <- varIdent(form=~rock_type_int+MAP_chelsa+range_tmp_chelsa)
M_3_19=gls(composition_genera ~ rock_type_int+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_19)
summary(M_3_19)
plot(M_3_19)

vP3_20 <- varIdent(form=~soil_thick_log+MAP_chelsa+range_tmp_chelsa)
M_3_20=gls(composition_genera ~ soil_thick_log+MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP3_20)
summary(M_3_20)
plot(M_3_20)


###Using 2 variables at a time###

vP4_1 <- varPower(form=~Latitude+dem_ele)
M_4_1=gls(composition_genera ~ Latitude+dem_ele, data=ecovar_data_GLS_no_categories, 
          weights=vP4_1)
summary(M_4_1)
plot(M_4_1)

vP4_2 <- varPower(form=~Latitude+rock_type_int)
M_4_2=gls(composition_genera ~ Latitude+rock_type_int, data=ecovar_data_GLS_no_categories, 
          weights=vP4_2)
summary(M_4_2)
plot(M_4_2)

vP4_3 <- varPower(form=~Latitude+soil_thick_log)
M_4_3=gls(composition_genera ~ Latitude+soil_thick_log, data=ecovar_data_GLS_no_categories, 
          weights=vP4_3)
summary(M_4_3)
plot(M_4_3)

vP4_4 <- varIdent(form=~Latitude+MAP_chelsa)
M_4_4=gls(composition_genera ~ Latitude+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP4_4)
summary(M_4_4)
plot(M_4_4)

vP4_5 <- varPower(form=~Latitude+range_tmp_chelsa)
M_4_5=gls(composition_genera ~ Latitude+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP4_5)
summary(M_4_5)
plot(M_4_5)

vP4_6 <- varPower(form=~dem_ele+rock_type_int)
M_4_6=gls(composition_genera ~ dem_ele+rock_type_int, data=ecovar_data_GLS_no_categories, 
          weights=vP4_6)
summary(M_4_6)
plot(M_4_6)

vP4_7 <- varIdent(form=~dem_ele+soil_thick_log)
M_4_7=gls(composition_genera ~ dem_ele+soil_thick_log, data=ecovar_data_GLS_no_categories, 
          weights=vP4_7)
summary(M_4_7)
plot(M_4_7)

vP4_8 <- varIdent(form=~dem_ele+MAP_chelsa)
M_4_8=gls(composition_genera ~ dem_ele+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP4_8)
summary(M_4_8)
plot(M_4_8)

vP4_9 <- varPower(form=~dem_ele+range_tmp_chelsa)
M_4_9=gls(composition_genera ~ dem_ele+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
          weights=vP4_9)
summary(M_4_9)
plot(M_4_9)

vP4_10 <- varIdent(form=~rock_type_int+soil_thick_log)
M_4_10=gls(composition_genera ~ rock_type_int+soil_thick_log, data=ecovar_data_GLS_no_categories, 
           weights=vP4_10)
summary(M_4_10)
plot(M_4_10)

vP4_11 <- varIdent(form=~rock_type_int+MAP_chelsa)
M_4_11=gls(composition_genera ~ rock_type_int+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP4_11)
summary(M_4_11)
plot(M_4_11)

vP4_12 <- varPower(form=~rock_type_int+range_tmp_chelsa)
M_4_12=gls(composition_genera ~ rock_type_int+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP4_12)
summary(M_4_12)
plot(M_4_12)

vP4_13 <- varExp(form=~soil_thick_log+MAP_chelsa)
M_4_13=gls(composition_genera ~ soil_thick_log+MAP_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP4_13)
summary(M_4_13)
plot(M_4_13)

vP4_14 <- varIdent(form=~soil_thick_log+range_tmp_chelsa)
M_4_14=gls(composition_genera ~ soil_thick_log+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP4_14)
summary(M_4_14)
plot(M_4_14)

vP4_15 <- varIdent(form=~MAP_chelsa+range_tmp_chelsa)
M_4_15=gls(composition_genera ~ MAP_chelsa+range_tmp_chelsa, data=ecovar_data_GLS_no_categories, 
           weights=vP4_15)
summary(M_4_15)
plot(M_4_15)


#Best Model selection 

#best_Model_prediction<- function(Latitude, MAP_chelsa){
  
#  composition_genera<- -3.944688-0.204889* Latitude + 0.007444 *MAP_chelsa 
  
#  return(composition_genera)
#}

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


#Predicitons


x_pred = read_parquet("/data/x_pred.parquet")

mod_pred<-predict(M_4_4, x_pred)

write.csv(mod_pred, file='../results/pred_gls_richness.csv')
