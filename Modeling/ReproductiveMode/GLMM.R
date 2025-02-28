##GLM reproductive mode ##
#elevation and latitude are tested according to geographical parthenogenesis hypothesis,
#MAP is tested given its effect on predicting genera richness 
#Range of temperature is tested given its effect on predicting genera richness 


library(arrow)

#load necessary data

asex_sex_data_GLM2 <- read.csv2("/data/EcoVarData_reprod_transect_testing_fixed.csv")

# Load necessary library
library(lme4)

# Model 1: Full model with all predictors
glM0_glmer <- glmer(Reproductive.mode ~ dem_ele + MAP_chelsa + Latitude + range_tmp_chelsa + (1 | Transect), 
                    family = binomial(link = "logit"), 
                    data = asex_sex_data_GLM2)

summary(glM0_glmer)

# Model 2: Model with only dem_ele
glM1_glmer <- glmer(Reproductive.mode ~ dem_ele + (1 | Transect), 
                    family = binomial(link = "logit"), 
                    data = asex_sex_data_GLM2)
summary(glM1_glmer)


# Model 3: Model with only MAP_chelsa
glM2_glmer <- glmer(Reproductive.mode ~ MAP_chelsa + (1 | Transect), 
                    family = binomial(link = "logit"), 
                    data = asex_sex_data_GLM2)
summary(glM2_glmer)

# Model 4: Model with only Latitude
glM3_glmer <- glmer(Reproductive.mode ~ Latitude + (1 | Transect), 
                    family = binomial(link = "logit"), 
                    data = asex_sex_data_GLM2)
summary(glM3_glmer)

# Model 5: Model with Latitude + dem_ele
glM4_glmer <- glmer(Reproductive.mode ~ Latitude + dem_ele + (1 | Transect), 
                    family = binomial(link = "logit"), 
                    data = asex_sex_data_GLM2)

summary(glM4_glmer)

# Model 6: Model with range of temperature
glM5_glmer <- glmer(Reproductive.mode ~ range_tmp_chelsa + (1 | Transect), 
                    family = binomial(link = "logit"), 
                    data = asex_sex_data_GLM2)

summary(glM5_glmer)


BICsglm<-sapply(mget(ls(pattern="glM*"),.GlobalEnv),BIC)
Dbic<-c()
Wm<-c()

for(i in 1:length(BICsglm)){
  tBIC=BICsglm[i]-min(BICsglm)
  Dbic<-c(Dbic,tBIC)
}

for(p in 1:length(BICsglm)){
  tWm=exp(-Dbic[p]/2)/sum(exp(-Dbic/2))
  Wm<-c(Wm,tWm)
}

BICsglm
Wm
max_Wm_value <- max(Wm)
max_Wm_value
min_Wm_model_name <- names(Wm)[which.max(Wm)]
min_Wm_model_name
