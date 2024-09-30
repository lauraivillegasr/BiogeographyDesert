#### Reproductive mode GLM with aditional analysis incorporating vegetation information #### 


asex_sex_data_GLM2 <- read.csv2("/data/asex_sex_data_GLM_nocategory.csv")

#Test the different models, null model now includes vegetation coverage around 100m and percentage of phanerophytes in sampling spot

glM0<-glm(Reproductive.mode ~ dem_ele+MAP_chelsa+Latitude+veg_cover_100m+per_Phanerophytes, family=binomial(link="logit"), 
          data=asex_sex_data_GLM2)
summary(glM0)


glM1<-glm(Reproductive.mode ~ dem_ele, family=binomial(link="logit"), 
          data=asex_sex_data_GLM2)
summary(glM1)


glM2<-glm(Reproductive.mode ~ MAP_chelsa, family=binomial(link="logit"), 
          data=asex_sex_data_GLM2)
summary(glM2)


glM3<-glm(Reproductive.mode ~ Latitude, family=binomial(link="logit"), 
          data=asex_sex_data_GLM2)
summary(glM3)

glM4<-glm(Reproductive.mode ~ Latitude+dem_ele, family=binomial(link="logit"), 
          data=asex_sex_data_GLM2)
summary(glM4)


glM5<-glm(Reproductive.mode ~ veg_cover_100m, family=binomial(link="logit"), 
          data=asex_sex_data_GLM2)
summary(glM5)


glM6<-glm(Reproductive.mode ~ per_Phanerophytes, family=binomial(link="logit"), 
          data=asex_sex_data_GLM2)
summary(glM6)



BICsglm<-sapply(mget(ls(pattern="glM*"),.GlobalEnv),BIC)
Dbic<-c()
Wm<-c()
