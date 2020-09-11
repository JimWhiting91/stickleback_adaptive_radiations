#######################################################################################################################
# This script compares compares environmental measures to phenotypes to explore env x pheno associations...
#######################################################################################################################

# Load packages
lib<-as.vector(c("multcomp","effects","lmSupport","MuMIn","lsmeans","ggsignif","MASS","data.table","patchwork","VennDiagram","ggplot2","parallel","effects","lme4","lmerTest"))
lapply(lib,library,character.only=TRUE)

# Read in data
dd<-read.csv("data/All_Lakes_Env_and_Pheno_v2.csv")
colnames(dd)[1]<-"RADIATION"

##### Calculate the degree of change for each variable, for each radiation ######

# List out the variables
variable_vec<-c("Ca", "Gyro",    "Na","pH",   "Schisto","Zn",  "Lake_Area",
                "SHAPE.IND.PC1",            "SHAPE.IND.PC2",            "SHAPE.IND.PC3",
                "Avg.resDS1",  "Avg.resDS2",  "Avg.resPS" ,  "Avg.resLP",   "Avg.resHP",   "Avg.resBAP" , "Avg.n.plate",
                "Gill.Raker.N", "Resid.Raker.L")

colnames(dd)[colnames(dd)=="lake_area_km"]<-"Lake_Area"

group_list<-list(variable_vec[1:7],variable_vec[8:10],variable_vec[11:17],variable_vec[8:19])
group_vec<-c("ENV","SHAPE","ARMOUR","PHENO")
rad_vec<-c("ALASKA","BC","ICELAND","SCOTLAND")


# Run through each phenotype variable, build models per variable ----------
phenos <- c("ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2",
            "SHAPE.IND.PC1","SHAPE.IND.PC2",
            "Gill.Raker.N","Resid.Raker.L")
envs <- variable_vec[1:7]

# Make the models
model_out <- data.frame(rbindlist(lapply(phenos,function(pheno){
print(pheno)
  # Run by each env variable...
  env_out <- lapply(envs,function(env){
    print(env)
    tmp<-data.frame(var=dd[,pheno],
                    rad=dd$RADIATION,
                    env=dd[,env])
    
    # Build each of the models
    null <- glm(var ~ 1,tmp,family="gaussian")
    mod1 <- glm(var ~ rad,tmp,family="gaussian")
    mod2 <- glm(var ~ env,tmp,family="gaussian")
    mod3 <- glm(var ~ env + rad,tmp,family="gaussian")
    mod4 <- glm(var ~ env * rad,tmp,family="gaussian")
    model_list<-list(null,mod1,mod2,mod3,mod4)
    
    # Get AIC and reward low complexity...
    aic_vec<-AIC(null,mod1,mod2,mod3,mod4)

    # Manually chooe best model...
    best_model<-5
    for(i in 4:1){
      if(aic_vec$AIC[i] - aic_vec$AIC[best_model] < 2){
      best_model<-i
      }
    }

    # Calculate model...
    res <- anova(model_list[[best_model]],test="F")
    
    # Tidy up
    rownames(res)<-gsub("rad","Radiation",rownames(res))
    rownames(res)<-gsub("env",env,rownames(res))
    
    # Write the model basics
    out1<-matrix(nrow=1,ncol=7)
    out1[1,1]<-pheno
    out1[1,2]<-env
    for(i in 1:length(aic_vec$AIC)){
      out1[1,i+2]<-round(aic_vec$AIC[i],1)
    }
    
    # We don't write the results for null models...
    if(best_model != 1){
    
    # Write the model res
    out2<-matrix(nrow=nrow(res)-1,ncol=5)
    for(i in 1:nrow(out2)){
      out2[i,1]<-rownames(res)[i+1]
      out2[i,2]<-paste0(res$Df[i+1],",",res$`Resid. Df`[i+1])
      out2[i,3]<-round(res$F[i+1],3)
      out2[i,4]<-as.numeric(res$`Pr(>F)`[i+1])
      out2[i,5]<-round(r.squaredGLMM(model_list[[best_model]])[1],3)
    }
    } else {
      out2<-matrix(NA,nrow=1,ncol=5)
    }
    
    # Fill blanks
    if(nrow(out1) != nrow(out2)){
      diff<-nrow(out2)-nrow(out1)
      for(i in 1:diff){
        out1<-rbind(out1,
                    matrix("",nrow=1,ncol=7))
      }
    }
    
    # cbind them
    out<-cbind(out1,out2)
    
    # Convert to data/frame and save
    out_dd<-data.frame(out)
    colnames(out_dd)<-c("Phenotype","Environment","Null","Intercept (Radiation)","Slope (Env)","Intercept + Slope","Intercept * Slope",
                        "Best Model Variables","df","F","p")
    return(data.frame(out))
})
  
# Merge them
return(data.frame(rbindlist(env_out)))
  
})))

# Save these to a csv
colnames(model_out)<-c("Phenotype","Environment","Null","Intercept (Radiation)","Slope (Env)","Intercept + Slope","Intercept * Slope",
                    "Best Model Variables","df","F","p","Variance Explained")

# Add an fdr
model_out$fdr<-p.adjust(model_out$p,method = "fdr")

# Tidy up
new_phenos<-c("Armour PC1","Armour PC2","Shape PC1","Shape PC2","Gill Raker N","Gill Raker L")
for(i in 1:length(new_phenos)){
model_out$Phenotype<-gsub(phenos[i],new_phenos[i],model_out$Phenotype)
}

write.csv(model_out,"tables/TableSX_All_EnvXPheno_models_AIC_based.csv",
          row.names = F)


# Describing models ------------------------------------------------
# Armour PC1 - Calcium
ca <- glm(ARMOUR.IND.WITHPLATE.PC1 ~ Ca*RADIATION,dd,family="gaussian")
plot(allEffects(ca))
pairs(lstrends(ca, "RADIATION", var="Ca"))

# Armour PC1 - Sodium
na <- glm(ARMOUR.IND.WITHPLATE.PC1 ~ Na*RADIATION,dd,family="gaussian")
plot(allEffects(na))
pairs(lstrends(na, "RADIATION", var="Na"))

# Armour PC1 - Lake Area
area <- glm(ARMOUR.IND.WITHPLATE.PC1 ~ Lake_Area*RADIATION,dd,family="gaussian")
plot(allEffects(area))
pairs(lstrends(area, "RADIATION", var="Lake_Area"))

# Armour PC2 - Calcium
ca <- glm(ARMOUR.IND.WITHPLATE.PC2 ~ Ca*RADIATION,dd,family="gaussian")
plot(allEffects(ca))
pairs(lstrends(ca, "RADIATION", var="Ca"))

# Armour PC2 - Sodium
na <- glm(ARMOUR.IND.WITHPLATE.PC2 ~ Na*RADIATION,dd,family="gaussian")
plot(allEffects(na))
pairs(lstrends(na, "RADIATION", var="Na"))

# Armour PC2 - pH
pH <- glm(ARMOUR.IND.WITHPLATE.PC2 ~ pH*RADIATION,dd,family="gaussian")
plot(allEffects(pH))
pairs(lstrends(area, "RADIATION", var="pH"))

########################################
# Shape PC1 - Gyro
gyro <- glm(SHAPE.IND.PC1 ~ Gyro*RADIATION,dd,family="gaussian")
plot(allEffects(gyro))
pairs(lstrends(gyro, "RADIATION", var="Gyro"))

# Shape PC1 - Calcium
na <- glm(SHAPE.IND.PC1 ~ Na + RADIATION,dd,family="gaussian")
plot(allEffects(na))
pairs(lstrends(na, "RADIATION", var="Na"))

# Shape PC2 - Schisto
schisto <- glm(SHAPE.IND.PC2 ~ Schisto + RADIATION,dd,family="gaussian")
plot(allEffects(schisto))
pairs(lstrends(schisto, "RADIATION", var="Gyro"))

# Shape PC2 - Calcium
zn <- glm(SHAPE.IND.PC2 ~ Zn * RADIATION,dd,family="gaussian")
plot(allEffects(zn))
pairs(lstrends(zn, "RADIATION", var="Zn"))

########################################
zn <- glm(Gill.Raker.N ~ Zn + RADIATION,dd,family="gaussian")
plot(allEffects(zn))

# Visualise relationship within each radiation...
ggplot(dd,aes(x=Zn,y=Gill.Raker.N,colour=RADIATION))+
  geom_smooth(method="lm")+
  geom_point()




