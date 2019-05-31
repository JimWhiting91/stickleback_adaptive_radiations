#################
# This script calculates effect sizes for genomic parallelism and environental/phenotypic similarity
# It then models these to investigate the relationships alongside radiation-level FST
#################


lib<-as.vector(c("data.table","VennDiagram","ggplot2","parallel","effsize","lme4","lmerTest"))
lapply(lib,library,character.only=TRUE)

# ------
# This script will be repeated over windows and genes so we put in function
# ------

# Define the grouping
grouping_vector<-as.vector(c("50k","75k","100k","200k","genes"))
# Define the variables
variables_vector<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn",
                              "Shape_PC1","Shape_PC2","Shape_PC3",
                              "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                              "Gill_Raker_L","Gill_Raker_N"))
# Define radiations
radiations_vector<-as.vector(c("Alaska","BC","Iceland","Nuist"))

# Prior to analysis, we can clean and edit the Env and Pheno data as it will be consistent across window sizes
env_dd<-read.csv("data/All_Lakes_Env_and_Pheno_v2.csv")

# Tidy up
env_dd2<-data.frame(env_dd$Ca,env_dd$Gyro,env_dd$Na,env_dd$pH,env_dd$Schisto,env_dd$Zn,
                    env_dd$SHAPE.IND.PC1,env_dd$SHAPE.IND.PC2,env_dd$SHAPE.IND.PC3,
                    env_dd$Avg.resDS1,env_dd$Avg.resDS2,env_dd$Avg.resPS,env_dd$Avg.resLP,env_dd$Avg.resHP,env_dd$Avg.resBAP,env_dd$Avg.n.plate,
                    env_dd$Gill.Raker.N,env_dd$Resid.Raker.L)

# Scale variables
env_dd2<-data.frame(cbind(as.character(env_dd[,1]),as.character(env_dd$LAKE),
                          scale(env_dd2)))
colnames(env_dd2)<-c("RADIATION","LAKE",variables_vector)

# For each of the 18 variables we want the effect size for each 2 way comparison
env_pheno_cohen_calcs<-lapply(3:20,function(i){
  cohen_calc_AB<-abs(mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ALASKA",i])))-mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "BC",i]))))/(((sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ALASKA",i])))^2+sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "BC",i])))^2)/2)^0.5)
  cohen_calc_AI<-abs(mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ALASKA",i])))-mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ICELAND",i]))))/(((sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ALASKA",i])))^2+sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ICELAND",i])))^2)/2)^0.5)
  cohen_calc_AN<-abs(mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ALASKA",i])))-mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "SCOTLAND",i]))))/(((sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ALASKA",i])))^2+sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "SCOTLAND",i])))^2)/2)^0.5)
  cohen_calc_BI<-abs(mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "BC",i])))-mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ICELAND",i]))))/(((sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "BC",i])))^2+sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ICELAND",i])))^2)/2)^0.5)
  cohen_calc_BN<-abs(mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "BC",i])))-mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "SCOTLAND",i]))))/(((sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "BC",i])))^2+sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "SCOTLAND",i])))^2)/2)^0.5)
  cohen_calc_IN<-abs(mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ICELAND",i])))-mean(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "SCOTLAND",i]))))/(((sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "ICELAND",i])))^2+sd(as.numeric(as.character(env_dd2[env_dd2$RADIATION == "SCOTLAND",i])))^2)/2)^0.5)

  cohen_vec<-abs(data.table(t(c(cohen_calc_AB,
                                cohen_calc_AI,
                                cohen_calc_AN,
                                cohen_calc_BI,
                                cohen_calc_BN,
                                cohen_calc_IN))))
})

# rbindlist to get all cohen scores
env_pheno_cohen<-data.frame(rbindlist(env_pheno_cohen_calcs))
rownames(env_pheno_cohen)<-variables_vector
colnames(env_pheno_cohen)<-c("AB","AI","AN","BI","BN","IN")

# Convert Matrix into data frame to be used in models by appending each row
clean_env_cohen<-data.frame(rbindlist(lapply(1:length(colnames(env_pheno_cohen)),function(i){
  tmp_dd<-data.frame(env_pheno_cohen[,i])
  colnames(tmp_dd)<-"Env_Eff"
  tmp_dd$var<-rownames(env_pheno_cohen)
  tmp_dd$comparison<-rep(colnames(env_pheno_cohen)[i],length(rownames(tmp_dd)))
  data.table(tmp_dd)
})))

# We also want the absolute difference in coefficient of the regression to gauge slope
slope_calcs<-lapply(3:20,function(x){
  comparisons_list<-list(c("ALASKA","BC"),
                         c("ALASKA","ICELAND"),
                         c("ALASKA","SCOTLAND"),
                         c("BC","ICELAND"),
                         c("BC","SCOTLAND"),
                         c("ICELAND","SCOTLAND"))
  # We make models for each comparison and acquire the coefficient of the gradient for sorted env/pheno variables
  for (i in 1:length(comparisons_list)){
    vec1<-sort(as.numeric(as.character(env_dd2[env_dd2$RADIATION == comparisons_list[[i]][1],x])))
    vec2<-sort(as.numeric(as.character(env_dd2[env_dd2$RADIATION == comparisons_list[[i]][2],x])))
    length1<-1:length(vec1)
    length2<-1:length(vec2)
    
    mod1<-lm(vec1~length1)
    mod2<-lm(vec2~length2)
    slope_diff<-abs(mod1$coefficients[2]-mod2$coefficients[2])
    assign(paste0("slope_calc_",i),slope_diff)
  }
  slope_vec<-data.table(t(c(slope_calc_1,
                            slope_calc_2,
                            slope_calc_3,
                            slope_calc_4,
                            slope_calc_5,
                            slope_calc_6)))
  
})

# rbindlist to get all slopes scores
env_pheno_slopes<-data.frame(rbindlist(slope_calcs))
rownames(env_pheno_slopes)<-variables_vector
colnames(env_pheno_slopes)<-c("AB","AI","AN","BI","BN","IN")

# Convert Matrix into data frame to be used in models by appending each row
clean_slope<-data.frame(rbindlist(lapply(1:length(colnames(env_pheno_slopes)),function(i){
  tmp_dd<-data.frame(env_pheno_slopes[,i])
  colnames(tmp_dd)<-"Slope_Eff"
  tmp_dd$var<-rownames(env_pheno_slopes)
  tmp_dd$comparison<-rep(colnames(env_pheno_slopes)[i],length(rownames(tmp_dd)))
  data.table(tmp_dd)
})))

# Cbind Slope scores into clean_env_cohen
clean_env_cohen<-cbind(clean_env_cohen,clean_slope$Slope_Eff)
colnames(clean_env_cohen)[length(colnames(clean_env_cohen))]<-"Slope_Eff"

#-------
# We are also interested in environmental and phenotypic theta and deltaL
#-------
# Read in from outputs
theta_delta<-read.table("outputs/theta_deltaL_for_all_pairings_PC-based_scaled.txt",header=T,sep="\t")

# Tidy and get into right format
thetas<-NULL
for(i in 1:6){
  thetas<-unlist(list(thetas,
                      c(rep(theta_delta$ENV_theta[i],6),
                               rep(theta_delta$SHAPE_theta[i],3),
                               rep(theta_delta$ARMOUR_theta[i],7),
                               rep(theta_delta$GILL_theta[i],2))))
}
clean_env_cohen$theta<-thetas

# Scale deltaL because they are not comparable across variables
theta_delta$ENV_deltaL_scaled<-scale(theta_delta$ENV_deltaL)
theta_delta$SHAPE_deltaL_scaled<-scale(theta_delta$SHAPE_deltaL)
theta_delta$ARMOUR_deltaL_scaled<-scale(theta_delta$ARMOUR_deltaL)
theta_delta$GILL_deltaL_scaled<-scale(theta_delta$GILL_deltaL)

deltas<-NULL
for(i in 1:6){
  deltas<-unlist(list(deltas,
                      c(rep(theta_delta$ENV_deltaL_scaled[i],6),
                        rep(theta_delta$SHAPE_deltaL_scaled[i],3),
                        rep(theta_delta$ARMOUR_deltaL_scaled[i],7),
                        rep(theta_delta$GILL_deltaL_scaled[i],2))))
}
clean_env_cohen$deltaL<-deltas


# Function will be run over 50kb windows
x<-1
  
  # Calculte Cohen's D for genomic parallelism across variables
  para_cohen_list<-lapply(1:length(variables_vector),function(y){
    
    # For each variable, we read in the simulated output file from 03
    sim_dd<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_10000runs_FULLSIM.txt"),header=T)
    
    # We are only interested in 2-way comparisons, so we need to take 2-way comparisons and sum across 3 and 4 ways
    sim_dd$AB_all<-sim_dd$ABIN+sim_dd$ABI+sim_dd$ABN+sim_dd$AB
    sim_dd$AI_all<-sim_dd$ABIN+sim_dd$ABI+sim_dd$AIN+sim_dd$AI
    sim_dd$AN_all<-sim_dd$ABIN+sim_dd$AIN+sim_dd$ABN+sim_dd$AN
    sim_dd$BI_all<-sim_dd$ABIN+sim_dd$ABI+sim_dd$BIN+sim_dd$BI
    sim_dd$BN_all<-sim_dd$ABIN+sim_dd$ABN+sim_dd$BIN+sim_dd$BN
    sim_dd$IN_all<-sim_dd$ABIN+sim_dd$AIN+sim_dd$BIN+sim_dd$IN
    
    # We also want to read in the observed counts
    obs_dd<-read.table(paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_",variables_vector[y],"Venn_overlaps.txt"),header=T)
    colnames(obs_dd)<-colnames(sim_dd)[1:15]
    
    # 2 way comparisons again
    obs_dd$AB_all<-obs_dd$ABIN+obs_dd$ABI+obs_dd$ABN+obs_dd$AB
    obs_dd$AI_all<-obs_dd$ABIN+obs_dd$ABI+obs_dd$AIN+obs_dd$AI
    obs_dd$AN_all<-obs_dd$ABIN+obs_dd$AIN+obs_dd$ABN+obs_dd$AN
    obs_dd$BI_all<-obs_dd$ABIN+obs_dd$ABI+obs_dd$BIN+obs_dd$BI
    obs_dd$BN_all<-obs_dd$ABIN+obs_dd$ABN+obs_dd$BIN+obs_dd$BN
    obs_dd$IN_all<-obs_dd$ABIN+obs_dd$AIN+obs_dd$BIN+obs_dd$IN
    
    # So now we want to calculate cohen's D for effect sizes
    cohen_d<-vector()
    for (i in 16:21){
      cohen_d[i-15]<-abs((obs_dd[,i]-mean(sim_dd[,i]))/sd(sim_dd[,i]))
    }
    cohen_d
  })
  
  # Make dataframe from output
  para_cohen_dd<-data.frame(matrix(ncol = 6,nrow=length(variables_vector)))
  for(i in 1:nrow(para_cohen_dd)){
    para_cohen_dd[i,]<-para_cohen_list[[i]]
  }
  rownames(para_cohen_dd)<-variables_vector
  colnames(para_cohen_dd)<-c("AB","AI","AN","BI","BN","IN")
  
  # Convert data frame to be used in models by appending each row
  clean_para_cohen<-data.frame(rbindlist(lapply(1:length(colnames(para_cohen_dd)),function(i){
    tmp_dd<-data.frame(para_cohen_dd[,i])
    colnames(tmp_dd)<-"Para_Eff"
    tmp_dd$var<-rownames(para_cohen_dd)
    tmp_dd$comparison<-rep(colnames(para_cohen_dd)[i],length(rownames(tmp_dd)))
    data.table(tmp_dd)
  })))
  
  # Now we can combine parallelism effect sizes with env/pheno
  eff_test_dd<-cbind(clean_env_cohen,clean_para_cohen$Para_Eff)
  colnames(eff_test_dd)[length(eff_test_dd)]<-"Para_Eff"
  
  # We also want to include FST in the data frame to model genetic differentiation
  eff_test_dd$FST<-rep(0,length(eff_test_dd$Env_Eff))
  eff_test_dd[eff_test_dd$comparison == "AB","FST"]<-0.19803061374269
  eff_test_dd[eff_test_dd$comparison == "AI","FST"]<-0.314056804093567
  eff_test_dd[eff_test_dd$comparison == "AN","FST"]<-0.329127038011696
  eff_test_dd[eff_test_dd$comparison == "BI","FST"]<-0.322677358024691
  eff_test_dd[eff_test_dd$comparison == "BN","FST"]<-0.337744700617284
  eff_test_dd[eff_test_dd$comparison == "IN","FST"]<-0.193848475617284
  
  # Add a column to specify which variables are environmental and which are phenotypic
  eff_test_dd$env_pheno<-rep(0,length(eff_test_dd$Env_Eff))
  eff_test_dd[eff_test_dd$var %in% variables_vector[1:6],"env_pheno"]<-"Env"
  eff_test_dd[eff_test_dd$var %in% variables_vector[7:18],"env_pheno"]<-"Pheno"
  
  # Specify intra-inter continent groupings
  eff_test_dd$continent<-"Inter"
  eff_test_dd[eff_test_dd$comparison %in% c("AB","IN"),"continent"]<-"Intra"
  
  # Also bring in just env_theta
  eff_test_dd$env_theta<-NA
  for(i in 1:nrow(eff_test_dd)){
    eff_test_dd$env_theta[i]<-eff_test_dd[eff_test_dd$comparison == eff_test_dd$comparison[i] & eff_test_dd$env_pheno == "Env","theta"][1]
  }
  
  # We are also interested in a binomial model effect for whether or not parallelism is greater than expected (FDR < 0.05)
  # Read in data for simulation FDR
  FDR_dd<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_ALL_byRegion_10000runs_pvals.txt"),header=T)
  
  # Retain only the last 6 columns for 2 way FDRs
  FDR_bin<-FDR_dd[,(length(colnames(FDR_dd))-5):(length(colnames(FDR_dd)))]
  
  # Convert to binomial
  FDR_bin[FDR_bin > 0.05] <- NA
  FDR_bin[FDR_bin < 0.05] <- 1
  FDR_bin[is.na(FDR_bin)] <- 0
  
  # Convert data frame to be used in models by appending each row
  clean_para_bin<-data.frame(rbindlist(lapply(1:length(colnames(FDR_bin)),function(i){
    tmp_dd<-data.frame(FDR_bin[,i])
    colnames(tmp_dd)<-"Para_Bin"
    tmp_dd$var<-rownames(FDR_bin)
    tmp_dd$comparison<-rep(colnames(FDR_bin)[i],length(rownames(tmp_dd)))
    data.table(tmp_dd)})))
  
  # Take Para_Bin column and match onto previous data
  bin_test_dd<-cbind(eff_test_dd,clean_para_bin$Para_Bin)
  colnames(bin_test_dd)[length(colnames(bin_test_dd))]<-"Para_Bin"

  
  ###############################################################################
  
  # Build model
  #1
  lmm1<-lmer(Para_Bin~Env_Eff+theta+continent+(1|var)+(1|comparison),data=bin_test_dd)
  drop1(lmm1,test="Chisq")
  #2
  lmm2<-lmer(Para_Bin~Env_Eff*theta+continent+(1|var)+(1|comparison),data=bin_test_dd)
  drop1(lmm2,test="Chisq")
  #3
  lmm3<-lmer(Para_Bin~Env_Eff+theta*continent+(1|var)+(1|comparison),data=bin_test_dd)
  drop1(lmm3,test="Chisq")
  #4
  lmm4<-lmer(Para_Bin~Env_Eff*continent+theta+(1|var)+(1|comparison),data=bin_test_dd)
  drop1(lmm4,test="Chisq")
  #5 
  lmm5<-lmer(Para_Bin~theta+continent+(1|var)+(1|comparison),data=bin_test_dd)
  drop1(lmm5,test="Chisq")
  #6 
  lmm6<-lmer(Para_Bin~continent+(1|var)+(1|comparison),data=bin_test_dd)
  drop1(lmm6,test="Chisq")
