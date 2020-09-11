########################################################
# Do PCs within radiations correlate with env variables
lib<-c("ggplot2","dplyr","data.table","corrplot","Hmisc")
lapply(lib,library,character.only=T)

rads<-c("Alaska","BC","Iceland","Scotland")
short<-c("Al","BC","Ic","Sc")

# Read in PC res
PC_data<-lapply(1:length(rads),function(x){
  tmp<-read.table(paste0("outputs/",rads[x],"_plink_out_NoLD_PCA.eigenvec"))
  colnames(tmp)<-c("ID","ID2",paste0("PC",1:(ncol(tmp)-2)))
  tmp$pop<-substr(tmp$ID, 1, 4)
  return(tmp)
})

# Only retain vals > 10%
eigenvals<-lapply(1:4,function(x){
  tmp<-read.table(paste0("outputs/",rads[x],"_plink_out_NoLD_PCA.eigenval"))
  tmp$prop<-tmp$V1/sum(tmp$V1)
  return(tmp[1:2,])
})

# Trim PCs to keep major eigenvectors and average
PC_avg<-lapply(1:length(PC_data),function(x){
  
  tmp<-PC_data[[x]]
  
  tmp<-tmp[,c(3:(3+nrow(eigenvals[[x]])-1),ncol(tmp))]
  tmp2<-data.frame(
    tmp %>%
      group_by(pop) %>%
      summarise_all(mean))
})

# Now bind with the relevant environmental data
env_data<-read.csv("data/All_Lakes_Env_and_Pheno_v2.csv")
env_data$LAKE<-as.character(env_data$LAKE)
env_data[env_data$LAKE == "TROF","LAKE"]<-"TROS"

PCs_envs<-lapply(1:4,function(x){
  PCs<-PC_avg[[x]]
  
  envs<-data.frame(rbindlist(lapply(PCs$pop,function(pop){
  env_data[env_data$LAKE == pop,
                 c("Ca","Na","pH","Zn","Gyro","Schisto","lake_area_km",
                   "SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
                   "Avg.resDS1","Avg.resDS2","Avg.resPS","Avg.resLP","Avg.resHP","Avg.resBAP","Avg.n.plate",
                   "Gill.Raker.N","Resid.Raker.L")]
  })))
  tmp<-cbind(PCs,envs)
  
  # Correlation matrix
  cors<-rcorr(as.matrix(tmp[,2:ncol(tmp)]),type="spearman")
  cors.r<-cors$r[1:2,3:ncol(cors$r)]
  cors.r<-data.frame(cors.r,
                   PC=paste0(rads[x]," SNP ",rownames(cors.r)))
  cors.p<-cors$P[1:2,3:ncol(cors$P)]
  cors.p<-data.frame(cors.p,
                   PC=paste0(rads[x]," SNP ",rownames(cors.p)))
  return(list(cors.r,cors.p))
  
})

cors_out<-rbind(PCs_envs[[1]][[1]],PCs_envs[[2]][[1]],PCs_envs[[3]][[1]],PCs_envs[[4]][[1]])
cor_mat<-as.matrix(cors_out[,1:(ncol(cors_out)-1)])
rownames(cor_mat)<-cors_out$PC

p_out<-rbind(PCs_envs[[1]][[2]],PCs_envs[[2]][[2]],PCs_envs[[3]][[2]],PCs_envs[[4]][[2]])
p_mat<-as.matrix(p_out[,1:(ncol(p_out)-1)])
rownames(p_mat)<-p_out$PC

# Plot
pdf("figs/Pop_structure_variable_correlations.pdf",width=12,height=5)
corrplot(cor_mat,
         method="color",
         p.mat = p_mat,
         sig.level=0.05,
         insig="blank",
         addCoef.col = "black",
         tl.srt=30)
dev.off()
