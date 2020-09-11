############################
# Pheno & Env distance matrices

# Get packages
lib<-c("data.table","ggplot2","dplyr","parallel","viridis","ggpubr","cowplot")
lapply(lib,library,character.only=T)

# Read in data
dd<-read.csv("data/All_Lakes_Env_and_Pheno_v2.csv",header=T)
env_var<-c("Ca","Na","Zn","pH","Gyro","Schisto","lake_area_km")
shape_var<-c("SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3")
armour_var<-c("Avg.resDS1","Avg.resDS2","Avg.resPS","Avg.resLP","Avg.resHP","Avg.resBAP","Avg.n.plate")
gill_var<-c("Gill.Raker.N","Resid.Raker.L")

##########################################
# Matrix of env similarity
env_PC<-prcomp(dd[,env_var],scale. = T,center = T)
env_loadings<-(env_PC$sdev^2)/sum(env_PC$sdev^2)
env_scores_trans<-env_PC$x
for(i in 1:ncol(env_scores_trans)){
  env_scores_trans[,i]<-env_scores_trans[,i]*env_loadings[i]
}
env_mat<-as.matrix(dist(env_scores_trans,method = "euclidean"))
rownames(env_mat)<-dd$LAKE
colnames(env_mat)<-dd$LAKE

# Plot
plot_env<-reshape2::melt(env_mat)
plot_env$similarity<-max(plot_env$value)-plot_env$value
env_fig<-ggplot(plot_env,aes(x=Var1,y=Var2,fill=similarity))+
  geom_tile()+
  scale_fill_viridis()+
  theme(axis.text.x=element_text(angle=90,hjust=1),
        axis.title=element_blank(),
        legend.position = "none")+
  labs(fill="Environmental Similarity")

# Plot
pdf("figs/fresh_env_similartiy_matrix.pdf",width=12,height=11)
env_fig
dev.off()

# Write to table
write.table(env_mat,"outputs/Freshwater_env_distance_mat.txt",row.names = F,quote = F,sep="\t")

##########################################
# Matrix of Pheno Similarity
pheno_PC<-prcomp(dd[,c(shape_var,armour_var,gill_var)],scale. = T,center = T)
pheno_loadings<-(pheno_PC$sdev^2)/sum(pheno_PC$sdev^2)
pheno_scores_trans<-pheno_PC$x
for(i in 1:ncol(pheno_scores_trans)){
  pheno_scores_trans[,i]<-pheno_scores_trans[,i]*pheno_loadings[i]
}
pheno_mat<-as.matrix(dist(pheno_scores_trans,method="euclidean"))
rownames(pheno_mat)<-dd$LAKE
colnames(pheno_mat)<-dd$LAKE

# Plot
plot_pheno<-melt(pheno_mat)
plot_pheno$similarity<-max(plot_pheno$value)-plot_pheno$value

pheno_fig<-ggplot(plot_pheno,aes(x=Var1,y=Var2,fill=similarity))+
  geom_tile()+
  scale_fill_viridis()+
  theme(axis.text.x=element_text(angle=90,hjust=1),
        axis.title=element_blank(),
        legend.position = "none")+
  labs(fill="Phenotypic Similarity")

# Plot
pdf("figs/fresh_pheno_similartiy_matrix.pdf",width=12,height=11)
pheno_fig
dev.off()

# Write to table
write.table(pheno_mat,"outputs/Freshwater_pheno_distance_mat.txt",row.names = F,quote = F,sep="\t")


##########################################
# Vars that show parallelism
parallel_vars<-c("Ca","pH","Avg.resDS2","Avg.resPS","Avg.resLP","Avg.n.plate","Gill.Raker.N")

parallel_PC<-prcomp(dd[,parallel_vars],scale. = T,center = T)
parallel_loadings<-(parallel_PC$sdev^2)/sum(parallel_PC$sdev^2)
parallel_scores_trans<-parallel_PC$x
for(i in 1:ncol(parallel_scores_trans)){
  parallel_scores_trans[,i]<-parallel_scores_trans[,i]*parallel_loadings[i]
}

parallel_mat<-as.matrix(dist(parallel_scores_trans,method="euclidean"))
rownames(parallel_mat)<-dd$LAKE
colnames(parallel_mat)<-dd$LAKE

# Plot
plot_parallel<-reshape2::melt(parallel_mat)
plot_parallel$similarity<-max(plot_parallel$value)-plot_parallel$value

parallel_fig<-ggplot(plot_parallel,aes(x=Var1,y=Var2,fill=similarity))+
  geom_tile()+
  scale_fill_viridis()+
  theme(axis.text.x=element_text(angle=90,hjust=1),
        axis.title=element_blank(),
        legend.position = "none")+
  labs(fill="Similarity")

# Write to table
write.table(parallel_mat,"outputs/Freshwater_parallel_distance_mat.txt",row.names = F,quote = F,sep="\t")

########################################### 
# Read in Genetic Distance matrix
FST_mat<-read.csv("outputs/All_Radiations_All_Pops_Fst.csv")
FST_mat<-as.matrix(FST_mat[,2:ncol(FST_mat)])

# Change TROU to TROUT
cols<-colnames(FST_mat)
cols[cols == "TROS"]<-"TROF"
colnames(FST_mat)<-cols
rownames(FST_mat)<-cols

# Reorder
FST_mat<-FST_mat[as.character(dd$LAKE),as.character(dd$LAKE)]

# Plot
plot_fst<-reshape2::melt(FST_mat)
plot_fst$similarity<-max(plot_fst$value)-plot_fst$value

fst_fig<-ggplot(plot_fst,aes(x=Var1,y=Var2,fill=similarity))+
  geom_tile()+
  scale_fill_viridis()+
  theme(axis.text.x=element_text(angle=90,hjust=1),
        axis.title=element_blank(),
        legend.position = "none")+
  labs(fill="Similarity")


# Plot
pdf("figs/fresh_fst_distance_similartiy_matrix.pdf",width=12,height=11)
fst_fig
dev.off()

########################################### 
# Read back in MxF parallelism matrix
MxF_mat<-as.matrix(read.table("outputs/MxF_FST_overlap_props.txt",header=T))
rownames(MxF_mat)<-colnames(MxF_mat)
MxF_melt<-reshape2::melt(MxF_mat)

# Plot for later
mxf_fig<-ggplot(MxF_melt,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  scale_fill_viridis()+
  theme(axis.text.x=element_text(angle=90,hjust=1,size=14),
        axis.text.y=element_text(size=14),
        axis.title=element_blank(),
        legend.position = "none")+
  labs(fill="Similarity")


# Rearrange and remove NAs
colnames(MxF_melt)<-c("pop1","pop2","MxF_overlap")

# Change TROU to TROUT
MxF_melt$pop1<-as.character(MxF_melt$pop1)
MxF_melt$pop2<-as.character(MxF_melt$pop2)

MxF_melt[MxF_melt$pop1 == "TROS","pop1"]<-"TROF"
MxF_melt[MxF_melt$pop2 == "TROS","pop2"]<-"TROF"

MxF_melt<-na.omit(MxF_melt)
MxF_melt$comparison<-paste0(MxF_melt$pop1,"-",MxF_melt$pop2)

# Now add the similarity for each one
MxF_melt$env_distance<-NA
MxF_melt$pheno_distance<-NA
MxF_melt$parallel_distance<-NA
MxF_melt$genetic_distance<-NA

for(i in 1:nrow(MxF_melt)){
  print(i)
  # Get env
  MxF_melt$env_distance[i]<-plot_env[as.character(plot_env$Var1) == as.character(MxF_melt$pop1[i]) &
                                                    as.character(plot_env$Var2) == as.character(MxF_melt$pop2[i]),"value"]
  # Get pheno
  MxF_melt$pheno_distance[i]<-plot_pheno[as.character(plot_pheno$Var1) == as.character(MxF_melt$pop1[i]) &
                                       as.character(plot_pheno$Var2) == as.character(MxF_melt$pop2[i]),"value"]
  # Get parralel
  MxF_melt$parallel_distance[i]<-plot_parallel[as.character(plot_parallel$Var1) == as.character(MxF_melt$pop1[i]) &
                                       as.character(plot_parallel$Var2) == as.character(MxF_melt$pop2[i]),"value"]
  # Get parallel
  MxF_melt$genetic_distance[i]<-plot_fst[as.character(plot_fst$Var1) == as.character(MxF_melt$pop1[i]) &
                                       as.character(plot_fst$Var2) == as.character(MxF_melt$pop2[i]),"value"]
}

# Remove duplicate comparisons
MxF_unique<-MxF_melt[!(duplicated(MxF_melt$pheno_distance)),]

# Add radiations
MxF_unique$rad1<-NA
MxF_unique$rad2<-NA
for(i in 1:nrow(MxF_unique)){
  MxF_unique$rad1[i]<-as.character(dd[dd$LAKE==MxF_unique$pop1[i],"RADIATION"])
  MxF_unique$rad2[i]<-as.character(dd[dd$LAKE==MxF_unique$pop2[i],"RADIATION"])
}

# Add variable for intra/inter radiation comparisons
MxF_unique$radiation_level<-"Across"
MxF_unique[MxF_unique$rad1 == MxF_unique$rad2,"radiation_level"]<-"Within"

# Add variable for continents
MxF_unique$continent_level<-"Across"
MxF_unique[MxF_unique$rad1 == "Alaska" &
             MxF_unique$rad2 == "BC" , "continent_level"]<-"Within"
MxF_unique[MxF_unique$rad1 == "BC" &
             MxF_unique$rad2 == "Alaska" , "continent_level"]<-"Within"
MxF_unique[MxF_unique$rad1 == "Scotland" &
             MxF_unique$rad2 == "Iceland" , "continent_level"]<-"Within"
MxF_unique[MxF_unique$rad1 == "Iceland" &
             MxF_unique$rad2 == "Scotland" , "continent_level"]<-"Within"
MxF_unique[MxF_unique$radiation_level == "Within","continent_level"]<-"Within"

# Plot histogram of Fst overlap
Fst_hist<-ggplot(MxF_unique,aes(x=MxF_overlap))+
  geom_histogram()+
  geom_vline(xintercept = 0.05,linetype="dashed")+
  labs(y="Count",x=expression(F[ST]~Outlier~Overlap))+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))

Fst_hist2<-ggplot(MxF_unique,aes(x=MxF_overlap))+
  geom_histogram()+
  geom_vline(xintercept = 0.05,linetype="dashed")+
  labs(y="Count",x=expression(F[ST]~Outlier~Overlap))+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text = element_text(size=18))+
  facet_wrap(~radiation_level)

# Plot together
pdf("figs/Fst_overlap_histograms.pdf",width=12,height=12)
ggarrange(Fst_hist,Fst_hist2,ncol=1,nrow=2,labels = "auto")
dev.off()

# Plot effect of env and mxf
env_effect<-ggplot(MxF_unique,aes(x=env_distance,y=log(MxF_overlap),colour=radiation_level))+
  geom_point(alpha=0.15)+
  geom_smooth(method="lm",se = F)+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  labs(colour="Radiation\nLevel",x="Environmental Distance",y=expression(MxF~F[ST]~Overlap~(log)))+
  scale_colour_manual(values=c("Across"="royalblue4",
                               "Within"="red2"))

pdf("figs/Env_effect_on_fst_overlap.pdf",width=10,height=8)
env_effect
dev.off()

# Plot effect of pheno and mxf
pheno_effect<-ggplot(MxF_unique,aes(x=pheno_distance,y=log(MxF_overlap),colour=radiation_level))+
  geom_point(alpha=0.15)+
  geom_smooth(method="lm",se = F)+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  labs(colour="Radiation\nLevel",x="Phenotypic Distance",y=expression(MxF~F[ST]~Overlap~(log)))+
  scale_colour_manual(values=c("Across"="royalblue4",
                               "Within"="red2"))

# Plot genetic effect of pheno and mxf
fst_effect<-ggplot(MxF_unique,aes(x=genetic_distance,y=log(MxF_overlap),colour=radiation_level))+
  geom_point(alpha=0.15)+
  geom_smooth(method="lm",se = F)+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  labs(colour="Radiation\nLevel",x="Genetic Distance",y=expression(MxF~F[ST]~Overlap~(log)))+
  scale_colour_manual(values=c("Across"="royalblue4",
                               "Within"="red2"))

pdf("figs/Pheno_effect_on_fst_overlap.pdf",width=10,height=8)
pheno_effect
dev.off()


# Save analysis file
write.table(MxF_unique,"outputs/Parallelism_distance_analysis_data.txt",
            row.names = F,quote = F,sep = "\t")

######################
# Final Matrix Plot
col1<-plot_grid(mxf_fig+ggtitle("a")+theme(plot.title = element_text(size = 40, face = "bold")))
col2<-plot_grid(env_fig+theme(axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks = element_blank(),
                              axis.title = element_blank(),
                              plot.title = element_text(size = 40, face = "bold"))+ggtitle("b"),
                pheno_fig+theme(axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks = element_blank(),
                              axis.title = element_blank(),
                              plot.title = element_text(size = 40, face = "bold"))+ggtitle("c"),
                fst_fig+theme(axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks = element_blank(),
                              axis.title = element_blank(),
                              plot.title = element_text(size = 40, face = "bold"))+ggtitle("d"),
                align = "hv",
                ncol=1)
col3<-plot_grid(env_effect+theme(legend.position = "none",plot.title = element_text(size = 40, face = "bold"))+ggtitle("e"),
                pheno_effect+theme(legend.position = "none",plot.title = element_text(size = 40, face = "bold"))+ggtitle("f"),
                fst_effect+theme(legend.position = "none",plot.title = element_text(size = 40, face = "bold"))+ggtitle("g"),
                ncol=1,
                align="hv",
                axis="tblr")

common_leg<-get_legend(fst_effect+theme(legend.position="right",
                                        legend.title = element_text(size=24),
                                        legend.text = element_text(size=22)))


all_figs<-plot_grid(col1,col2,col3,
                    ncol=3,align="hv",
                    rel_widths = c(3,1,1.5,0.3),axis = "tblr")

all_figs_with_legend<-plot_grid(all_figs, common_leg, ncol = 2, rel_widths = c(9, 1))


pdf("figs/Figure4_all_matrix_figures.pdf",width=24,height=12)
all_figs_with_legend
dev.off()
