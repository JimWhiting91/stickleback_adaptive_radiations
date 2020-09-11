#################
# This script reads in env and pheno data to explore multivariate associations and produces figure 2.
#################

# Stuart Functions
source("R/Genomic_Parallelism_Rcode_functions.R")

# Load packages
lib<-c("ggbiplot","cowplot","ggpubr","ggrepel","data.table","VennDiagram","ggplot2","parallel","dplyr","tidyverse")
lapply(lib,library,character.only=TRUE)

# Read in data
dd<-read.csv("data/All_Lakes_Env_and_Pheno_v2.csv")
colnames(dd)[1]<-"RADIATION"

##### PCAs #####
# Env
env_vec<-c("Ca","Gyro","Na","pH","Schisto","Zn","lake_area_km")
env_pca<-prcomp(dd[,env_vec],scale. = T,center = T)
env_pca
summary(env_pca)

env_plot<-plot_pca(env_pca,env_pca$x)
env_plot<-env_plot+ggtitle("Environment",
                           subtitle = "PC1 = pH, Ca, Gyro (+)\nPC2 = Zn, Schisto (+) & Ca (-)")+
  theme(title = element_text(size=18),
        legend.position = "none")+
  xlab("PC1 31.4%")+
  ylab("PC2 22.4%")

# Write PC matrix to output
write.table(env_pca$rotation,"outputs/PCA_stats/PCA_loading_matrix_Environment.txt",
            row.names=F,quote = F,sep="\t")

# Average by Lake and export
Lake_env<-data.frame(LAKE=dd$LAKE,
                    env_pca$x)
Lake_env2<-data.frame(Lake_env %>%
                           group_by(LAKE) %>%
                           summarise(PC1=mean(na.omit(PC1)),
                                     PC2=mean(na.omit(PC2)),
                                     PC3=mean(na.omit(PC3)),
                                     PC4=mean(na.omit(PC4))))

Lake_env2<-Lake_env2[match(dd$LAKE, Lake_env2$LAKE),]
Lake_env2$RAD<-dd$RADIATION

write.table(Lake_env2,
            "outputs/PCA_stats/PCA_scores_lake_avg_Env.txt",
            row.names=F,quote = F,sep="\t")

# Also Plot only PC Loading axes
env_loading<-ggbiplot(env_pca,alpha = 0)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=16))

pdf("figs/Supp_Env_loading_axes.pdf")
env_plot
env_loading
dev.off()

# Armour (Residuals calculated in separate, PCA is done on all individuals and avg plotted)
# Read in pheno data
pheno_dd<-read.csv("data/All_Shape_and_Pheno_data.csv",header=T)
armour_vec<-colnames(pheno_dd)[c(17:23)]
armour_pca<-prcomp(na.omit(pheno_dd[,armour_vec]),scale. = T,center = T)
armour_pca$x[,1]<--1*armour_pca$x[,1]
armour_pca$x[,2]<--1*armour_pca$x[,2]
armour_pca_scores<-data.frame(data.frame(lake=as.character(pheno_dd$Lake_code)[1:1301],
                              PC1=armour_pca$x[,1],
                              PC2=armour_pca$x[,2]) %>%
group_by(lake) %>%
summarise_each(funs(mean)))

# Reorder
armour_pca_scores<-armour_pca_scores[match(dd$LAKE, armour_pca_scores$lake),]

# Plot
armour_plot<-plot_pca(armour_pca,armour_pca_scores[,2:3])
armour_plot<-armour_plot+ggtitle("Armour",subtitle="PC1 = All armour traits (+)\nPC2 = Pelvis Length (+) & Plate N (-)")+
  theme(title = element_text(size=18))+
              xlab("PC1 53.7%")+
              ylab("PC2 14.5%")

# Write PC matrix to output
write.table(armour_pca$rotation,"outputs/PCA_stats/PCA_loading_matrix_Armour.txt",
            row.names=F,quote = F,sep="\t")
write.table(data.frame(cbind(as.character(pheno_dd$sample_name),armour_pca$x)),
            "outputs/PCA_stats/PCA_scores_individuals_Armour.txt",
            row.names=F,quote = F,sep="\t")

# Average by Lake and export
Lake_armour<-data.frame(LAKE=pheno_dd$Lake_code,
                     armour_pca$x)
Lake_armour2<-data.frame(Lake_armour %>%
                           group_by(LAKE) %>%
                           summarise(PC1=mean(na.omit(PC1)),
                                     PC2=mean(na.omit(PC2)),
                                     PC3=mean(na.omit(PC3)),
                                     PC4=mean(na.omit(PC4))))

Lake_armour2<-Lake_armour2[match(dd$LAKE, Lake_armour2$LAKE),]
Lake_armour2$RAD<-dd$RADIATION

# Write
write.table(Lake_armour2,
            "outputs/PCA_stats/PCA_scores_lake_avg_Armour.txt",
            row.names=F,quote = F,sep="\t")

# Shape
# PCA is calculated externally, note scale = F here because all same unit. Again calculated over individuals with pop avg plotted
shape_plot<-ggplot(dd,aes(x=SHAPE.IND.PC1,y=SHAPE.IND.PC2,colour=RADIATION,shape=RADIATION))+
  geom_point(size=3)+
  stat_ellipse(level=0.95,show.legend = F)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(size=16),
        legend.text = element_text(size=16),
        legend.position = "right",
        title = element_text(size=18))+
  labs(colour="Radiation",shape="Radiation",
       x=paste0("PC1 29.0%"),
       y=paste0("PC2 20.6%"))+
  scale_colour_manual(values = c("red2","blue2","gold2","forestgreen"))+
  scale_shape_manual(values=c(15:18))+
  ggtitle("Shape",
          subtitle="PC1 = Posterior Body (+) & Head Length (-)\nPC2 = Body/Head Depth (+) & Head Elongation (+)")

# Shape Lake Avg
Shape_Lake_avgs<-dd[,c("LAKE","RADIATION","SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3")]
write.table(Shape_Lake_avgs,
            "outputs/PCA_stats/PCA_scores_lake_avg_Shape.txt",
            row.names=F,quote = F,sep="\t")

# Gill Rakers
gill_plot<-ggplot(dd,aes(x=Gill.Raker.N,y=Resid.Raker.L,colour=RADIATION,shape=RADIATION))+
  geom_point(size=3)+
  stat_ellipse(level=0.95,show.legend = F)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(size=16),
        legend.text = element_text(size=16),
        legend.position = "right",
        title = element_text(size=18))+
  labs(colour="Radiation",shape="Radiation",
       x=paste0("Gill Raker N"),y=paste0("Gill Raker L"))+
  scale_colour_manual(values = c("red2","blue2","gold2","forestgreen"))+
  scale_shape_manual(values=c(15:18))+
  ggtitle("Trophic")

# All together
PC_plots<-ggarrange(env_plot,armour_plot,shape_plot,gill_plot,
                    ncol=2,nrow=2,common.legend = T)


##### Also plot overlap of PCs as histogram ########
all_scores<-data.frame(scores=c(dd$Env.PC1,dd$Env.PC2,
                                dd$ARMOUR.IND.WITHPLATE.PC1,dd$ARMOUR.IND.WITHPLATE.PC2,
                                dd$SHAPE.IND.PC1,dd$SHAPE.IND.PC2,
                                dd$Gill.Raker.N,dd$Resid.Raker.L),
                       PC=c(rep(rep(c("PC1","PC2"),each = nrow(dd)),3),rep(c("N","Length"),each=nrow(dd))),
                       var=rep(c("Environment","Armour","Shape","Gill Raker"),each=nrow(dd)*2),
                       rad=rep(dd$RADIATION,8))

# Set factors
all_scores$var_F<-factor(all_scores$var,levels=c("Environment","Armour","Shape","Gill Raker"))

# Change rad labels
all_scores[all_scores$rad == "ALASKA","rad"] <- "Alaska"
all_scores[all_scores$rad == "BC","rad"] <- "BC"
all_scores[all_scores$rad == "ICELAND","rad"] <- "Iceland"
all_scores[all_scores$rad == "SCOTLAND","rad"] <- "Scotland"

score_hists<-
  ggplot(all_scores[all_scores$var != "Gill",],aes(x=scores,fill=rad))+
  geom_density(alpha=0.3)+
  scale_fill_manual(values = c("red2","blue2","gold2","forestgreen"))+
  facet_wrap(var_F~PC,scales = "free",ncol=2)+
  labs(y="Density",fill="Radiation")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(size=20),
        legend.position = "top",
        legend.text = element_text(size=20),
        legend.title = element_blank())

# pdf("figs/FigureSX_PC_overlap_histograms.pdf",height=10,width=6)
# score_hists
# dev.off()


###################################################
# And Plot with original PCs and Vectors - May need to run script 10 first...
vector_figs<-readRDS("figs/All_vector_figures.RDS")
# 
# pdf("figs/Figure2_NatEvoEco_PCs_Vectors.pdf",height=10,width=32)
# plot_grid(PC_plots,
#           score_hists+theme(strip.text = element_text(size=16),
#                             axis.text = element_text(size=15)),
#           vector_figs[[1]],
#           vector_figs[[2]],
#           ncol=4,rel_widths=c(1.7,1,0.7,0.7),
#           labels = "AUTO",
#           label_size = 24,align = "h",axis="tblr")
# dev.off()

####################################################
# Also plot with marine individuals projected...
projected<-read.table("outputs/MxF/AveragePCs_with_projected_marine.txt",header=T)
colnames(projected)[c(2,4)]<-c("POP","bodyPC2")

marines<-c("LICA",
           "MUDL",
           "NYPS",
           "OBSM")

# Add Label column
projected$label<-as.character(projected$Radiation)
projected[as.character(projected$POP) %in% marines,"label"]<-"Marine"
projected$type<-"Freshwater"
projected[as.character(projected$POP) %in% marines,"type"]<-"Marine"

# Plot  Armour
armour_project<-ggplot(projected,aes(x=armourPC1,y=armourPC2,colour=label,shape=label))+
  geom_point(size=3)+
  stat_ellipse(level=0.95,show.legend = F)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(size=16),
        legend.text = element_text(size=16),
        legend.position = "right",
        title = element_text(size=18))+
  labs(colour="Radiation",shape="Radiation")+
  scale_colour_manual(values = c("red2","blue2","gold2","black","forestgreen"))+
  ggtitle("Armour",subtitle="PC1 = All armour traits (+)\nPC2 = Pelvis Length (+) & Plate N (-)")+
  theme(title = element_text(size=18))+
  xlab("PC1 53.7%")+
  ylab("PC2 14.5%")

# Plot  Armour
shape_project<-ggplot(projected,aes(x=bodyPC1,y=bodyPC2,colour=label,shape=label))+
  geom_point(size=3)+
  stat_ellipse(level=0.95,show.legend = F)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(size=16),
        legend.text = element_text(size=16),
        legend.position = "right",
        title = element_text(size=18))+
  labs(colour="Radiation",shape="Radiation",
       x=paste0("PC1 29.0%"),
       y=paste0("PC2 20.6%"))+
  scale_colour_manual(values = c("red2","blue2","gold2","black","forestgreen"))+
  #scale_shape_manual(values=c(15:18))+
  ggtitle("Shape",
          subtitle="PC1 = Posterior Body (+) & Head Length (-)\nPC2 = Body/Head Depth (+) & Head Elongation (+)")

# Also plot Gills with marine data
gill_dd<-read.csv("data/phenotype_data/gill_rakers_allpops.csv")
gill_dd$Location<-as.character(gill_dd$Location)
gill_dd$Habitat<-"Freshwater"
gill_dd[gill_dd$Pop %in% c("NYPS","MUD","LICA"),"Habitat"]<-"Marine"  
gill_dd[gill_dd$Habitat == "Marine","Location"]<-"Marine"

# Get Pop avgs to plot
gill_avg<-data.frame(gill_dd %>% group_by(Pop) %>% dplyr::summarise(Habitat=mean(Habitat),
                                                                    Rad=mean(Location),
                                                                    N=mean(GillRakerN),
                                                                    L=mean(GillRakerL)))

# Fill gaps
for(i in 1:nrow(gill_avg)){
  gill_avg$Habitat[i]<-gill_dd[gill_dd$Pop == gill_avg$Pop[i],"Habitat"][1]
  gill_avg$Rad[i]<-gill_dd[gill_dd$Pop == gill_avg$Pop[i],"Location"][1]
}

gill_avg<-gill_avg[2:nrow(gill_avg),]

# Plot
gill_plot<-ggplot(gill_avg,aes(x=N,y=L,colour=Rad,shape=Rad))+
  geom_point(size=3)+
  stat_ellipse(level=0.95,show.legend = F)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(size=16),
        legend.text = element_text(size=20),
        legend.title = element_text(size=22),
        legend.position = "top",
        title = element_text(size=18))+
  labs(colour="Radiation",shape="Radiation",
       x=paste0("Gill Raker N"),y=paste0("Gill Raker L"))+
  scale_colour_manual(values = c("red2","blue2","gold2","black","forestgreen"))+
  #scale_shape_manual(values=c(15:18))+
  ggtitle("Trophic")

# Get legend
all_plot_legend<-get_legend(gill_plot + theme(legend.box.margin = margin(0, 0, 0, 0)))

# New all togetherther
PC_plots<-ggarrange(env_plot,armour_project,shape_project,gill_plot,
                    ncol=2,nrow=2,legend = "none")
# Add legend
PC_plots2<-plot_grid(all_plot_legend,
                     PC_plots,ncol=1,nrow=2,rel_heights=c(1,20))

pdf("figs/Figure2_NatEvoEco_PCs_Vectors.pdf",height=10,width=32)
plot_grid(PC_plots2,
          score_hists+theme(strip.text = element_text(size=16),
                            axis.text = element_text(size=15)),
          vector_figs[[1]],
          vector_figs[[2]],
          ncol=4,rel_widths=c(1.7,1,0.7,0.7),
          labels = "auto",
          label_size = 30,align = "h",axis="tblr")
dev.off()


