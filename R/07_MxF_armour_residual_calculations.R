# This script takes armour data from freshwater populations and calculates regressions to estimate marine residuals for all measures

source("R/Genomic_Parallelism_Rcode_functions.R")

lib<-as.vector(c("dplyr","data.table","ggplot2","ggpubr","cowplot"))
lapply(lib,library,character.only=TRUE)

# Read in data
dd<-read.table("outputs/PCA_stats/Armour_fresh_marine.txt",header=T)
dd[dd$Country == "Canada_BC","Country"]<-"BC"

armour_vec<-c("DS1","DS2","PS","LP","HP","BAP")
rads_vec<-unique(dd$Country)

# Find the regressions

# Loop over radiations
regressions_out<-data.frame(rbindlist(lapply(1:4,function(x){
  
  tmp<-dd[dd$Country == rads_vec[x],]
  tmp_fresh<-tmp[tmp$system == "freshwater",]
  tmp_marine<-tmp[tmp$system == "marine",]
  
  # Now loop over armours
  armour_out<-lapply(1:length(armour_vec),function(y){
    
    # Get data
    regress<-data.frame(SL=tmp_fresh$SL,
                        armour=tmp_fresh[,armour_vec[y]])
    
    # Remove any individuals that have 0, mainly for pelvis
    regress<-regress[regress$armour !=0,]
    
    # Do linear regression
    lm1<-lm(armour~SL,regress)
    
    # Get data to predict
    to_predict<-data.frame(SL=tmp$SL,
                           armour=tmp[,armour_vec[y]],
                           fish=tmp$sample_name,
                           country=tmp$Country,
                           pop=tmp$Lake_code,
                           armour_type=armour_vec[y])
    
    to_predict$predicted<-predict(lm1,newdata=to_predict)
    to_predict$residuals<-to_predict$armour-to_predict$predicted
    
    return(to_predict)
  })
  
  # Formulate an output table
  output_matrix<-matrix(nrow = nrow(tmp),
                        ncol = 6+length(armour_vec))
  output_matrix[,1]<-as.character(tmp$sample_name)
  output_matrix[,2]<-as.character(tmp$Country)
  output_matrix[,3]<-as.character(tmp$Lake_code)
  output_matrix[,4]<-as.character(tmp$system)
  output_matrix[,5]<-tmp$SL
  output_matrix[,12]<-tmp$n.plate
  
  for(i in 1:length(armour_out)){
    output_matrix[,5+i]<-armour_out[[i]]$residuals
  }
  
  colnames(output_matrix)<-c("Fish","Radiation","Pop","Type","SL",paste0("res.",armour_vec),"Plate_N")
  
  return(data.frame(output_matrix))
})))

# Reorder
regressions_out<-regressions_out[order(as.character(regressions_out$Type)),]

write.table(regressions_out,
            "outputs/MxF_armour_residuals.txt",row.names = F,quote = F,sep="\t")


################################################################################
# Read in averaged projected armour and shape residuals and plot

# Read in
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
        axis.title = element_text(size=9),
        axis.text = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = "right",
        title = element_text(size=11))+
  labs(colour="Radiation",shape="Radiation")+
  xlab("All armour traits[+] (PC1 53.7%)")+
  ylab("Pelvis Length[+] & Plate N[-] (PC2 14.5%)")+
  scale_colour_manual(values = c("red2","blue2","gold2","black","forestgreen"))+
 # scale_shape_manual(values=c(15:18))+
  ggtitle("Armour")

# Plot  Armour
shape_project<-ggplot(projected,aes(x=bodyPC1,y=bodyPC2,colour=label,shape=label))+
  geom_point(size=3)+
  stat_ellipse(level=0.95,show.legend = F)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=8),
        axis.text = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = "right",
        title = element_text(size=11))+
  labs(colour="Radiation",shape="Radiation")+
  xlab("Posterior Body[+] & Head Length[-] (PC1 29.0%)")+
  ylab("Body/Head Depth[+] & Head Elongation[+] (PC2 20.6%)")+
  scale_colour_manual(values = c("red2","blue2","gold2","black","forestgreen"))+
  # scale_shape_manual(values=c(15:18))+
  ggtitle("Shape")

# Also plot Environmental PCA
# Env
dd<-read.csv("data/All_Lakes_Env_and_Pheno_v2.csv")
colnames(dd)[1]<-"RADIATION"
env_vec<-c("Ca","Gyro","Na","pH","Schisto","Zn")
env_pca<-prcomp(dd[,env_vec],scale. = T,center = T)
#env_pca$x[,1]<--1*env_pca$x[,1]
env_plot<-plot_pca(env_pca,env_pca$x)
env_plot<-env_plot+ggtitle("Environment")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = "none",
        title = element_text(size=11))+
  scale_shape_manual(values = c(16,17,15,7))+
  #xlab("pH (-) & Ca (-) (PC1 35.3%)")+
  xlab("pH[+] & Ca[+] (PC1 35.3%)")+
  ylab("Zn[+] & Ca[-] (PC2 26.1%)")

# Also plot Gills with marine data
gill_dd<-read.csv("data/phenotype_data/gill_rakers_allpops.csv")
gill_dd$Location<-as.character(gill_dd$Location)
gill_dd$Habitat<-"Freshwater"
gill_dd[gill_dd$Pop %in% c("NYPS","MUD","LICA"),"Habitat"]<-"Marine"  
gill_dd[gill_dd$Habitat == "Marine","Location"]<-"Marine"

# Get Pop avgs to plot
gill_avg<-data.frame(gill_dd %>% group_by(Pop) %>% summarise(Habitat=mean(Habitat),
                                                  Rad=mean(Location),
                                                  N=mean(GillRakerN),
                                                  L=mean(GillRakerL)))

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
        axis.title = element_text(size=11),
        axis.text = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = "right",
        title = element_text(size=11))+
  labs(colour="Radiation",shape="Radiation",
       x=paste0("Gill Raker N"),y=paste0("Gill Raker L"))+
  scale_colour_manual(values = c("red2","blue2","gold2","black","forestgreen"))+
  #scale_shape_manual(values=c(15:18))+
  ggtitle("Trophic")

# pdf("figs/temporary_MxF_projections.pdf",width=16)
# ggarrange(armour_project,shape_project,nrow=1,ncol = 2,common.legend = T,legend = "right")
# dev.off()

# Get legend
common_leg<-get_legend(armour_project+
                         theme(legend.position="top"))

# All together
# arrange the three plots in a single row
prow <- plot_grid(
  env_plot + theme(legend.position="none"),
  armour_project + theme(legend.position="none"),
  shape_project + theme(legend.position="none"),
  gill_plot + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 2,
  ncol=2
)

# Add the legend
PC_plots<-plot_grid(common_leg, prow, ncol = 1, rel_heights = c(.1, 1))

pdf("figs/Figure4b_All_PCAs_and_projections.pdf",width=7,height=8)
PC_plots
dev.off()
