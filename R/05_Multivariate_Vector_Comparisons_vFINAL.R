#################
# This script examines parallelism of phenotypic and environmental multivariate variation across radiation comparisons using
# methods similar to Stuart et al 2017 (Nat EcoEvo)

# PCs scaled by PC loading
#################

# Stuart Functions
source("R/Stuart_Nature_Rcode_functions.R")
source("R/Genomic_Parallelism_Rcode_functions.R")

# Load packages
lib<-as.vector(c("ggbiplot","pbapply","ggpubr","ggrepel","data.table","VennDiagram","ggplot2","parallel","dplyr"))
lapply(lib,library,character.only=TRUE)

# Read in data
dd<-read.csv("data/All_Lakes_Env_and_Pheno_v2.csv")
colnames(dd)[1]<-"RADIATION"

##### PCAs #####
# Env
env_vec<-c("Ca","Gyro","Na","pH","Schisto","Zn")
env_pca<-prcomp(dd[,env_vec],scale. = T,center = T)
#env_pca$x[,1]<--1*env_pca$x[,1]
env_plot<-plot_pca(env_pca,env_pca$x)
env_plot<-env_plot+ggtitle("Environment")+
  theme(title = element_text(size=18),
        legend.position = "top")+
  #xlab("pH (-) & Ca (-) (PC1 35.3%)")+
  xlab("pH (+) & Ca (+) (PC1 35.3%)")+
  ylab("Zn (+) & Ca (-) (PC2 26.1%)")

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
armour_plot<-armour_plot+ggtitle("Armour")+theme(title = element_text(size=18))+
              xlab("All armour traits (+) (PC1 53.7%)")+
              ylab("Pelvis Length (+) & Plate N (-) (PC2 14.5%)")

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
       x=paste0("Posterior Body (+) & Head Length (-) (PC1 29.0%)"),
       y=paste0("Body/Head Depth (+) & Head Elongation (+) (PC2 20.6%)"))+
  scale_colour_manual(values = c("red2","blue2","gold2","forestgreen"))+
  scale_shape_manual(values=c(15:18))+
  ggtitle("Shape")

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
PC_plots<-ggarrange(env_plot,shape_plot,armour_plot,gill_plot,
                    ncol=2,nrow=2,common.legend = T)
pdf("figs/Figure4_All_PCAs.pdf",width=12,height=12)
PC_plots
dev.off()

##### Calculate the degree of change for each variable, for each radiation ######

# Input scalings as PC weighting
shape_scales<-c(29.0,20.6,16.5)
env_scales<-c(35.3, 26.1,18.6,11.4)
armour_scales<-c(53.7,14.5,11.0,10.3)


# Perform scaling
# Env
dd$Env.PC1_scaled<-dd$Env.PC1*env_scales[1]
dd$Env.PC2_scaled<-dd$Env.PC2*env_scales[2]
dd$Env.PC3_scaled<-dd$Env.PC3*env_scales[3]
dd$Env.PC4_scaled<-dd$Env.PC4*env_scales[4]
# Armour
dd$ARMOUR.IND.WITHPLATE.PC1_scaled<-dd$ARMOUR.IND.WITHPLATE.PC1*armour_scales[1]
dd$ARMOUR.IND.WITHPLATE.PC2_scaled<-dd$ARMOUR.IND.WITHPLATE.PC2*armour_scales[2]
dd$ARMOUR.IND.WITHPLATE.PC3_scaled<-dd$ARMOUR.IND.WITHPLATE.PC3*armour_scales[3]
dd$ARMOUR.IND.WITHPLATE.PC4_scaled<-dd$ARMOUR.IND.WITHPLATE.PC4*armour_scales[4]
# Env
dd$SHAPE.IND.PC1_scaled<-dd$SHAPE.IND.PC1*shape_scales[1]
dd$SHAPE.IND.PC2_scaled<-dd$SHAPE.IND.PC2*shape_scales[2]
dd$SHAPE.IND.PC3_scaled<-dd$SHAPE.IND.PC3*shape_scales[3]

# Scale Gills
dd$Gill.Raker.N_scaled<-scale(dd$Gill.Raker.N)
dd$Resid.Raker.L_scaled<-scale(dd$Resid.Raker.L)

##### Plot scaled PC Plots ######
env_scale_fig<-plot_scaled(dd$Env.PC1_scaled,
                              dd$Env.PC2_scaled,
                              env_scales[1],
                              env_scales[2])

armour_scale_fig<-plot_scaled(dd$ARMOUR.IND.WITHPLATE.PC1_scaled,
                              dd$ARMOUR.IND.WITHPLATE.PC2_scaled,
                              armour_scales[1],
                              armour_scales[2])

shape_scale_fig<-plot_scaled(dd$SHAPE.IND.PC1_scaled,
                              dd$SHAPE.IND.PC2_scaled,
                              shape_scales[1],
                              shape_scales[2])

gill_scale_fig<-plot_scaled(dd$Gill.Raker.N_scaled,
                             dd$Resid.Raker.L_scaled,
                             1,
                             1)+
  xlab("Number")+ylab("Length")

pdf("figs/scaled_PC_plots.pdf",width=8,height=4)
env_scale_fig
armour_scale_fig
shape_scale_fig
dev.off()


#################################

# List out the variables

variable_vec<-c("Env.PC1_scaled","Env.PC2_scaled","Env.PC3_scaled","Env.PC4_scaled",
                "SHAPE.IND.PC1_scaled","SHAPE.IND.PC2_scaled","SHAPE.IND.PC3_scaled",
                "ARMOUR.IND.WITHPLATE.PC1_scaled","ARMOUR.IND.WITHPLATE.PC2_scaled","ARMOUR.IND.WITHPLATE.PC3_scaled","ARMOUR.IND.WITHPLATE.PC4_scaled",
                "Gill.Raker.N_scaled", "Resid.Raker.L_scaled")

group_list<-list(variable_vec[1:4],variable_vec[5:7],variable_vec[8:11],variable_vec[12:13])
group_vec<-c("ENV","SHAPE","ARMOUR","GILL")
rad_vec<-c("ALASKA","BC","ICELAND","SCOTLAND")

# Calculate all vectors of change across radiations  
change_out<-lapply(1:length(group_list),function(x){
  
  # Subset
  tmp_var<-variable_vec[variable_vec %in% group_list[[x]]]
  
  # Across tmp_var calculate vectors for all radiations
  out<-data.frame(rbindlist(lapply(1:length(tmp_var),function(y){
    A_tmp<-data.frame(vec=calculate_change(dd,tmp_var[y],"ALASKA"),
                      rad="ALASKA",
                      var=tmp_var[y])
    B_tmp<-data.frame(vec=calculate_change(dd,tmp_var[y],"BC"),
                      rad="BC",
                      var=tmp_var[y])
    I_tmp<-data.frame(vec=calculate_change(dd,tmp_var[y],"ICELAND"),
                      rad="ICELAND",
                      var=tmp_var[y])
    N_tmp<-data.frame(vec=calculate_change(dd,tmp_var[y],"SCOTLAND"),
                      rad="SCOTLAND",
                      var=tmp_var[y])
    return(rbind(A_tmp,B_tmp,I_tmp,N_tmp))
  })))
  
  # Add label
  out$group<-group_vec[x]
  
  # Return
  out
})

##### CREATE FINAL VECTORS AND OUTPUT THETA L, DELTA_L #####
angles_out<-lapply(1:4,function(x){
  
  # List out vectors
  vec_dd<-t(data.frame(lapply(1:length(rad_vec),function(y){
    return(change_out[[x]][change_out[[x]]$rad==rad_vec[y],"vec"])
  })))
  
  avg_vec_diffs<-data.frame(vec_dd)
  colnames(avg_vec_diffs)<-group_list[[x]]
  
  # Stuart calc functions
  # Use Stuart et al function to calculate theta and L
  rad.theta <- f.theta.jim(avg_vec_diffs, unique.id = rad_vec, select.col =1:length(avg_vec_diffs))
  rad.theta.radians <- rad.theta[[1]]
  rad.theta.degrees <- rad.theta[[2]]
  
  rad.theta.degrees_vec<-rad.theta.degrees[upper.tri(rad.theta.degrees, diag = FALSE)]
  
  rad.deltaL <- f.deltaL(avg_vec_diffs, unique.id = rad_vec, select.col =1:length(avg_vec_diffs[1,]))
  rad.theta.deltaL_vec<-abs(rad.deltaL[[1]][upper.tri(rad.deltaL[[1]], diag = FALSE)])
  
  # Combine the two
  out_mat<-matrix(nrow = 2,ncol = length(combn(rad_vec,2)[1,]))
  out_mat[1,]<-rad.theta.degrees_vec
  out_mat[2,]<-rad.theta.deltaL_vec
  out_dd<-data.frame(out_mat,
                     measure=c("Theta_Deg","deltaL"))
  
  # Rename columns
  colnames(out_dd)[1:6]<-c("Alaska_BC",
                           "Alaska_Iceland",
                           "BC_Iceland",
                           "Alaska_Scotland",
                           "BC_Scotland",
                           "Iceland_Scotland")
  
  # Reorder columns alphabetically
  out_dd<-out_dd[,order(colnames(out_dd))]
  
  # Give rad Info
  out_dd$group<-group_vec[x]
  
  # Also get length of vectors
  out_dd_L<-rad.deltaL[[2]]
  colnames(out_dd_L)<-c("Vector_L","Radiation")
  out_dd_L$group<-group_vec[x]
  
  # Return from CHR
  return(list(out_dd,out_dd_L,rad.theta.radians,rad.deltaL[[1]]))
  
})

# Extract the outputs
# out_dd
out_dd<-data.frame(rbindlist(lapply(1:length(angles_out),function(x){
  return(angles_out[[x]][[1]])})))

# Modify out_dd so that we only have angles less than 90
for(i in 1:6){
  tmp<-data.frame(new=180-out_dd[out_dd$measure == "Theta_Deg",i],
                  old=out_dd[out_dd$measure == "Theta_Deg",i])
  out_dd[out_dd$measure == "Theta_Deg",i]<-do.call(pmin,tmp)
}

# Transform out_dd and export for modelling
out_dd_2<-data.frame(t(out_dd))
colnames(out_dd_2)<-c("ENV_theta","ENV_deltaL",
                      "SHAPE_theta","SHAPE_deltaL",
                      "ARMOUR_theta","ARMOUR_deltaL",
                      "GILL_theta","GILL_deltaL")
out_dd_2<-out_dd_2[1:6,]
out_dd_2$rads<-rownames(out_dd_2)

# out_dd_L
out_dd_L<-data.frame(rbindlist(lapply(1:length(angles_out),function(x){
  return(angles_out[[x]][[2]])})))

# Write
write.table(out_dd_2,
            "outputs/05_Observed_theta_delta_vectors.txt",
            quote = F,row.names = F,sep="\t")

##### MANTEL TESTS OF CORRELATION BETWEEN ENV, SHAPE and ARMOUR for theta and delta L #####

# THETA
env_theta<-as.dist(angles_out[[1]][[3]],diag = F,upper=F)
shape_theta<-as.dist(angles_out[[2]][[3]],diag = F,upper=F)
armour_theta<-as.dist(angles_out[[3]][[3]],diag = F,upper=F)
gill_theta<-as.dist(angles_out[[4]][[3]],diag = F,upper=F)

all.mantel.theta.env_shape <- mantel.rtest(env_theta, shape_theta, nrepet = 9999) # Not significant
all.mantel.theta.env_armour <- mantel.rtest(env_theta, armour_theta, nrepet = 9999) # Not significant
all.mantel.theta.env_pheno <- mantel.rtest(env_theta, gill_theta, nrepet = 9999) # Not significant

# delta L
env_deltaL<-abs(as.dist(angles_out[[1]][[4]],diag = F,upper=F))
shape_deltaL<-abs(as.dist(angles_out[[2]][[4]],diag = F,upper=F))
armour_deltaL<-abs(as.dist(angles_out[[3]][[4]],diag = F,upper=F))
gill_deltaL<-abs(as.dist(angles_out[[4]][[4]],diag = F,upper=F))

all.mantel.deltaL.env_shape <- mantel.rtest(env_deltaL, shape_deltaL, nrepet = 9999) # Not significant
all.mantel.deltaL.env_armour <- mantel.rtest(env_deltaL, armour_deltaL, nrepet = 9999) # Significant: Extent of environmental variation is correlated with extent of armour variation
all.mantel.deltaL.env_gill <- mantel.rtest(env_deltaL, gill_deltaL, nrepet = 9999) # Not significant

##### Spearman's Rhos of correlations between env theta delta and phenos #####
# Theta corrs
shape_cor<-cor.test(as.vector(env_theta),as.vector(shape_theta),method = "spearman")$estimate
armour_cor<-cor.test(as.vector(env_theta),as.vector(armour_theta),method = "spearman")$estimate
gill_cor<-cor.test(as.vector(env_theta),as.vector(gill_theta),method = "spearman")$estimate

# delta corrs
shape_cor2<-cor.test(as.vector(env_deltaL),as.vector(shape_deltaL),method = "spearman")$estimate
armour_cor2<-cor.test(as.vector(env_deltaL),as.vector(armour_deltaL),method = "spearman")$estimate
gill_cor2<-cor.test(as.vector(env_deltaL),as.vector(gill_deltaL),method = "spearman")$estimate

# Make plot dataframe
plot_dd<-data.frame(rad_labs=rep(c("Alaska & BC","Alaska & Iceland","Alaska & Scotland","BC & Iceland","BC & Scotland","Iceland & Scotland"),6),
                    env=c(rep(as.vector(env_theta),3),rep(as.vector(env_deltaL),3)),
                    var=c(as.vector(shape_theta),
                          as.vector(armour_theta),
                          as.vector(gill_theta),
                          as.vector(shape_deltaL),
                          as.vector(armour_deltaL),
                          as.vector(gill_deltaL)),
                    measure=rep(c("theta","Delta~L"),each=18),
                    pheno=rep(rep(c("Shape","Armour","Trophic"),each=6),2),
                    cor_val=rep(c(shape_cor,armour_cor,gill_cor,
                          shape_cor2,armour_cor2,gill_cor2),each=6))

plot_dd$pheno<-factor(plot_dd$pheno,levels=c("Shape","Armour","Trophic"))
plot_dd$measure<-factor(plot_dd$measure,levels=c("theta","Delta~L"))

labels<-paste0("rho==",round(plot_dd$cor_val[seq(1,length(plot_dd$cor_val),6)],2)[c(1,4,2,5,3,6)])

spearman_plots<-ggplot(plot_dd,aes(x=env,y=var))+
  geom_point(size=5)+
  geom_label_repel(aes(label=rad_labs),
                   point.padding = 0.5,
                   size =5)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title = element_text(size=22),
        axis.text = element_text(size=20),
        strip.text=element_text(size=18))+
  xlab(expression("Environment"))+
  ylab(expression("Phenotype"))+
  facet_wrap(pheno~measure,scales="free",labeller = label_parsed,ncol=2)+
  annotate(geom="text",x=Inf,y=Inf,label=labels,parse=T,size=10, vjust=4, hjust=1)
                    
pdf("figs/Spearman_cors_envpheno2.pdf",width=12,height=16)                    
print(spearman_plots)
dev.off()

##### Correlation Dendrograms of env v pheno #####
# Individual Env Vars
Alaska_mat<-dd[dd$RADIATION == "ALASKA",
               c("Ca","Gyro","Na","pH","Schisto","Zn",
                 "SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
                 "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
                 "Gill.Raker.N","Resid.Raker.L")]

BC_mat<-dd[dd$RADIATION == "BC",
               c("Ca","Gyro","Na","pH","Schisto","Zn",
                 "SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
                 "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
                 "Gill.Raker.N","Resid.Raker.L")]

Ice_mat<-dd[dd$RADIATION == "ICELAND",
               c("Ca","Gyro","Na","pH","Schisto","Zn",
                 "SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
                 "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
                 "Gill.Raker.N","Resid.Raker.L")]

Scot_mat<-dd[dd$RADIATION == "SCOTLAND",
               c("Ca","Gyro","Na","pH","Schisto","Zn",
                 "SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
                 "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
                 "Gill.Raker.N","Resid.Raker.L")]

mat_list<-list(Alaska_mat,BC_mat,Ice_mat,Scot_mat)
title_vec<-c("Alaska","BC","Iceland","Scotland")

# For each matrix compute the correlation matrix and reduce down so x is env and y is pheno
cor_figs<-lapply(1:length(mat_list),function(x){
  
  # Get Corrs
  tmp_mat<-mat_list[[x]]
  res_tmp <- round(cor(tmp_mat),2)
  
  # Filter
  res_tmp2<-res_tmp[c("Ca","Gyro","Na","pH","Schisto","Zn"),
          c("SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
            "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
            "Gill.Raker.N","Resid.Raker.L")]
  
  # Rename
  colnames(res_tmp2)<-c("Shape PC1","Shape PC2","Shape PC3",
                        "Armour PC1","Armour PC2","Armour PC3","Armour PC4",
                        "Gill Raker N","Gill Raker L")

  # Plot
  melted_cormat <- melt(res_tmp2)
  g_cor <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +  
    theme_bw()+
    geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
    theme(axis.text = element_text(size=18),
          axis.title = element_blank(),
          title = element_text(size=20),
          legend.position = "none")+
    ggtitle(title_vec[x])
  
  return(list(g_cor,res_tmp2))
})

# Plot figure
pdf("figs/Cor_matrices_env_pheno_allRads.pdf",width=12,height=8)
print(ggarrange(cor_figs[[1]][[1]],
          cor_figs[[2]][[1]],
          cor_figs[[3]][[1]],
          cor_figs[[4]][[1]],ncol=2,nrow=2))
dev.off()

# -------------------------
# Env PCs
# -------------------------

Alaska_mat<-dd[dd$RADIATION == "ALASKA",
               c("Env.PC1","Env.PC2","Env.PC3","Env.PC4",
                 "SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
                 "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
                 "Gill.Raker.N","Resid.Raker.L")]

BC_mat<-dd[dd$RADIATION == "BC",
           c("Env.PC1","Env.PC2","Env.PC3","Env.PC4",
             "SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
             "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
             "Gill.Raker.N","Resid.Raker.L")]

Ice_mat<-dd[dd$RADIATION == "ICELAND",
            c("Env.PC1","Env.PC2","Env.PC3","Env.PC4",
              "SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
              "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
              "Gill.Raker.N","Resid.Raker.L")]

Scot_mat<-dd[dd$RADIATION == "SCOTLAND",
             c("Env.PC1","Env.PC2","Env.PC3","Env.PC4",
               "SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
               "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
               "Gill.Raker.N","Resid.Raker.L")]

mat_list<-list(Alaska_mat,BC_mat,Ice_mat,Scot_mat)
title_vec<-c("Alaska","BC","Iceland","Scotland")

# For each matrix compute the correlation matrix and reduce down so x is env and y is pheno
cor_figs<-lapply(1:length(mat_list),function(x){
  
  # Get Corrs
  tmp_mat<-mat_list[[x]]
  res_tmp <- round(cor(tmp_mat),2)
  
  # Filter
  res_tmp2<-res_tmp[c("Env.PC1","Env.PC2","Env.PC3","Env.PC4"),
                    c("SHAPE.IND.PC1","SHAPE.IND.PC2","SHAPE.IND.PC3",
                      "ARMOUR.IND.WITHPLATE.PC1","ARMOUR.IND.WITHPLATE.PC2","ARMOUR.IND.WITHPLATE.PC3","ARMOUR.IND.WITHPLATE.PC4",
                      "Gill.Raker.N","Resid.Raker.L")]
  
  # Rename
  colnames(res_tmp2)<-c("Shape PC1","Shape PC2","Shape PC3",
                        "Armour PC1","Armour PC2","Armour PC3","Armour PC4",
                        "Gill Raker N","Gill Raker L")
  
  rownames(res_tmp2)<-c("Env PC1","Env PC2","Env PC3","Env PC4")
  
  # Plot
  melted_cormat <- melt(res_tmp2)
  g_cor <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +  
    theme_bw()+
    geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
    theme(axis.text = element_text(size=18),
          axis.title = element_blank(),
          title = element_text(size=20),
          legend.position = "none")+
    ggtitle(title_vec[x])
  
  return(list(g_cor,res_tmp2))
})

# Plot figure
pdf("figs/Cor_matrices_env_pheno_allRads_EnvPCs.pdf",width=12,height=8)
print(ggarrange(cor_figs[[1]][[1]],
                cor_figs[[2]][[1]],
                cor_figs[[3]][[1]],
                cor_figs[[4]][[1]],ncol=2,nrow=2))
dev.off()


