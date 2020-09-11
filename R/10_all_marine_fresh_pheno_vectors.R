######################################################################################
# Calculating Vectors from Marine phenotypes to each freshwater population

# We need these...
lib<-c("ggplot2","data.table","Morpho","cowplot","dplyr","parallel","wesanderson","viridis")
lapply(lib,library,character.only=T)
source("R/Genomic_Parallelism_Rcode_functions.R")

# Read in phenotype data
phenos<-read.csv("data/files isabel/All_Shape_and_Pheno_data.csv")
pop_phenos<-read.csv("data/All_Lakes_Env_and_Pheno_v2.csv")

marines<-c("MUDL","LICA","NYPS","OBSM")

####################################################################################
# Also use this information to calculate vectors between marine and all fresh pops

# Assmeble raw armour data
# Read in data
armour_dd<-read.table("outputs/PCA_stats/Armour_fresh_marine.txt",header=T)
armour_dd[armour_dd$Country == "Canada_BC","Country"]<-"BC"
armour_vec<-c("DS1","DS2","PS","LP","HP","BAP","n.plate")

# Summarise by pop
armour_dd<-armour_dd[,c("Lake_code","Country","system",
                        armour_vec)]
colnames(armour_dd)[1:3]<-c("pop","rad","type")

# Scale
for(i in 4:ncol(armour_dd)){
  armour_dd[,i]<-scale(armour_dd[,i])
}

# Avg per pop
armour_dd_avg <- data.frame(armour_dd %>% 
  group_by(pop) %>%
  summarise_all(mean))
for(i in 1:nrow(armour_dd_avg)){
  armour_dd_avg$rad[i] <- as.character(armour_dd[armour_dd$pop == armour_dd_avg$pop[i],"rad"])
  armour_dd_avg$type[i] <- as.character(armour_dd[armour_dd$pop == armour_dd_avg$pop[i],"type"])
}

####### Shape
shape_dd<-read.table("data/phenotype_data/body_shape_residuals_separate_radiations_and_marine.txt",header=T)
shape_dd$country <- as.character(shape_dd$country)
shape_dd[shape_dd$country == "Canada","country" ]<- "BC"
shape_vec<-paste0("RegResid",1:26)

# Summarise by pop
shape_dd<-shape_dd[,c("LAKE","country","system",
                        shape_vec)]
colnames(shape_dd)[1:3]<-c("pop","rad","type")

# # Don't scale here, because all are on same unit
# for(i in 4:ncol(shape_dd)){
#   shape_dd[,i]<-scale(shape_dd[,i])
# }

# Avg per pop
shape_dd_avg <- data.frame(shape_dd %>% 
                              group_by(pop) %>%
                              summarise_all(mean))
for(i in 1:nrow(shape_dd_avg)){
  shape_dd_avg$rad[i] <- as.character(shape_dd[shape_dd$pop == shape_dd_avg$pop[i],"rad"])
  shape_dd_avg$type[i] <- as.character(shape_dd[shape_dd$pop == shape_dd_avg$pop[i],"type"])
}

####### Gill Data
gill_dd<-read.csv("data/phenotype_data/gill_rakers_allpops.csv")
gill_dd$Location<-as.character(gill_dd$Location)
gill_dd$Habitat<-"freshwater"
gill_dd[gill_dd$Pop %in% c("NYPS","MUD","LICA"),"Habitat"]<-"marine"  
#gill_dd[gill_dd$Habitat == "marine","Location"]<-"marine"

# Summarise by pop
gill_dd<-gill_dd[,c("Pop","Location","Habitat",
                      "GillRakerL","GillRakerN")]
colnames(gill_dd)[1:3]<-c("pop","rad","type")

## Don't scale here, because all are on same unit
 for(i in 4:ncol(gill_dd)){
   gill_dd[,i]<-scale(gill_dd[,i])
}

gill_avg<-data.frame(gill_dd %>% 
                               group_by(pop) %>% 
                               summarise_all(mean))
gill_avg<-gill_avg[2:nrow(gill_avg),]

for(i in 1:nrow(gill_avg)){
  gill_avg$rad[i] <- as.character(gill_dd[gill_dd$pop == gill_avg$pop[i],"rad"])
  gill_avg$type[i] <- as.character(gill_dd[gill_dd$pop == gill_avg$pop[i],"type"])
}

gill_avg$pop<-as.character(gill_avg$pop)
gill_avg[gill_avg$pop == "MUD","pop"]<-"MUDL"

####### GET VECTORS ##########
phenos_list<-list(armour_dd_avg,shape_dd_avg,gill_avg)
var_vec<-c("Armour","Shape","Gill")
# vectors_out<-lapply(1:length(phenos_list),function(x){
# 
#   # Get vals
#   tmp<-phenos_list[[x]]
# 
#   # What are we using
#   vars<-colnames(tmp)[!(colnames(tmp) %in% c("pop","rad","type"))]
# 
#   # Which marines?
#   marines<-tmp[tmp$type=="marine","pop"]
# 
#   # Calculate vectors for each marine pop to every other
#   marine_out<-data.frame(rbindlist(lapply(marines,function(y){
# 
#     freshs<-as.character(tmp[tmp$type=="freshwater","pop"])
# 
#     fresh_vecs<-data.frame(rbindlist(lapply(freshs,function(z){
# 
#       # Set up exit df
#       vec_out<-data.frame(marine=y,
#                           fresh=z,
#                           marine_rad=tmp[tmp$pop==y,"rad"],
#                           fresh_rad=tmp[tmp$pop==z,"rad"],
#                           tmp[tmp$pop==z,vars]-tmp[tmp$pop==y,vars])
#     })))
# 
#     return(fresh_vecs)
#   })))
# 
#   # Now we need to compare ALL vectors and get angles and length
#   comps<-combn(1:nrow(marine_out),2)
# 
#   vector_comps<-data.frame(rbindlist(mclapply(1:ncol(comps),function(comp){
# 
#     # Get vecs
#     avg_vec_diffs<-data.frame(vec1=t(marine_out[comps[1,comp],vars]),
#                               vec2=t(marine_out[comps[2,comp],vars]))
# 
#     # Calc angle
#     angle.deg<-angle.calc(avg_vec_diffs[,1],avg_vec_diffs[,2])*(180/pi)
# 
#     # Calc length diff
#     get.vectlength <- function(vect){
#       return(sqrt(sum(vect^2, na.rm = TRUE)))
#     }
#     deltaL<-abs(round(get.vectlength(avg_vec_diffs[,1]) - get.vectlength(avg_vec_diffs[,2]), 3))
# 
#     # Out df
#     out<-data.frame(marine1=marine_out[comps[1,comp],"marine"],
#                     marine2=marine_out[comps[2,comp],"marine"],
#                     fresh1=marine_out[comps[1,comp],"fresh"],
#                     fresh2=marine_out[comps[2,comp],"fresh"],
#                     theta_degrees=angle.deg,
#                     deltaL=deltaL)
# 
#     out$comp1<-paste0(out$marine1,"_",out$fresh1)
#     out$comp2<-paste0(out$marine2,"_",out$fresh2)
# 
#     # Return
#     return(out)
#   },mc.cores=detectCores()-1)))
# 
#   # Add label
#   vector_comps$type<-var_vec[x]
#   return(vector_comps)
# })
# 
# # Save R object
# saveRDS(vectors_out,"outputs/raw_variables_marine_fresh_vectors.rds")

# Read back in 
vectors_out<-readRDS("outputs/raw_variables_marine_fresh_vectors.rds")

######### Analyse angles etc and plot ########
plot_dd<-data.frame(rbindlist(vectors_out))

# Fill a matrix, each row represents marine source, each column fresh
alaska_fresh<-as.character(unique(pop_phenos[pop_phenos$RADIATION == "ALASKA","LAKE"]))
bc_fresh<-as.character(unique(pop_phenos[pop_phenos$RADIATION == "BC","LAKE"]))
iceland_fresh<-as.character(unique(pop_phenos[pop_phenos$RADIATION == "ICELAND","LAKE"]))
scotland_fresh<-as.character(unique(pop_phenos[pop_phenos$RADIATION == "SCOTLAND","LAKE"]))
fresh_vec<-list(alaska_fresh,bc_fresh,iceland_fresh,scotland_fresh)

# Only keep vectors where marine > fresh are in the same radiation
to_keep<-c(paste0("MUDL_",alaska_fresh),
           paste0("LICA_",bc_fresh),
           paste0("NYPS_",iceland_fresh),
           paste0("OBSM_",scotland_fresh))

plot_dd<-plot_dd[plot_dd$comp1 %in% to_keep &
                   plot_dd$comp2 %in% to_keep,]

# Plot out histograms
vector_hists<-ggplot(plot_dd,aes(x=theta_degrees))+
  geom_histogram(bins=100)+
  facet_wrap(~type,ncol=1,scales="free_y")

# And separate within/between radiation comparisons
plot_dd2<-plot_dd
plot_dd2$comp_grouping<-"Different"
plot_dd2[plot_dd2$marine1 == plot_dd2$marine2,"comp_grouping"]<-"Same"
plot_dd2$type_F<-factor(plot_dd2$type,levels=c("Armour","Shape","Gill"))

vector_hists2<-ggplot(plot_dd2,aes(x=theta_degrees,fill=comp_grouping))+
  geom_density(alpha=0.5)+
  facet_wrap(~type_F,ncol=1,scales="free_y",strip.position = "right")+
  labs(y="Density",x=expression(Vector~Angle~(Theta)),fill="Same Radiation?")+
  scale_fill_manual(values=wes_palette("Rushmore1")[4:5])+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        legend.position = "top",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        strip.text = element_text(size=20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=20))

vector_hists3<-ggplot(plot_dd2,aes(x=deltaL,fill=comp_grouping))+
  geom_density(alpha=0.5)+
  facet_wrap(~type_F,ncol=1,scales="free",strip.position = "right")+
  labs(y="Density",x=expression(Vector~Length~Difference~(Delta~italic(L))),fill="Same Radiation?")+
  scale_fill_manual(values=wes_palette("Rushmore1")[4:5])+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        legend.position = "top",
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        strip.text = element_text(size=15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=14))

# Do for each variable - these are ordered alaska - scotland
marines<-c("MUDL","LICA","NYPS","OBSM")
angle_mats<-lapply(var_vec,function(x){
  print(x)
  # Get data
  tmp<-plot_dd[plot_dd$type == x,]
  
  # Fill
  to_fill1<-matrix(ncol=4,nrow=4)
  to_fill2<-matrix(ncol=4,nrow=4)
  
  for(i in 1:4){

    # Compare to each
    for (j in 1:4){
    tmp2<-tmp[tmp$marine1 == marines[i] &
                tmp$fresh1 %in% fresh_vec[[i]],]
      
    vec1<-tmp2[tmp2$marine2 == marines[i] &
                 tmp2$fresh2 %in% fresh_vec[[i]],"theta_degrees"]
    vec2<-tmp2[tmp2$marine2 == marines[j] &
                 tmp2$fresh2 %in% fresh_vec[[j]],"theta_degrees"]
    
    # Catch occassions where it is the other way around
    if(sum(vec2)==0 | sum(vec1)==0){
       tmp2<-tmp[tmp$marine2 == marines[i] &
                   tmp$fresh2 %in% fresh_vec[[i]],]
      
      vec1<-tmp2[tmp2$marine1 == marines[i] &
                   tmp2$fresh1 %in% fresh_vec[[i]],"theta_degrees"]
      vec2<-tmp2[tmp2$marine1 == marines[j] &
                   tmp2$fresh1 %in% fresh_vec[[j]],"theta_degrees"]
    }
    
    if(sum(vec1) > 0 & sum(vec2) > 0){
    # Compare
    test <- t.test(vec1,vec2)
    to_fill1[i,j] <- mean(vec2)
    to_fill2[i,j] <- test$p.value
    } else {
      to_fill1[i,j] <- NA
      to_fill2[i,j] <- NA
    }
    }
  }

# to_fill1[lower.tri(to_fill1)]<-NA
# to_fill2[lower.tri(to_fill2)]<-NA
 
rads<-c("Alaska","BC","Iceland","Scotland")
colnames(to_fill1)<-rads
colnames(to_fill2)<-rads
rownames(to_fill1)<-rads
rownames(to_fill2)<-rads

# Return the outputs
return(list(to_fill1,to_fill2))
})

# Show results
angle_mats[[1]]
angle_mats[[2]]
angle_mats[[3]]

# Plot each as a heatmap
angle_plots<-lapply(1:3,function(x){
  
  # Get to plot
  to_plot<-melt(angle_mats[[x]][[1]])
  p_vals<-melt(angle_mats[[x]][[2]])
  p_vals$value<-round(p_vals$value,3)
  
  to_plot$type<-var_vec[x]
  p_vals$type<-var_vec[x]
  
  # retunr
return(list(to_plot,p_vals))
  
})

to_plot<-data.frame(rbindlist(lapply(angle_plots,function(x){return(x[[1]])})))
to_plot$type_F<-factor(to_plot$type,levels =unique(to_plot$type))
p_vals<-data.frame(rbindlist(lapply(angle_plots,function(x){return(x[[2]])})))
p_vals$type_F<-factor(p_vals$type,levels =unique(p_vals$type))

# Also plot "X" where 
p_vals2<-p_vals
p_vals2[is.na(p_vals2$value),"value"]<-0
p_vals2[p_vals2$value > 0.05,"value"]<-"*"
p_vals2[p_vals2$value != "*","value"]<-""
p_vals2[p_vals2$X1 == p_vals2$X2,"value"]<-"-"
#p_vals2<-na.omit(p_vals2)
colnames(p_vals2)[3]<-"asterisk"
p_vals2$value<-p_vals$value

heats<-ggplot(to_plot,aes(to_plot[,1],to_plot[,2],fill=value))+
  geom_tile()+
 scale_fill_gradient(low = "red2",high="white",na.value = "white")+
  #scale_fill_gradientn(na.value="white")+
  facet_wrap(~type_F,strip.position = "right",ncol=1)+
  theme_bw()+
 geom_text(data=p_vals2,aes(label=asterisk),size=8)+
  theme(axis.title=element_blank(),
       axis.text=element_text(size=20),
       panel.grid = element_blank(),
       strip.text = element_text(size=20),
       axis.ticks = element_blank(),
       legend.title=element_text(size=20),
       legend.text=element_text(size=20),
       legend.position="right",
       axis.text.x=element_text(size=20,angle=30,hjust=1))+
  labs(fill=expression(Avg~Theta))

# Plot together with vectors
pdf("figs/vector_results_with_analysis_tests.pdf",width=10,height=7)
plot_grid(vector_hists2,
          heats,ncol=2,
          align="vh",
          labels = "auto",
          axis="tblr")
dev.off()

# And save
vector_figs<-plot_grid(vector_hists2,
          heats,ncol=2,
          align="vh",
          labels = "auto",
          axis="tblr")
saveRDS(list(vector_hists2,heats),
        "figs/All_vector_figures.RDS")




###################################
# Get avgs
plot_dd2 %>% 
  group_by(type_F,comp_grouping) %>%
  dplyr::summarise(avg_theta = mean(theta_degrees))
