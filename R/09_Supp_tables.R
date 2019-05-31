################################################################
# Making Supplementary Tables
#################################################################

# Load packages
lib<-as.vector(c("dplyr","lsmeans","ggsignif","MASS","data.table","VennDiagram","ggplot2","parallel","effects","lme4","lmerTest"))
lapply(lib,library,character.only=TRUE)

# Source functions
source("R/Genomic_Parallelism_Rcode_functions.R")

# Run over all the window sizes
winds<-c("50k","75k","100k","200k","cM")
wind_size<-c(50000,75000,100000,200000,0.1)

rads<-c("Alaska","BC","Iceland","Scotland")

# Run over variables
variables_vector<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn",
                              "Shape_PC1","Shape_PC2","Shape_PC3",
                              "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                              "Gill_Raker_L","Gill_Raker_N"))

##### TABLE S4 ######
# This table shows all associated windows and the variable and region associated

# We make a table per window
S4_list<-lapply(1:length(winds),function(x){
  
  # Read in all the data for the window
  outliers<-data.frame(rbindlist(lapply(1:length(variables_vector),function(y){
    
    # Read in variable
    tmp<-read.table(paste0("outputs/",winds[x],"/",winds[x],"_",variables_vector[y],"_top_outliers_by_rad.txt"),header=T)
    
    # Replace Nuist with Scotland
    tmp$radiation<-as.character(tmp$radiation)
    tmp[tmp$radiation == "Nuist","radiation"]<-"Scotland"
    
    # Format
    tmp2<-data.frame(Chr=tmp$chr,
                     Start=tmp$BP1,
                     End=tmp$BP2,
                     ID=tmp$window_id)
    
    # Concat within variable
    tmp3<-data.frame(tmp %>%
      group_by(window_id) %>%
      summarise(radiation=paste(unique(radiation),collapse = '/')))
    
    # Retain only unique tmp2 rows
    tmp2<-tmp2[!duplicated(tmp2$ID),]
    tmp2$ID<-as.character(tmp2$ID)
    
    # Merge
    tmp2$Associated<-NA
    for(i in 1:nrow(tmp3)){
      tmp2[tmp2$ID == as.character(tmp3$window_id[i]),"Associated"]<-paste0(variables_vector[y]," (",tmp3$radiation[i],")")
    }
    
    # Return
    return(tmp2)
    # End read
  })))
  
  # Sort and concatenate the associated column
  outliers_sort<-outliers[order(outliers$Chr,as.integer(outliers$Start)),]
  outliers_cat<-data.frame(outliers_sort %>%
            group_by(ID) %>%
            summarise(Associated=paste(unique(Associated),collapse = ' / ')))
  
  # Return to outliers_sort
  outliers_sort<-outliers_sort[!duplicated(outliers_sort$ID),]
  outliers_sort$Associated<-NA
  for(i in 1:nrow(outliers_cat)){
    outliers_sort[outliers_sort$ID == as.character(outliers_cat$ID[i]),"Associated"]<-outliers_cat$Associated[i]
  }
  
  # Add on window label
  outliers_sort$Window<-winds[x]
  
  # Compare to bed file of previous Marine x Freshwater outliers
  MxF<-read.csv("data/marine_fresh_windows.csv")
  
  # Do they overlap?
  outliers_sort$MxF<-"No"
  for(i in 1:nrow(outliers_sort)){
    MxF_tmp<-MxF[as.character(MxF$group)==as.character(outliers_sort$Chr[i]) &
                   MxF$Begin <= as.integer(as.character(outliers_sort$Start[i])) &
                   MxF$End >= as.integer(as.character(outliers_sort$End[i])),"Source"]
    MxF_tmp<-unique(as.character(MxF_tmp))
    if(length(MxF_tmp) > 0){
      outliers_sort$MxF[i]<-paste(MxF_tmp,collapse="/")
    }
  }
  
  # End
  return(outliers_sort)
})

# Cat all outputs together
S4_final<-data.frame(rbindlist(S4_list))

# Save
write.table(S4_final,
            "tables/TableS4_All_associated_windows_al_window_sizes.txt",
            quote = F,row.names = F,sep = "\t")
write.csv(S4_final,
            "tables/TableS4_All_associated_windows_al_window_sizes.csv")

##### Table S5 #####
# Table S5 shows how many associated windows are associated across different variables

# We make a matrix per window
S5_list<-lapply(1:length(winds),function(x){
  
  # Read in all the data for the window
  outliers<-data.frame(rbindlist(lapply(1:length(variables_vector),function(y){
    
    # Read in variable
    tmp<-read.table(paste0("outputs/",winds[x],"/",winds[x],"_",variables_vector[y],"_top_outliers_by_rad.txt"),header=T)
    
    # Replace Nuist with Scotland
    tmp$radiation<-as.character(tmp$radiation)
    tmp[tmp$radiation == "Nuist","radiation"]<-"Scotland"
    
   
    return(tmp)
    # End read
  })))
  
  # Also fetch MxF results
  MxF<-read.csv("tables/TableS6_MxF_outliers_50kb.csv")
  MxF<-MxF[MxF$N > 1,]
  outliers<-outliers[,c("variable","window_id")]
  MxF$variable<-"MxF"
  MxF<-MxF[,c("variable","window_id")]
  outliers<-rbind(outliers,MxF)
  variables_vector2<-c("MxF",variables_vector)
  
  # For each variable we need to perform all pairwise comps
  comps<-combn(1:19,2)
  out_mat<-matrix(nrow=19,ncol=19,0)
  
  
  # Fill the matrix
  for(i in 1:length(comps[1,])){
    vec1<-as.character(outliers[outliers$variable == variables_vector2[comps[1,i]],"window_id"])
    vec2<-as.character(outliers[outliers$variable == variables_vector2[comps[2,i]],"window_id"])
  
    # Only keep duplicates in each to get parallel regions
    if(comps[1,i] != 1){
    vec1<-vec1[duplicated(vec1)]
    }
    vec2<-vec2[duplicated(vec2)]
    
    # Overlap
    overlap<-length(Reduce(intersect,list(vec1,vec2)))
    
    # Return to matrix
    out_mat[comps[1,i],comps[2,i]]<-overlap
  }
  
  # Rename cols etc
  colnames(out_mat)<-variables_vector2
  rownames(out_mat)<-variables_vector2
  
  out_mat_tri<-out_mat[upper.tri(out_mat)]
  
  # End
  return(out_mat)
})

S5_list[[1]]

##### Table S6 #######
# Table with expected results fsrom 10k sims

# Just 50kb

  # Read in the sims
  sims_out<-lapply(1:18,function(y){
    
    # Read in and calculate means
    sim_dd<-read.table(paste0("outputs/",winds[x],"/03/",winds[x],"_",variables_vector[y],"_byREGION_10000runs_FULLSIM.txt"),header=T)

    return(sim_dd[,1:11])
  })

# For each element of list, return the colmean
expecteds<-lapply(1:length(sims_out),function(x){
  tmp<-sims_out[[x]]
  exp_out<-colMeans(tmp)
  return(exp_out)
})

# Make table
out_mat<-matrix(ncol=11,nrow = 18)
for(i in 1:18){
  out_mat[i,]<-expecteds[[i]]
}

# Tidy
rownames(out_mat)<-variables_vector

# Save
write.table(out_mat,
            "tables/S6_expected_vals.txt",quote = F,sep = "\t")

##### Table S12 #####
# Table shows total N of SNPs, windows, and overlap of windows

# Perform over windows
S12_list<-data.frame(rbindlist(lapply(1:length(winds),function(x){
  
  # Read in the data
  dd_tmp<-data.frame(rbindlist(lapply(1:length(rads),function(y){
    
    # Data for SNP counts
    tmp<-read.table(paste0("data/outlier_beds/",rads[y],"_",winds[x],"_windows_SNPcount.bed"))
    
    # Tidy
    colnames(tmp)<-c("chr","BP1","BP2","SNP_N")
    tmp$window_id<-paste0(tmp$chr,":",tmp$BP1,"-",tmp$BP2)
    tmp$rad<-rads[y]
    return(tmp)
  })))
  
  # Get Window Count for each Rad
  window_N<-data.frame(dd_tmp %>%
                         group_by(rad) %>%
                         summarise(window_N=length(unique(window_id))))
  
  # We also want to know the overlap of windows
  window_venn<-calculate.overlap(list(dd_tmp[dd_tmp$rad==rads[1],"window_id"],
                                      dd_tmp[dd_tmp$rad==rads[2],"window_id"],
                                      dd_tmp[dd_tmp$rad==rads[3],"window_id"],
                                      dd_tmp[dd_tmp$rad==rads[4],"window_id"]))
  
  # Turn calculate.overlap into human readable output
  overlap_counts<-list()
  for (i in seq(1,15,by=1)){
    overlap_counts[[i]]<-length((as.vector(window_venn[[i]])))
  }
  overlap_counts<-as.data.frame(overlap_counts)
  colnames(overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
  
  overlap_counts$window_size<-winds[x]
  
  return(overlap_counts)
  
})))

# Write to a table
write.table(S12_list,
            "tables/TableS12_window_counts.txt",
            row.names = F,quote = F,sep="\t")

##### Pairwise FST Heatmap #####
library(grid)
library(ggplot2)
library(gridExtra)
library(gtable)
library(reshape)

#Set working directory 
dd<-read.csv("outputs/All_Radiations_All_Pops_Fst.csv")

# First calculate Pairwise FST of each matrix
dd2<-dd[,-1]
Scot_Ice<-mean(as.matrix(dd2[1:18,19:36]))
Scot_BC<-mean(as.matrix(dd2[1:18,37:54]))
Scot_Alas<-mean(as.matrix(dd2[1:18,55:73]))
Ice_BC<-mean(as.matrix(dd2[19:36,37:54]))
Ice_Alas<-mean(as.matrix(dd2[19:36,55:73]))
BC_Alas<-mean(as.matrix(dd2[37:54,55:73]))

Fst_out<-data.frame(Comps=c("Alaska-BC",
                            "Alaska-Iceland",
                            "Alaska-Scotland",
                            "BC-Iceland",
                            "BC-Scotland",
                            "Iceland-Scotland"),
                    FST=c(BC_Alas,
                          Ice_Alas,
                          Scot_Alas,
                          Ice_BC,
                          Scot_BC,
                          Scot_Ice))

# Write
write.table(Fst_out,
            "tables/TableSX_Avg_Radiation_Fst_comparions.txt",
            row.names = F,quote = F,sep="\t")

#Convert to Matrix
dd.m<-melt(dd)

#Turn your 'treatment' column into a character vector
dd.m$X <- as.character(dd.m$X)
#Then turn it back into an ordered factor
dd.m$X <- factor(dd.m$X, levels=unique(dd.m$X))


#Base Plot
p <- ggplot(dd.m, aes(variable, X)) + geom_tile(aes(fill = value),
                                                colour = "white") + scale_fill_gradient(low = "white",
                                                                                        high = "red")+
  labs(fill="Fst")+
  xlab("")+
  ylab("")+
  coord_cartesian(xlim=c(0,90,ylim=c(0,100)))+
  theme_bw() +
  theme(axis.text.x=element_text(size=6,angle = 90),
        axis.text.y=element_text(size=6),
        panel.border = element_rect(colour = "white", fill=NA, size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none")

# Create the text Grobs
Text1 = textGrob("Scotland")
Text2 = textGrob("Iceland")
Text4 = textGrob("Alaska")
Text3 = textGrob("BC")

# Text 1
p1 = p + annotation_custom(Text1,  xmin = 81, xmax = 81, ymin = 1, ymax = 18) +
  annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 18, ymax = 18) + #side line
  annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 1, ymax = 1) + #side line
  annotation_custom(linesGrob(), xmin = 75, xmax = 75, ymin = 1, ymax = 18) #vert line

# # Text 1
# p1 = p1 + annotation_custom(Text1,  ymin = 76, ymax = 76, xmin = 1, xmax = 18) +
#   annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 18, xmax = 18) + #side line
#   annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 1, xmax = 1) + #side line
#   annotation_custom(linesGrob(), xmin = 75, xmax = 75, ymin = 1, xmax = 18) #vert line


# Text 2
p1 = p1 + annotation_custom(Text2,  xmin = 80, xmax = 80, ymin = 19, ymax = 36) +
  annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 36, ymax = 36) +
  annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 19, ymax = 19) +
  annotation_custom(linesGrob(), xmin = 75, xmax = 75, ymin = 19, ymax = 36)

# Text 3
p1 = p1 + annotation_custom(Text3,  xmin = 79, xmax = 79, ymin = 37, ymax = 54) +
  annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 54, ymax = 54) +
  annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 37, ymax = 37) +
  annotation_custom(linesGrob(), xmin = 75, xmax = 75, ymin = 37, ymax = 54)


# Text 4
p1 = p1 + annotation_custom(Text4,  xmin = 80, xmax = 80, ymin = 59, ymax = 73) +
  annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 73, ymax = 73) +
  annotation_custom(linesGrob(), xmin = 74, xmax = 75, ymin = 55, ymax = 55) +
  annotation_custom(linesGrob(), xmin = 75, xmax = 75, ymin = 55, ymax = 73)

# Write to output
pdf("figs/All_radiation_Fst_heatmap.pdf",width=8,height=6)
p1
dev.off()

