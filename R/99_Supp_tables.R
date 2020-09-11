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
variables_vector<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn","Lake_Area",
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
  
  # Get env and pheno outliers, how many in both?
  env_vars <- variables_vector[1:7]
  pheno_vars <- variables_vector[!(variables_vector %in% env_vars)]
  env_outliers <- unlist(lapply(env_vars,function(env){return(outliers[grep(env,outliers$Associated),"ID"])}))
  pheno_outliers <- unlist(lapply(pheno_vars,function(phen){return(outliers[grep(phen,outliers$Associated),"ID"])}))
  env_and_pheno <- unique(env_outliers)[unique(env_outliers) %in% unique(pheno_outliers)]
  
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
  MxF<-read.csv("data/marine_fresh_windows_v2.csv")
  MxF$group<-gsub("chr","group",MxF$group)
  
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
            "tables/TableSX_All_associated_windows_al_window_sizes.txt",
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
  MxF<-read.csv("tables/TableSX_MxF_outliers_50kb.csv")[,1:3]
  MxF<-MxF[MxF$N > 1,]
  outliers<-outliers[,c("variable","window_id")]
  MxF$variable<-"MxF"
  MxF<-MxF[,c("variable","window_id")]
  outliers2<-rbind(outliers,MxF)
  variables_vector2<-c("MxF",variables_vector)
  
  # For each variable we need to perform all pairwise comps
  comps<-combn(1:20,2)
  out_mat<-matrix(nrow=20,ncol=20,0)
  
  
  # Fill the matrix
  for(i in 1:length(comps[1,])){
    vec1<-as.character(outliers2[outliers2$variable == variables_vector2[comps[1,i]],"window_id"])
    vec2<-as.character(outliers2[outliers2$variable == variables_vector2[comps[2,i]],"window_id"])
    
    # Only keep duplicates in each to get parallel regions
    if(comps[1,i] != 1){
      vec1<-vec1[duplicated(vec1)]
    }
    vec2<-vec2[duplicated(vec2)]
    
    # Overlap
    overlap<-length(Reduce(intersect,list(vec1,vec2)))
    
    # Return to matrix
    out_mat[comps[1,i],comps[2,i]]<-overlap
    
    # Also for each add the overall parallel window count along diagognal
    out_mat[comps[1,i],comps[1,i]]<-length(vec1)
    out_mat[comps[2,i],comps[2,i]]<-length(vec2)
  }
  
  # Rename cols etc
  colnames(out_mat)<-variables_vector2
  rownames(out_mat)<-variables_vector2
  
  out_mat_tri<-out_mat[upper.tri(out_mat)]
  
  # End
  return(out_mat)
})

# Actually we only want 50kb
S5_list[[1]]

##### Table S6 #######
# Table with expected results fsrom 10k sims

# Just 50kb
x<-1
# Read in the sims
sims_out<-lapply(1:19,function(y){
  
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
out_mat<-matrix(ncol=11,nrow = 19)
for(i in 1:19){
  out_mat[i,]<-expecteds[[i]]
}

# Tidy
rownames(out_mat)<-variables_vector

# Save
write.table(out_mat,
            "tables/TableSX_expected_vals.txt",quote = F,sep = "\t")

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
