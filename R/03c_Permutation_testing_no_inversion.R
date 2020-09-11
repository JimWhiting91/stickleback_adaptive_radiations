#################
# This script permutes expected distributions for convergence across radiations according to randomised
# attainment of "randomly" associated SNPs
#################

lib<-as.vector(c("data.table","VennDiagram","ggplot2","parallel","dplyr"))
lapply(lib,library,character.only=TRUE)

# ------
# This script will be repeated over windows and genes so we put in function
# ------

# Define the grouping
grouping_vector<-as.vector(c("50k","75k","100k","200k","cM"))

# Define the variables
variables_vector<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn",
                              "Shape_PC1","Shape_PC2","Shape_PC3",
                              "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                              "Gill_Raker_L","Gill_Raker_N"))
env_pheno_vec<-c(rep("Env",6),
                 rep("Shape",3),
                 rep("Armour",7),
                 rep("Gill",2))
# Define radiations
radiations_vector<-as.vector(c("Alaska","BC","Iceland","Scotland"))

# Define iterations
iterations<-10000

# Run over all window/gene sizes
lapply(1:length(grouping_vector),function(x){
  
  # Read in data for SNP counts for each radiation and save
  for (i in 1:length(radiations_vector)){
    dd_SNP_count<-read.table(paste0("data/outlier_beds/",radiations_vector[i],"_",grouping_vector[x],"_windows_SNPcount.bed"),header=F,sep='\t',fill = TRUE)
    dd_SNP_count$window_id<-paste0(dd_SNP_count$V1,":",dd_SNP_count$V2,"-",dd_SNP_count$V3)
    colnames(dd_SNP_count)<-c("chr","BP1","BP2","SNP_N","window_id")
    assign(paste0(radiations_vector[i],"_dd_SNP_count"),dd_SNP_count)
  }
  
  # Save each SNP count file to a list to reference elements
  SNP_count_list<-list(Alaska_dd_SNP_count,
                       BC_dd_SNP_count,
                       Iceland_dd_SNP_count,
                       Scotland_dd_SNP_count)
  
  # Make output for saving
  SNP_based_output<-list()
  region_based_output<-list()

  # Permute across variables
  variable_function<-function(y){
    
    # For the variable we need to know the N of associated SNPs per radiation
    # Read in data for outlier SNPs
    for (i in 1:length(radiations_vector)){
      dd_outliers<-read.table(paste0("data/outlier_beds/",grouping_vector[x],"/",grouping_vector[x],"_windows_",radiations_vector[i],"_",variables_vector[y],"_cleaned_outliers.bed_sorted.bed"),
                              header=F,sep='\t',fill = T)
      dd_outliers$window_id<-paste0(dd_outliers$V1,":",dd_outliers$V2,"-",dd_outliers$V3)
      colnames(dd_outliers)<-c("chr","BP1","BP2","SNP_N","window_id")
      assign(paste0(radiations_vector[i],"_dd_outliers"),dd_outliers)
    }
    
    # Save each outliers file to a list to be referenced later
    outlier_list<-list(Alaska_dd_outliers,
                         BC_dd_outliers,
                         Iceland_dd_outliers,
                         Scotland_dd_outliers)
    
    # We also want the number of associated regions/genes
    # Read in data for associated region/gene counts
    dd_regions<-read.table(paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_",variables_vector[y],"_top_outliers_by_rad.txt"),header=T)
    
    #############################################################
    # UNIQUE TO THIS SCRIPT, REMOVE THE CHR1 INVERSION
    if(grouping_vector[x] != "cM"){
      to_remove<-as.integer(rownames(dd_regions[dd_regions$chr=="groupI" &
                              dd_regions$BP1 < 21800000 &
                              dd_regions$BP2 > 21550000,]))
      to_keep<-1:nrow(dd_regions)
      to_keep<-to_keep[!(to_keep %in% to_remove)]
      dd_regions<-dd_regions[to_keep,]
    } else {
      to_remove<-as.integer(rownames(dd_regions[dd_regions$chr=="groupI" &
                                                  dd_regions$BP1 < 59.72684 &
                                                  dd_regions$BP2 > 55.43045,]))
      to_keep<-1:nrow(dd_regions)
      to_keep<-to_keep[!(to_keep %in% to_remove)]
      dd_regions<-dd_regions[to_keep,]
    }
    
    
    # Split data file up by radiation
    for (i in 1:length(levels(dd_regions$radiation))){
      regions_tmp<-dd_regions[dd_regions$radiation == levels(dd_regions$radiation)[i],]
      assign(paste0(levels(dd_regions$radiation)[i],"_dd_regions"),regions_tmp)
    }
    
    # Save each outliers file to a list to be referenced later
    regions_list<-list(Alaska_dd_regions,
                       BC_dd_regions,
                       Iceland_dd_regions,
                       Scotland_dd_regions)
    
    # Run the permutation with calculate.overlap each time
    # The desired output is a list of lists of 11 elements, each showing overlap in all 4 (NIBA), 3 and 2.
    # Only do this part if there is not already an output
    if(file.exists(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_bySNP_",iterations,"runs_FULLSIM_NO_INVERSION.txt"))=="FALSE"){
      
    perm_output_SNP_based<-mclapply(1:iterations,function(z){
      
      # This gets run for each iteration of the permutation
        # For each radiation, we need to take the total SNP count file, expand windows out by their SNP N
        # Randomise order, take head = to associated SNP N, count occurrence of each window in head, and return value to original SNP count
        
        randomised_regions_list<-list()
        for (i in 1:length(SNP_count_list)){
          
          # Fetch our inputs
          SNP_count_tmp<-SNP_count_list[[i]]
          outlier_tmp<-outlier_list[[i]]
          
          # Find the sum of outlier SNPs
          total_outliers<-sum(outlier_tmp$SNP_N)
          
          # Expand out SNP_Count by SNP N
          expanded_windows<-rep(SNP_count_tmp$window_id,SNP_count_tmp$SNP_N)
          
          # Randomise and take 'head' of length "total_outliers"
          random_associated<-head(sample(expanded_windows),total_outliers)
          
          # Count occurrence of each window in random_associated and return data.frame
          random_associated_df<-as.data.frame(table(random_associated))
          colnames(random_associated_df)<-c("window_id","Associated_N")
          
          # Combine random_associated_df with original SNP_count file
          SNP_count_tmp<-merge(SNP_count_tmp,random_associated_df,by=c("window_id"),all.x=T)
          SNP_count_tmp[is.na(SNP_count_tmp)] <- 0
          
          # Calculate the binomial expectation
          p<-total_outliers/sum(SNP_count_tmp$SNP_N)
          SNP_count_tmp$exp<-qbinom(0.99,SNP_count_tmp$SNP_N,p)
          
          # Extract regions that exceed expectation and assign to list
          randomised_regions_list[[i]]<-SNP_count_tmp[SNP_count_tmp$Associated_N > SNP_count_tmp$exp,"window_id"]
        }
        
        # Now calculate overlap for each radiation
        # Calculate overlap of windows across radiations
        perm_overlap<-calculate.overlap(randomised_regions_list)
        names(perm_overlap)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
        
        # Turn calculate.overlap into human readable output
        perm_overlap_counts<-list()
        for (i in seq(1,15,by=1)){
          perm_overlap_counts[[i]]<-length((as.vector(perm_overlap[[i]])))
        }
        perm_overlap_counts<-as.data.table(perm_overlap_counts)
        colnames(perm_overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
        
        return(perm_overlap_counts)
      },mc.cores = detectCores()-1)
    
    # Combine the outputs of each permutation into a single data.table  
    perm_output_SNP_all<-data.frame(rbindlist(perm_output_SNP_based))
    
    # Save all SIM for future analysis
    dir.create(paste0("outputs/",grouping_vector[x],"/03"),showWarnings = F)
    write.table(perm_output_SNP_all,
                paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_bySNP_",iterations,"runs_FULLSIM_NO_INVERSION.txt"),
                row.names = F,
                quote = F)
    }
    
    # Calculate real overlap
    rads<-unique(dd_regions$radiation)
    rad_outliers<-lapply(rads,function(rad){return(dd_regions[dd_regions$radiation==rad,"window_id"])})
    real_overlap<-calculate.overlap(rad_outliers)
    
    # Turn calculate.overlap into human readable output
    real_overlap_counts<-list()
    for (i in seq(1,15,by=1)){
      real_overlap_counts[[i]]<-length((as.vector(real_overlap[[i]])))
    }
    real_overlap_counts<-as.data.table(real_overlap_counts)
    colnames(real_overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    obs_dd<-data.frame(real_overlap_counts)
    
    # Calculate p-values by moving across columns 
    p_val_list<-list()
    for (i in 1:11){
      column_interest<-perm_output_SNP_all[,i]
      p_val<-length(perm_output_SNP_all[column_interest >= obs_dd[,i],i])
      p_val_list[[i]]<-(p_val/iterations)
    }
    
    # Unlist to get vector
    p_val_dd<-t(data.frame(unlist(p_val_list)))
    colnames(p_val_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN")
    rownames(p_val_dd)<-variables_vector[y]
    
    # Save output
    write.table(p_val_dd,
                paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_bySNP_",iterations,"runs_pvals_NO_INVERSION.txt"),
                quote = F)
    
  # We also want to repeat the analysis using region/gene numbers rather than associated SNPs
  if(file.exists(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_",iterations,"runs_FULLSIM_NO_INVERSION.txt"))=="FALSE"){
      
    perm_output_region_based<-mclapply(1:iterations,function(z){
    
    randomised_regions_list<-list()
    for (i in 1:length(SNP_count_list)){
      
      # Fetch our inputs
      SNP_count_tmp<-SNP_count_list[[i]]
      regions_tmp<-regions_list[[i]]
      
      # Find the sum of outlier regions
      total_outliers<-length(regions_tmp$chr)

      # Randomise and take 'head' of length "total_outliers"
      randomised_regions_list[[i]]<-head(sample(SNP_count_tmp$window_id),total_outliers)
    }
    
    # Now calculate overlap for each radiation
    # Calculate overlap of windows across radiations
    perm_overlap<-calculate.overlap(randomised_regions_list)
    names(perm_overlap)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    # Turn calculate.overlap into human readable output
    perm_overlap_counts<-list()
    for (i in seq(1,15,by=1)){
      perm_overlap_counts[[i]]<-length((as.vector(perm_overlap[[i]])))
    }
    perm_overlap_counts<-as.data.table(perm_overlap_counts)
    colnames(perm_overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    return(perm_overlap_counts)
    
  },mc.cores=detectCores()-1)
  
  # Rbindlist the permuted outputs
  perm_output_region_all<-data.frame(rbindlist(perm_output_region_based))
  
  # Save all SIM for future analysis
  write.table(perm_output_region_all,
              paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_",iterations,"runs_FULLSIM_NO_INVERSION.txt"),
              row.names = F,
              quote = F)
  }
  
    # Calculate real overlap
    rads<-unique(dd_regions$radiation)
    rad_outliers<-lapply(rads,function(rad){return(dd_regions[dd_regions$radiation==rad,"window_id"])})
    real_overlap<-calculate.overlap(rad_outliers)
    
    # Turn calculate.overlap into human readable output
    real_overlap_counts<-list()
    for (i in seq(1,15,by=1)){
      real_overlap_counts[[i]]<-length((as.vector(real_overlap[[i]])))
    }
    real_overlap_counts<-as.data.table(real_overlap_counts)
    colnames(real_overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    obs_dd<-data.frame(real_overlap_counts)
    
  
  # Calculate p-values by moving across columns 
  p_val_list<-list()
  for (i in 1:11){
    column_interest<-perm_output_region_all[,i]
    p_val<-length(perm_output_region_all[column_interest >= obs_dd[,i],i])
    p_val_list[[i]]<-(p_val/iterations)
  }
  
  # Unlist to get vector
  p_val_dd<-t(data.frame(unlist(p_val_list)))
  colnames(p_val_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN")
  rownames(p_val_dd)<-variables_vector[y]
  
  write.table(p_val_dd,
              paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byRegion_",iterations,"runs_pvals_NO_INVERSION.txt"),
              quote = F)
  }

# Run the function
mclapply(1:length(variables_vector),variable_function,mc.cores=detectCores()-1)

# Read in all files for byRegion and bySNP and combine
region_save_list<-list()
SNP_save_list<-list()
for (i in 1:length(variables_vector)){
  region_save_list[[i]]<-data.table(read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[i],"_byRegion_",iterations,"runs_pvals_NO_INVERSION.txt"),header=T))
  SNP_save_list[[i]]<-data.table(read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[i],"_bySNP_",iterations,"runs_pvals_NO_INVERSION.txt"),header=T))
}

# Rbindlist each output and save as a final output
region_save_dd<-data.frame(rbindlist(region_save_list))
rownames(region_save_dd)<-variables_vector

# Calculate FDR values and save alongside
tmp_vector<-as.vector(as.matrix(region_save_dd))
tmp_FDR<-p.adjust(tmp_vector,method = 'fdr')
tmp_mat<-data.frame(matrix(tmp_vector,ncol=length(colnames(region_save_dd))))
region_save_dd<-cbind(region_save_dd,tmp_mat)

colnames(region_save_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN",
                            "ABIN_FDR","ABI_FDR","ABN_FDR","AIN_FDR","BIN_FDR","AB_FDR","AI_FDR","AN_FDR","BI_FDR","BN_FDR","IN_FDR")

write.table(region_save_dd,
            paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_ALL_byRegion_",iterations,"runs_pvals_NO_INVERSION.txt"),
            quote=F)

# Rbindlist each output and save as a final output
SNP_save_dd<-data.frame(rbindlist(SNP_save_list))
rownames(SNP_save_dd)<-variables_vector

# Calculate FDR values and save alongside
tmp_vector<-as.vector(as.matrix(SNP_save_dd))
tmp_FDR<-p.adjust(tmp_vector,method = 'fdr')
tmp_mat<-data.frame(matrix(tmp_vector,ncol=length(colnames(SNP_save_dd))))
SNP_save_dd<-cbind(SNP_save_dd,tmp_mat)

colnames(SNP_save_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN",
                            "ABIN_FDR","ABI_FDR","ABN_FDR","AIN_FDR","BIN_FDR","AB_FDR","AI_FDR","AN_FDR","BI_FDR","BN_FDR","IN_FDR")

write.table(SNP_save_dd,
            paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_ALL_bySNP_",iterations,"runs_pvals_NO_INVERSION.txt"),
            quote = F)
})

##### Read back in data and calculate the Expected for 2 or more (ie. Not grouping specific) #####
windows_out<-data.frame(rbindlist(lapply(1:length(grouping_vector),function(x){
  
  variables_out<-data.frame(rbindlist(lapply(1:length(variables_vector),function(y){
    
    # Read in and calculate means
    sim_dd<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_10000runs_FULLSIM_NO_INVERSION.txt"),header=T)
    sim_dd<-sim_dd[,1:11]
    
    # Read in data for associated region/gene counts
    dd_regions<-read.table(paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_",variables_vector[y],"_top_outliers_by_rad.txt"),header=T)
    
    #############################################################
    # UNIQUE TO THIS SCRIPT, REMOVE THE CHR1 INVERSION
    if(grouping_vector[x] != "cM"){
      to_remove<-as.integer(rownames(dd_regions[dd_regions$chr=="groupI" &
                                                  dd_regions$BP1 < 21800000 &
                                                  dd_regions$BP2 > 21550000,]))
      to_keep<-1:nrow(dd_regions)
      to_keep<-to_keep[!(to_keep %in% to_remove)]
      dd_regions<-dd_regions[to_keep,]
    } else {
      to_remove<-as.integer(rownames(dd_regions[dd_regions$chr=="groupI" &
                                                  dd_regions$BP1 < 59.72684 &
                                                  dd_regions$BP2 > 55.43045,]))
      to_keep<-1:nrow(dd_regions)
      to_keep<-to_keep[!(to_keep %in% to_remove)]
      dd_regions<-dd_regions[to_keep,]
    }
  
    # Calculate real overlap
    rads<-unique(dd_regions$radiation)
    rad_outliers<-lapply(rads,function(rad){return(dd_regions[dd_regions$radiation==rad,"window_id"])})
    real_overlap<-calculate.overlap(rad_outliers)
    
    # Turn calculate.overlap into human readable output
    real_overlap_counts<-list()
    for (i in seq(1,15,by=1)){
      real_overlap_counts[[i]]<-length((as.vector(real_overlap[[i]])))
    }
    real_overlap_counts<-as.data.table(real_overlap_counts)
    colnames(real_overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    obs_dd<-data.frame(real_overlap_counts)[,1:11]
    
    
    # Sum Exp and Obs with UL and LL across all
    out_dd<-data.frame(Exp=mean(rowSums(sim_dd)),
                       UL=quantile(rowSums(sim_dd),probs = 0.95),
                       Obs=rowSums(obs_dd),
                       var=variables_vector[y],
                       window=grouping_vector[x],
                       group=env_pheno_vec[y],
                       p.value=length(rowSums(sim_dd)[rowSums(sim_dd) >= sum(obs_dd)])/iterations)
    
    # Return
    return(out_dd)
})))
  
})))

# Save results
write.table(windows_out,
            "outputs/Pairwise_permutation_results_no_inversion.txt",
            row.names=F,quote = F,sep="\t")
