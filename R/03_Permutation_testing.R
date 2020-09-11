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
variables_vector<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn","Lake_Area",
                              "Shape_PC1","Shape_PC2","Shape_PC3",
                              "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                              "Gill_Raker_L","Gill_Raker_N"))
env_pheno_vec<-c(rep("Env",7),
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
      dd_outliers<-read.table(paste0("data/outlier_beds/",grouping_vector[x],"/",grouping_vector[x],"_windows_",radiations_vector[i],"_",variables_vector[y],"_cleaned_outliers_v2.bed"),
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

    # Save each outliers file to a list to be referenced later
    regions_list<-lapply(unique(dd_regions$radiation),function(rad){return(dd_regions[dd_regions$radiation == rad,])})
    
    # Run the permutation with calculate.overlap each time
    # The desired output is a list of lists of 11 elements, each showing overlap in all 4 (NIBA), 3 and 2.
    # Only do this part if there is not already an output
    if(file.exists(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_bySNP_",iterations,"runs_FULLSIM.txt"))=="FALSE"){
      
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
                paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_bySNP_",iterations,"runs_FULLSIM.txt"),
                row.names = F,
                quote = F)
    }
    
    # Read in sim data
    perm_output_SNP_all<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_bySNP_",iterations,"runs_FULLSIM.txt"),header=T)
    
    # Read in observed data and calculate p-values
    obs_dd<-read.table(paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_Venn_overlaps.txt"))
    colnames(obs_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    # Subset for variable of interest
    obs_dd<-obs_dd[y,1:11]
    
    # Calculate p-values by moving across columns 
    p_val_list<-list()
    for (i in 1:11){
      column_interest<-perm_output_SNP_all[,i]
      p_val<-length(perm_output_SNP_all[column_interest >= obs_dd[1,i],i])
      p_val_list[[i]]<-(p_val/iterations)
    }
    
    # Unlist to get vector
    p_val_dd<-t(data.frame(unlist(p_val_list)))
    colnames(p_val_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN")
    rownames(p_val_dd)<-variables_vector[y]
    
    # Save output
    write.table(p_val_dd,
                paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_bySNP_",iterations,"runs_pvals.txt"),
                quote = F)
    
  # We also want to repeat the analysis using region/gene numbers rather than associated SNPs
  if(file.exists(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_",iterations,"runs_FULLSIM.txt"))=="FALSE"){
      
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
              paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_",iterations,"runs_FULLSIM.txt"),
              row.names = F,
              quote = F)
  }
    
  # Read in perm data
  perm_output_region_all<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_",iterations,"runs_FULLSIM.txt"),header=T)
  
  # Read in observed data and calculate p-values
  obs_dd<-read.table(paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_Venn_overlaps.txt"))
  colnames(obs_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
  
  # Subset for variable of interest
  obs_dd<-obs_dd[y,1:11]
  
  # Calculate p-values by moving across columns 
  p_val_list<-list()
  for (i in 1:11){
    column_interest<-perm_output_region_all[,i]
    p_val<-length(perm_output_region_all[column_interest >= obs_dd[1,i],i])
    p_val_list[[i]]<-(p_val/iterations)
  }
  
  # Unlist to get vector
  p_val_dd<-t(data.frame(unlist(p_val_list)))
  colnames(p_val_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN")
  rownames(p_val_dd)<-variables_vector[y]
  
  write.table(p_val_dd,
              paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byRegion_",iterations,"runs_pvals.txt"),
              quote = F)
  }

# Run the function
mclapply(1:length(variables_vector),variable_function,mc.cores=detectCores())

# Read in all files for byRegion and bySNP and combine
region_save_list<-list()
SNP_save_list<-list()
for (i in 1:length(variables_vector)){
  region_save_list[[i]]<-data.table(read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[i],"_byRegion_",iterations,"runs_pvals.txt"),header=T))
  SNP_save_list[[i]]<-data.table(read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[i],"_bySNP_",iterations,"runs_pvals.txt"),header=T))
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
            paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_ALL_byRegion_",iterations,"runs_pvals.txt"),
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
            paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_ALL_bySNP_",iterations,"runs_pvals.txt"),
            quote = F)
})

##### Calculate the expected value for each window/variable/rad #######
# Requires to read back in simulated numbers and calculate

windows_out<-data.frame(rbindlist(lapply(1:length(grouping_vector),function(x){
  
variables_out<-data.frame(rbindlist(lapply(1:length(variables_vector),function(y){
    
    # Read in and calculate means
tmp_dd<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_10000runs_FULLSIM.txt"),header=T)
exps<-data.frame(t(colMeans(tmp_dd)))

# Remove single RADs
exps<-exps[,1:11]

# Add labels
exps$var<-variables_vector[y]
return(exps)

# End of variables_out
  })))

  variables_out$window<-grouping_vector[x]
  return(variables_out)
# End of windows_out lapply
})))

# Write the output
write.table(windows_out,
            "outputs/03_Expected_N_of_shared_windows_ALL.txt",
            quote = F,row.names=F)

##### Combine windows_out with observed values and FDR...
comparisons <- colnames(windows_out)[1:11]
comp_names <- c("All 4",
                "Alaska & BC & Iceland","Alaska & BC & Scotland","Alaska & Iceland & Scotland","BC & Iceland & Scotland",
                "Alaska & BC","Alaska & Iceland","Alaska & Scotland","BC & Iceland","BC & Scotland","Iceland & Scotland")

exp_obs_fdr <- data.frame(rbindlist(lapply(variables_vector,function(x){
  
  tmp <- data.frame(Variable = x,
                    Comparison = comp_names,
                    Expected = round(t(windows_out[windows_out$var == x & windows_out$window == "50k",1:11]),3))
  rownames(tmp)<-NULL
  colnames(tmp)[3]<-"Expected"
  
  # Fetch observed
  tmp$Observed <- t(read.table(paste0("outputs/50k/50k_",x,"Venn_overlaps.txt"),header=T)[,1:11])
  
  # Fetch FDRs...
  tmp_fdr <- t(read.table(paste0("outputs/50k/03/50k_ALL_byRegion_",iterations,"runs_pvals.txt"))[x,12:22])
  tmp$FDR <- tmp_fdr[,1]
  
  return(tmp)
})))

# Write supp table
write.csv(exp_obs_fdr,
          "tables/TableSX_exp_obs_fdr_50kb_windows_by_venn_segments.csv",
          quote = F,row.names = F)


##### Read back in data and calculate the Expected for 2 or more (ie. Not grouping specific) #####

# Set up PDF for output graphs
pdf("figs/Figure2_Exp_Obs_All_vars_all_windows.pdf",width=12,height=6)

windows_out<-data.frame(rbindlist(lapply(1:length(grouping_vector),function(x){
  
variables_out<-data.frame(rbindlist(lapply(1:length(variables_vector),function(y){
  
# Read in and calculate means
sim_dd<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_10000runs_FULLSIM.txt"),header=T)
sim_dd<-sim_dd[,1:11]

# Read in observed data and calculate p-values
obs_dd<-read.table(paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_Venn_overlaps.txt"))
colnames(obs_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")

# Subset for variable of interest
obs_dd<-obs_dd[y,1:11]

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

# FDR tests
variables_out$FDR<-p.adjust(variables_out$p.value,method = "fdr")
variables_out$FDR_asterisk<-rep(NA)

variables_out[variables_out$FDR <= 0.05,"FDR_asterisk"]<-"*"
variables_out[variables_out$FDR <= 0.01,"FDR_asterisk"]<-"**"
variables_out[variables_out$FDR <= 0.001,"FDR_asterisk"]<-"***"

# Modify for plotting
exp<-data.frame(variables_out[,!(colnames(variables_out) %in% c("Obs","UL"))],
                count_type="Expected")
colnames(exp)[1]<-"Count"
obs<-data.frame(variables_out[,!(colnames(variables_out) %in% c("Exp","UL"))],
                count_type="Observed")
colnames(obs)[1]<-"Count"

plot_dd<-rbind(exp,obs)

# Modify for plot
plot_dd$count_type<-as.character(plot_dd$count_type)
plot_dd[plot_dd$group == "Env" & plot_dd$count_type == "Observed","count_type"]<-"Observed Environment"
plot_dd[plot_dd$group == "Armour" & plot_dd$count_type == "Observed","count_type"]<-"Observed Armour"
plot_dd[plot_dd$group == "Shape" & plot_dd$count_type == "Observed","count_type"]<-"Observed Shape"
plot_dd[plot_dd$group == "Gill" & plot_dd$count_type == "Observed","count_type"]<-"Observed Gill Raker"

# Conf levels
conf<-unlist(list(c(variables_out$UL),rep(NA,19)))
# Labels
asterisk_labs<-c(rep(NA,19),plot_dd[plot_dd$count_type != "Expected","FDR_asterisk"])

# Plot
dodge<-position_dodge(.75)

# Set order
plot_dd$var<-gsub("_"," ",plot_dd$var)
plot_dd$var_F<-factor(plot_dd$var,levels=unique(plot_dd$var))
p1<-ggplot(plot_dd,aes(x=var_F,y=Count,fill=count_type))+
  geom_errorbar(data=plot_dd,aes(ymin=0.1,ymax=conf),width=.2,position=dodge) +
  geom_bar(position = dodge, stat="identity",width=0.75)+
  theme_bw()+
  theme(axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10),
        axis.text.x = element_text(angle=45,hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "right") +
  scale_fill_manual(values = c("grey70", "red4","forestgreen","blue4","gold2")) +
  labs(fill="")+
  ylim(0,(max(plot_dd$Count)+0.1*max(plot_dd$Count)))+
  xlab("Variable")+
  ylab(paste0("Count of Associated ",grouping_vector[x]," Windows Shared in >1 Radiation"))+
  geom_text(label = asterisk_labs,nudge_y = 1)

# Print to PDF
print(p1)

return(variables_out)
})))

dev.off()

##### Parallelism score by Rad pairing #####
pdf("figs/03_Parallelism_scores_allwindows.pdf",width=12,height=8)
windows_out<-lapply(1:length(grouping_vector),function(x){
  
  # Read in the pvals
tmp_dd<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_ALL_bySNP_",iterations,"runs_pvals.txt"))
tmp_dd<-tmp_dd[,c("AB","AI","AN","BI","BN","IN")]

  # Rearrange data
tmp_2<-data.frame(rbindlist(lapply(1:length(tmp_dd),function(y){
  col_dd<-data.frame(score=1-tmp_dd[,y],
                     var=rownames(tmp_dd),
                     rad=rep(colnames(tmp_dd)[y]),
                     type=c(rep("Env",7),rep("Pheno",12)))
  return(col_dd)
})))



  # Plot means and standard errors with violins
plot_dd<-tmp_2 %>% 
  group_by(rad,type) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())),score)

pd<-position_dodge(0.5)
score_plot<-ggplot(plot_dd, aes(rad,mean,shape=type))+
  geom_point(size=5,position = pd)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.4, position=pd)+
  theme_bw()+
  scale_x_discrete(labels=c("AB"="Alaska & BC",
                            "AI"="Alaska & Iceland",
                            "AN" ="Alaska & Scotland",
                            "BI"="BC & Iceland",
                            "BN"="BC & Scotland",
                            "IN"="Iceland & Scotland")) +
  theme(axis.title.y = element_text(size=24),
        axis.title.x  = element_blank(),
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=22,angle=45,hjust=1),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=18))+
  ylab("Parallelism Score")+
  labs(shape="")

print(score_plot)
})

dev.off()

##### Read back in data and calculate pooled expected and p-vals for env and pheno totals #####

env_exp_out<-data.frame(rbindlist(lapply(1:length(grouping_vector),function(x){
  
  # Env
  env_out<-lapply(1:7,function(y){
    
    # Read in and calculate means
    sim_dd<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_10000runs_FULLSIM.txt"),header=T)
    sim_dd$var<-variables_vector[y]
    
    # Read in observed data and calculate p-values
    obs_dd<-read.table(paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_Venn_overlaps.txt"))
    colnames(obs_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    # Subset for variable of interest
    obs_dd<-obs_dd[y,1:11]
    return(list(sim_dd[,c(1:11,16)],obs_dd))
  })
  
  # Calculate pooled expected
  env_all<-rbind(env_out[[1]][[1]],env_out[[2]][[1]],env_out[[3]][[1]],env_out[[4]][[1]],env_out[[5]][[1]],env_out[[6]][[1]],env_out[[7]][[1]])
  env_all_obs<-rbind(env_out[[1]][[2]],env_out[[2]][[2]],env_out[[3]][[2]],env_out[[4]][[2]],env_out[[5]][[2]],env_out[[6]][[2]],env_out[[7]][[2]])
  
  # All 4
  env4<-data.frame(Ca=env_all[env_all$var=="Ca","ABIN"],
                   Gyro=env_all[env_all$var=="Gyro","ABIN"],
                   Na=env_all[env_all$var=="Na","ABIN"],
                   pH=env_all[env_all$var=="pH","ABIN"],
                   Schisto=env_all[env_all$var=="Schisto","ABIN"],
                   Zn=env_all[env_all$var=="Zn","ABIN"],
                   Lake_Area=env_all[env_all$var=="Lake_Area","ABIN"])
  
  env_out_dd<-data.frame(Exp=mean(rowSums(env4)),
                     UL=quantile(rowSums(env4),probs = 0.95),
                     Obs=sum(env_all_obs$ABIN),
                     var="Env",
                     window=grouping_vector[x],
                     p.value=length(env4[rowSums(env4) >= sum(env_all_obs$ABIN),1])/iterations)
  
  # Any 3
  env3<-data.frame(Ca=rowSums(env_all[env_all$var=="Ca",2:5]),
                   Gyro=rowSums(env_all[env_all$var=="Gyro",2:5]),
                   Na=rowSums(env_all[env_all$var=="Na",2:5]),
                   pH=rowSums(env_all[env_all$var=="pH",2:5]),
                   Schisto=rowSums(env_all[env_all$var=="Schisto",2:5]),
                   Zn=rowSums(env_all[env_all$var=="Zn",2:5]),
                   Lake_Area=rowSums(env_all[env_all$var=="Lake_Area",2:5]))
  
  env_out_dd<-data.frame(Exp=mean(rowSums(env3)),
                         UL=quantile(rowSums(env3),probs = 0.95),
                         Obs=sum(unlist(list(env_all_obs[,2:5]))),
                         var="Env",
                         window=grouping_vector[x],
                         p.value=length(env3[rowSums(env3) >= sum(unlist(list(env_all_obs[,2:5]))),1])/iterations)
  
  # Any 2
  env2<-data.frame(Ca=rowSums(env_all[env_all$var=="Ca",6:11]))
  for(i in 1:7){
   # print(i)
    env2$new<-rowSums(env_all[env_all$var==variables_vector[i],6:11])
    colnames(env2)[i]<-variables_vector[i]
  }
  
  env_out_dd<-data.frame(Exp=mean(rowSums(env2)),
                           UL=quantile(rowSums(env2),probs = 0.95),
                           Obs=sum(rowSums(env_all_obs[,6:11])),
                           var="env",
                           window=grouping_vector[x],
                           p.value=length(env2[rowSums(env2) >= sum(rowSums(env_all_obs[,6:11])),1])/iterations)
  
  
  
})))

# Same for pheno
pheno_exp_out<-data.frame(rbindlist(lapply(1:length(grouping_vector),function(x){

  # Pheno
  pheno_out<-lapply(8:19,function(y){
    
    # Read in and calculate means
    sim_dd<-read.table(paste0("outputs/",grouping_vector[x],"/03/",grouping_vector[x],"_",variables_vector[y],"_byREGION_10000runs_FULLSIM.txt"),header=T)
    sim_dd$var<-variables_vector[y]
    
    # Read in observed data and calculate p-values
    obs_dd<-read.table(paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_Venn_overlaps.txt"))
    colnames(obs_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    # Subset for variable of interest
    obs_dd<-obs_dd[y,1:11]
    return(list(sim_dd[,c(1:11,16)],obs_dd))
  })
  
  # Calculate pooled expected
  pheno_all<-rbind(pheno_out[[1]][[1]],pheno_out[[2]][[1]],pheno_out[[3]][[1]],pheno_out[[4]][[1]],pheno_out[[5]][[1]],pheno_out[[6]][[1]],
                   pheno_out[[7]][[1]],pheno_out[[8]][[1]],pheno_out[[9]][[1]],pheno_out[[10]][[1]],pheno_out[[11]][[1]],pheno_out[[12]][[1]])
  pheno_all_obs<-rbind(pheno_out[[1]][[2]],pheno_out[[2]][[2]],pheno_out[[3]][[2]],pheno_out[[4]][[2]],pheno_out[[5]][[2]],pheno_out[[6]][[2]],
                       pheno_out[[7]][[2]],pheno_out[[8]][[2]],pheno_out[[9]][[2]],pheno_out[[10]][[2]],pheno_out[[11]][[2]],pheno_out[[12]][[2]])

  # All 4
  pheno4<-data.frame(Shape_PC1=pheno_all[pheno_all$var=="Shape_PC1","ABIN"])
  for(i in 8:19){
    pheno4$new<-pheno_all[pheno_all$var==variables_vector[i],"ABIN"]
    colnames(pheno4)[i-7]<-variables_vector[i]
  }
  
  pheno_out_dd<-data.frame(Exp=mean(rowSums(pheno4)),
                         UL=quantile(rowSums(pheno4),probs = 0.95),
                         Obs=sum(pheno_all_obs$ABIN),
                         var="pheno",
                         window=grouping_vector[x],
                         p.value=length(pheno4[rowSums(pheno4) >= sum(pheno_all_obs$ABIN),1])/iterations)
  
  # Any 3
  pheno3<-data.frame(Shape_PC1=rowSums(pheno_all[pheno_all$var=="Shape_PC1",2:5]))
  for(i in 8:19){
    pheno3$new<-rowSums(pheno_all[pheno_all$var==variables_vector[i],2:5])
    colnames(pheno3)[i-7]<-variables_vector[i]
  }
  
  pheno_out_dd<-data.frame(Exp=mean(rowSums(pheno3)),
                         UL=quantile(rowSums(pheno3),probs = 0.95),
                         Obs=sum(rowSums(pheno_all_obs[,2:5])),
                         var="pheno",
                         window=grouping_vector[x],
                         p.value=length(pheno3[rowSums(pheno3) >= sum(rowSums(pheno_all_obs[,2:5])),1])/iterations)
  
  # Any 2
  pheno2<-data.frame(Shape_PC1=rowSums(pheno_all[pheno_all$var=="Shape_PC1",6:11]))
  for(i in 8:19){
    pheno2$new<-rowSums(pheno_all[pheno_all$var==variables_vector[i],6:11])
    colnames(pheno2)[i-7]<-variables_vector[i]
  }
  
  pheno_out_dd<-data.frame(Exp=mean(rowSums(pheno2)),
                           UL=quantile(rowSums(pheno2),probs = 0.95),
                           Obs=sum(rowSums(pheno_all_obs[,6:11])),
                           var="pheno",
                           window=grouping_vector[x],
                           p.value=length(pheno2[rowSums(pheno2) >= sum(rowSums(pheno_all_obs[,6:11])),1])/iterations)
  
  
})))

# Get res
env_exp_out
pheno_exp_out
