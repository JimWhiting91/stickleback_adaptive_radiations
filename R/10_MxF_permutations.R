######################################################################################################################
# This script reads in the MxF outputs alongside those for Freshwater environmental variation to compare parallelism
######################################################################################################################

# Load packages
lib<-as.vector(c("PopGenome","dplyr","data.table","ggplot2","parallel","VennDiagram"))
lapply(lib,library,character.only=TRUE)

# Source functions
source("R/Genomic_Parallelism_Rcode_functions.R")
setwd("/Users/jw962/Google Drive/PNAS_Resubmission_Online-20180823")

# Analysis is grouped by window size
window_vec<-c("50k","75k","100k","200k")
window_size<-as.integer(c(50000,75000,100000,200000))

# Radiations
rad_vec<-c("Alaska","BC","Iceland","Scotland")
rad_vec2<-c("Al","Bc","Ic","Sc")

# Define the variables
variables_vector<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn",
                              "Shape_PC1","Shape_PC2","Shape_PC3",
                              "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                              "Gill_Raker_L","Gill_Raker_N"))
env_pheno_vec<-c(rep("Env",6),
                 rep("Shape",3),
                 rep("Armour",7),
                 rep("Gill",2))


# Loop analysis over windows
MxF_bars<-mclapply(1:length(window_vec),function(x){

##### First, what is the overlap of MxF outliers ######
  
  # Read in the data
  MxF<-lapply(1:4,function(y){
    
    # Data in
    tmp1<-read.table(paste0("data/outlier_beds/MxF/",rad_vec2[y],"_bayenv_XTX_windows_",window_size[x],".txt"),header=T)
    tmp2<-read.table(paste0("data/outlier_beds/MxF/",rad_vec2[y],"_bayenv_MxF_windows_",window_size[x],"_OUTLIERS.txt"),header=T)
    
    # Add Radiation column
    tmp1$rad<-rad_vec[y]
    tmp2$rad<-rad_vec[y]
    
    # Return
    return(list(tmp1,tmp2))
  })
  
  # Save as dataframes
  MxF_winds<-data.frame(rbindlist(lapply(1:4,function(y){return(MxF[[y]][[1]])})))
  MxF_outlier<-data.frame(rbindlist(lapply(1:4,function(y){return(MxF[[y]][[2]])})))
  MxF_outlier_list<-lapply(1:4,function(y){return(as.character(MxF[[y]][[2]]$window_id))})
  
  # What is overlap of outliers
  MxF_overlap<-calculate.overlap(MxF_outlier_list)
  names(MxF_overlap)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
  
  # Turn to table and export
  overlap_labels<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
  outlier_overlap<-na.omit(data.frame(rbindlist(lapply(1:length(MxF_overlap),function(venn){
    
    if(length(MxF_overlap[[venn]]) > 0){
    tmp=data.frame(window_id=MxF_overlap[[venn]],
                   overlap=overlap_labels[[venn]])
    return(tmp)
    } else {
      tmp=data.frame(window_id=NA,
                     overlap=overlap_labels[[venn]])
      return(tmp)    }
  }))))
  
  # Count
  MxF_overlap_counts<-list()
  for (i in 1:15){
    MxF_overlap_counts[[i]]<-length((as.vector(MxF_overlap[[i]])))
  }
  MxF_overlap_counts<-as.data.frame(MxF_overlap_counts)
  colnames(MxF_overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
  
  # Perform permutations
  iterations<-10000
  # MxF_perm<-data.frame(rbindlist(mclapply(1:iterations,function(y){
  #   
  #   # From each radiation element, draw a random cohort of windows equal to outlier N
  #   random_draw<-lapply(1:4,function(z){
  #     outlier_N<-nrow(MxF[[z]][[2]])
  #     random_windows<-as.character(sample(MxF[[z]][[1]]$window_id,outlier_N))
  #     return(random_windows)
  #   })
  #   
  #   # What's the overlap
  #   random_overlap<-calculate.overlap(random_draw)
  #   
  #   # Turn calculate.overlap into human readable output
  #   perm_overlap_counts<-list()
  #   for (i in seq(1,15,by=1)){
  #     perm_overlap_counts[[i]]<-length((as.vector(random_overlap[[i]])))
  #   }
  #   perm_overlap_counts<-as.data.frame(perm_overlap_counts)
  #   colnames(perm_overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
  #   
  #   return(perm_overlap_counts)
  # },mc.cores=detectCores())))
  # 
  # # Write output
  # write.table(MxF_perm,
  #             paste0("outputs/MxF/MxF_Permutation_iter",iterations,"_",window_vec[x],".txt"),row.names = F,quote = F,sep="\t")
  
  # Read in sims
  MxF_sims<-read.table(paste0("outputs/MxF/MxF_Permutation_iter",iterations,"_",window_vec[x],".txt"),header=T)[,1:11]
  
  # Calculate P-values for each level
  p_vals<-data.frame(overlap=t(MxF_overlap_counts)[1:11],
                     p=NA)
  for(i in 1:nrow(p_vals)){
    tmp<-MxF_sims[,i]
    tmp2<-tmp[tmp >= p_vals$overlap[i]]
    p_vals$p[i]<-length(tmp2)/iterations
  }
  p_vals$FDR<-p.adjust(p_vals$p,method = 'fdr')
  
  ##### Output outliers to table ######
  outlier_winds<-unique(MxF_outlier$window_id)
  MxF_outlier_out<-data.frame(rbindlist(lapply(1:length(outlier_winds),function(z){
    rads_to_merge<-MxF_outlier[MxF_outlier$window_id == outlier_winds[z],"rad"]
    rad_out<-paste(rads_to_merge, collapse=" & ")
    out<-data.frame(window_id=outlier_winds[z],
                    Radiations=rad_out,
                    N=length(rads_to_merge))
  })))
  
  # Re-order and write to tables
  MxF_outlier_out<-MxF_outlier_out[order(-MxF_outlier_out$N),]
  
  # Get outlier window measures
  tmp_chr<-strsplit(as.character(MxF_outlier_out$window_id),":")
  chrs<-unlist(lapply(tmp_chr,function(j){j[1]}))
  bps<-unlist(lapply(tmp_chr,function(j){j[2]}))
  tmp_bp<-strsplit(bps,"-")
  bp1<-unlist(lapply(tmp_bp,function(j){j[1]}))
  bp2<-unlist(lapply(tmp_bp,function(j){j[2]}))
  outlier_windows<-data.frame(chr=chrs,
                              BP1=bp1,
                              BP2=bp2)
  
  
  
  # Compare to bed file of previous Marine x Freshwater outliers
  MxF<-read.csv("data/marine_fresh_windows.csv")
  
  
  # Do they overlap?
  outlier_windows$MxF<-"No"
  for(i in 1:nrow(outlier_windows)){
    MxF_tmp<-MxF[as.character(MxF$group)==as.character(outlier_windows$chr[i]) &
                   MxF$Begin <= as.integer(as.character(outlier_windows$BP2[i])) &
                   MxF$End >= as.integer(as.character(outlier_windows$BP1[i])),"Source"]
    MxF_tmp<-unique(as.character(MxF_tmp))
    if(length(MxF_tmp) > 0){
      outlier_windows$MxF[i]<-paste(MxF_tmp,collapse="/")
    }
  }

  MxF_outlier_out$Previous_Studies<-outlier_windows$MxF
  
  
  write.csv(MxF_outlier_out,
              "tables/TableS6_MxF_outliers_50kb.csv",
              row.names = F)
  
  ####### Plot a new Figure 2 which includes the Marine Bar ##########
  
    sims_out<-lapply(1:length(variables_vector),function(y){
      
      # Read in and calculate means
      sim_dd<-read.table(paste0("outputs/",window_vec[x],"/03/",window_vec[x],"_",variables_vector[y],"_byREGION_10000runs_FULLSIM.txt"),header=T)
      sim_dd<-sim_dd[,1:11]
      
      return(sim_dd)
    })
    
    # Add marine x fresh sims
    sims_list<-list(MxF_sims)
    for(i in 1:length(sims_out)){
      sims_list[[i+1]]<-sims_out[[i]]
    }
      
      # Read in observed data and calculate p-values
      obs_dd<-read.table(paste0("outputs/",window_vec[x],"/",window_vec[x],"_Venn_overlaps.txt"))
      colnames(obs_dd)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
      obs_dd<-rbind(MxF_overlap_counts,obs_dd)
      rownames(obs_dd)[1]<-"M x F"
      
      # Subset for variable of interest
      obs_dd<-obs_dd[1:19,1:11]
      
      # Sum Exp and Obs with UL and LL across all
      marine_env_pheno_vec<-c("M x F",env_pheno_vec)
      out_dd<-data.frame(Exp=unlist(lapply(sims_list,function(sim){return(mean(rowSums(sim)))})),
                         UL=unlist(lapply(sims_list,function(sim){return(quantile(rowSums(sim),probs=0.95))})),
                         Obs=rowSums(obs_dd),
                         var=c("M x F",variables_vector),
                         window=rep(window_vec[x]),
                         group=marine_env_pheno_vec)
      
      # Fix the var columns
      out_dd$var<-rownames(out_dd)
      
      # Get p values
      p.value<-unlist(lapply(1:length(sims_list),function(z){
        
        return(length(rowSums(sims_list[[z]])[rowSums(sims_list[[z]]) >= sum(obs_dd[z,])])/iterations)
      }))
      out_dd$p<-p.value
      
    # FDR tests
      out_dd$FDR<-p.adjust(out_dd$p,method = "fdr")
      out_dd$FDR_asterisk<-rep(NA)
    
      out_dd[out_dd$FDR <= 0.05,"FDR_asterisk"]<-"*"
      out_dd[out_dd$FDR <= 0.01,"FDR_asterisk"]<-"**"
      out_dd[out_dd$FDR <= 0.001,"FDR_asterisk"]<-"***"
    
    # Modify for plotting
    exp<-data.frame(out_dd[,!(colnames(out_dd) %in% c("Obs","UL"))],
                    count_type="Expected")
    colnames(exp)[1]<-"Count"
    obs<-data.frame(out_dd[,!(colnames(out_dd) %in% c("Exp","UL"))],
                    count_type="Observed")
    colnames(obs)[1]<-"Count"
    
    plot_dd<-rbind(exp,obs)
    
    # Modify for plot
    plot_dd$count_type<-as.character(plot_dd$count_type)
    plot_dd[plot_dd$group == "M x F" & plot_dd$count_type == "Observed","count_type"]<-"Observed M x F"
    plot_dd[plot_dd$group == "Env" & plot_dd$count_type == "Observed","count_type"]<-"Observed Environment"
    plot_dd[plot_dd$group == "Armour" & plot_dd$count_type == "Observed","count_type"]<-"Observed Armour"
    plot_dd[plot_dd$group == "Shape" & plot_dd$count_type == "Observed","count_type"]<-"Observed Shape"
    plot_dd[plot_dd$group == "Gill" & plot_dd$count_type == "Observed","count_type"]<-"Observed Gill Raker"
    
    # Conf levels
    conf<-unlist(list(c(out_dd$UL),rep(NA,19)))
    # Labels
    asterisk_labs<-c(rep(NA,19),plot_dd[plot_dd$count_type != "Expected","FDR_asterisk"])
    
    
    # Plot
    dodge<-position_dodge(.75)
    p1<-ggplot(plot_dd,aes(x=var,y=Count,fill=count_type))+
      geom_errorbar(data=plot_dd,aes(ymin=0.1,ymax=conf),width=.2,position=dodge) +
      geom_bar(position = dodge, stat="identity",width=0.75)+
      theme_bw()+
      theme(axis.title.y = element_text(size=8),
            axis.text.x = element_text(angle=45,hjust = 1),
            panel.grid=element_blank(),
            legend.position = "none") +
      scale_fill_manual(values = c("grey70","red4","forestgreen","blue4","black","gold2")) +
      scale_x_discrete(limits=unique(plot_dd$var),
                       labels=unique(gsub("_"," ",plot_dd$var)))+
      labs(fill="")+
      ylim(0,(max(plot_dd$Count)+0.1*max(plot_dd$Count)))+
      xlab("")+
      ylab(paste0("N Parallel ",window_vec[x],"b Windows in >1 Radiation"))+
      geom_text(label = asterisk_labs,nudge_y = 1)

return(p1)
  # End of window vec
})

# Print figs  
pdf("figs/Figure2b_Marine_Fresh_comparsons.pdf",width=5,height=3)
for (i in 1:length(window_vec)){
  print(MxF_bars[[i]])
}
dev.off()

# # Print figs  
# pdf("resubmission_files_to_send/Figures/Figure2b_Marine_Fresh_comparsons_50kb.pdf",width=5,height=3)
# print(MxF_bars[[1]])
# dev.off()
