######################################################################################################################
# This script reads in linkage stats for each Radiation and uses the information to build windows of SNPs in linkage
######################################################################################################################

# Load packages
lib<-as.vector(c("data.table","PopGenome","dplyr","data.table","ggplot2","parallel"))
lapply(lib,library,character.only=TRUE)

# Source functions
source("R/Genomic_Parallelism_Rcode_functions.R")
setwd("/Users/jw962/Google Drive/PNAS_Resubmission_Online-20180823")

# Need to calculate over each radiation independently
rads<-c("Alaska","BC","Iceland","Scotland")

# IMPORTANT - HERE WE SET THE WINDOW SIZE AT 0.1 cM
window_size<-0.1

##################################################################################
# Calculate genetic map distances using Roesti 2013 Recombination file
##################################################################################

# Read in rates from Roesti
recomb<-read.table("data/roesti_BROADS1_recomb_rates.txt",header=T)
# Remove group 19 and remove incorrect rows
recomb<-recomb[recomb$lg != 19,]
recomb<-recomb[recomb$pos2 > recomb$pos1,]


# Run over rads
rad_out<-mclapply(1:4,function(x){
  
  # Read in the relevant SNP positions
  CHRBP<-read.table(paste0("data/raw_bayenv_outputs/",rads[x],"_CHRBP"),header=F)
  chrs<-paste0("group",as.roman(c(1:18,20:21)))
  LG_names<-c(1:18,20:21)
  
  # Run over LGs
  chr_out<-lapply(1:length(chrs),function(y){
    print(chrs[y])
    # Subset
    CHR_sub<-CHRBP[CHRBP$V2==chrs[y],]
    recomb_sub<-recomb[recomb$lg==LG_names[y],]
    recomb_sub<-recomb_sub[order(recomb_sub$pos1),]
    CHR_sub<-CHR_sub[order(CHR_sub$V3),]
    colnames(CHR_sub)<-c("N","CHR","BP")
    
    # We can only calculate for where we have estimates
    CHR_sub2<-data.frame(rbindlist(lapply(1:nrow(CHR_sub),function(n){
      out<-recomb_sub[recomb_sub$pos1 <= CHR_sub$BP[n] & recomb_sub$pos2 >= CHR_sub$BP[n],]
      if(nrow(out) > 0){
        return(CHR_sub[n,])
      }
    })))
    
    # Attach the starting point
    starter<-data.frame(N=1,
                        CHR=CHR_sub2$CHR[1],
                        BP=recomb_sub$pos1[1])
    CHR_sub2<-rbind(starter,CHR_sub2)
    
    # If they exist, remove dups
    CHR_sub2<-CHR_sub2[!(duplicated(CHR_sub2$BP)),]
    
    # To approximate cM, we want to know recombination rate between 2 SNPs in cM/Mb and multiply by the distance
    CHR_sub2$cM<-0
    for(i in 2:nrow(CHR_sub2)){
      # Get relevant recombination rates
      recomb_rates<-recomb_sub[recomb_sub$pos2 >= CHR_sub2$BP[i-1] &
                                 recomb_sub$pos1 <= CHR_sub2$BP[i],]
      recomb_rates$pos1[1]<-CHR_sub2$BP[i-1]
      recomb_rates$pos2[nrow(recomb_rates)]<-CHR_sub2$BP[i]
      
      # Remember this is Mb not bp so divide by 1,000,000
      recomb_rates$cM<-((recomb_rates$pos2-recomb_rates$pos1)/1000000)*recomb_rates$rate
      
      CHR_sub2$cM[i]<-CHR_sub2$cM[i-1]+sum(recomb_rates$cM)
    }
    
    # Add Rad IDs
    CHR_sub2$RAD<-rads[x]
    
    # Can use rad_out to generate per radiation beds with SNP count per window
    max_window<-max(CHR_sub2$cM)
    wind_start<-seq(0,max_window,window_size)
    wind_end<-wind_start+window_size
    window_dd<-data.frame(CHR=CHR_sub2$CHR[1],
                          RAD=rads[x],
                          wind_start=wind_start,
                          wind_end=wind_end)
    window_dd$window_id<-paste0(window_dd$CHR,":",as.numeric(window_dd$wind_start),"-",as.numeric(window_dd$wind_end))
    
    # How many SNPs are in each window
    window_dd$SNP_Count<-0
    for(i in 1:nrow(window_dd)){
      window_dd$SNP_Count[i]<-nrow(CHR_sub2[(CHR_sub2$cM)*1e+10 >= (window_dd$wind_start[i])*1e+10 & 
                                         (CHR_sub2$cM)*1e+10 <= (window_dd$wind_end[i])*1e+10,])
    }
    
    # Remove windows with no SNPs because they're useless
    covered_windows<-window_dd[window_dd$SNP_Count > 0,]
    
    # Return
    return(list(CHR_sub2,covered_windows))
  })

  # Extract each
  SNP_pos<-data.frame(rbindlist(lapply(chr_out,function(chr){return(chr[[1]])})))
  window_out<-data.frame(rbindlist(lapply(chr_out,function(chr){return(chr[[2]])})))

# Return  
return(list(SNP_pos,window_out))
},mc.cores=detectCores()-1)

# We now have a list for each radiation that includes a table of all SNP positions and covered windows
# We now need to fill this with the counts of outlier SNPs in each one for each variable
# We are aiming to match the output of bed files constructed in bedtools

# First we create a vector of variable names
variable_names<-c("Ca","Gyro","Na","pH","Schisto","Zn",
                            "Shape_PC1","Shape_PC2","Shape_PC3",
                            "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                            "Gill_Raker_L","Gill_Raker_N")
# Now we list the corresponding column identifiers
columns_ID_vector<-c(11,14,8,2,17,5,
                               20,23,26,
                               47,50,53,56,59,62,65,
                               71,74)

# Run over each radiation
mclapply(1:length(rads),function(x){
  
  # Get new tables
  cM_pos<-rad_out[[x]][[1]]
  LD_windows<-rad_out[[x]][[2]]
  dd<-data.frame(read.table(paste0("data/raw_bayenv_outputs/",rads[x],"_EnvPheno_sorted_AVG_v2.txt")))
  # Remove column of 0s
  dd<-dd[,-2]
  
  # Read in the relevant SNP positions
  dd_SNP<-data.frame(fread(paste0("data/",rads[x],"_CHRBP"),header=F))
  
  # Prior to looking at variables, we want to assign cM values
  cM_pos$SNP_ID<-paste0(cM_pos$CHR,":",cM_pos$BP)
  dd$SNP_ID<-paste0(dd_SNP$V2,":",dd_SNP$V3)
  
  # Get cM vector
  cM_out<-sapply(1:nrow(dd),function(SNP){
    cM<-cM_pos[cM_pos$SNP_ID == dd$SNP_ID[SNP],"cM"]
    if(length(cM)==0){
      return(NA)
    } else {
      return(cM)
    }
  })
  dd$cM<-cM_out
  
  # Run over each variable
  lapply(1:length(variable_names),function(i){
    
    # Make the beds
    # Get Chrom BP1 and BP2
    dd_var<-data.frame(cbind(as.character(dd_SNP$V2),dd_SNP$V3,dd_SNP$V3+1))
    # Get Index
    dd_var$INDEX<-seq(1:length(dd_var$X1))
    # Define output columns to take
    col1<-columns_ID_vector[i]
    col2<-col1+2
    dd_var<-cbind(dd_var,dd[,col1:col2])
    # Tidy
    colnames(dd_var)<-c("Chrom","BP1","BP2","INDEX","BF","Spear","Pear")
    dd_var$cM<-dd$cM
    
    # Outliers defined as logBF > 1.5 and top 5% of ABSOLUTE (direction of effect doesn't matter) spearman quantile
    dd_var_outliers<-dd_var[log10(dd_var$BF) > 1.5 & abs(dd_var$Spear) > quantile(abs(dd_var$Spear),probs = 0.95),]

    # Using this info, we can assign outlier counts to LD_windows
    LD_windows$outlier_N<-0
    for(j in 1:nrow(LD_windows)){
    outlier_sub<-dd_var_outliers[dd_var_outliers$Chrom == LD_windows$CHR[j],]
    LD_windows$outlier_N[j]<-nrow(na.omit(outlier_sub[as.numeric(outlier_sub$cM) >= (LD_windows$wind_start[j]) &
                                               as.numeric(outlier_sub$cM) <= (LD_windows$wind_end[j]),]))
    }
    
    outlier_windows<-LD_windows[LD_windows$outlier_N > 0,]
      
    # Now export to meet formatting standards of downstream analyses
    out_mat<-matrix(c(as.character(outlier_windows$CHR),
                      outlier_windows$wind_start,
                      outlier_windows$wind_end,
                      outlier_windows$outlier_N),
                      ncol=4)
    write.table(out_mat,
                paste0("data/outlier_beds/cM/cM_windows_",rads[x],"_",variable_names[i],"_cleaned_outliers.bed_sorted.bed"),
                quote = F,row.names = F,sep="\t",col.names = F)
    
  })
  
  # Output radiation SNP counts in agreeable formatting for downstream
  out_mat<-matrix(c(as.character(LD_windows$CHR),
                    LD_windows$wind_start,
                    LD_windows$wind_end,
                    LD_windows$SNP_Count),ncol=4)
  write.table(out_mat,
              paste0("data/outlier_beds/",rads[x],"_cM_windows_SNPcount.bed"),
              quote = F,row.names = F,sep="\t",col.names = F)
},mc.cores=detectCores()-1)
  
