################################################################################################
# This script reads back in "adjacent" convergent windows and repeats permutations based on these
lib<-as.vector(c("data.table","VennDiagram","ggplot2","parallel","dplyr","GenomicRanges","tidyverse"))
lapply(lib,library,character.only=TRUE)

source("R/Genomic_Parallelism_Rcode_functions.R")

# Run over each variable
variables_vector<-c("Ca","Gyro","Na","pH","Schisto","Zn","Lake_Area",
                    "Shape_PC1","Shape_PC2","Shape_PC3",
                    "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                    "Gill_Raker_L","Gill_Raker_N")

rads<-c("Alaska","BC","Iceland","Scotland")

var_out<-data.frame(rbindlist(lapply(variables_vector,function(var){
  
  # Read in outliers
  regions<-read.table(paste0("outputs/50k/50k_",var,"_top_outliers_by_rad.txt"),header=T)
  
  # Get list of outliers and order numerically
  outliers<-regions[order(regions[,1],regions[,2]),c("radiation","window_id")]
  outliers$radiation<-as.character(outliers$radiation)
  outliers[outliers$radiation == "Nuist","radiation"]<-"Scotland"
  
  outliers_adj<-data.frame(rbindlist(lapply(rads,function(rad){
    
    # Sort for rad  
    outliers2<-outliers[outliers$radiation==rad,]
    
    # Fix
    outliers2$chr=unlist(tstrsplit(outliers2$window_id,":",keep=1))
    BPs<-unlist(tstrsplit(outliers2$window_id,":",keep=2))
    outliers2$BP1<-as.integer(unlist(tstrsplit(BPs,"-",keep=1)))
    
    # Sort
    outliers2<-outliers2[order(outliers2[,3],outliers2[,4]),]
    
    # Group adjacents and remove duplicates
    outliers2$adj_window_id<-condense_windows(outliers2$window_id)
    outliers2<-outliers2[!(duplicated(outliers2$adj_window_id)),]
    
    return(outliers2)
  })))
  
  # Remake chr intervals
  outliers_adj$chr=unlist(tstrsplit(outliers_adj$adj_window_id,":",keep=1))
  outliers_adj$BPs<-unlist(tstrsplit(outliers_adj$adj_window_id,":",keep=2))
  outliers_adj<-separate(outliers_adj,BPs,into=c("BP1","BP2"),"-")
  
  # Now get overlap with GenomicRanges
  # Make ranges
  rad_ranges<-lapply(rads,function(rad){
    tmp<-outliers_adj[outliers_adj$radiation==rad,]
    ranges<-makeGRangesFromDataFrame(tmp,
                                     ignore.strand = T,
                                     strand.field = "chr",
                                     start.field = "BP1",
                                     end.field = "BP2")
  })
  
  # Overlap
  comparisons_to_make<-combn(1:4,2)
  overlap_res<-data.frame(rbindlist(lapply(1:ncol(comparisons_to_make),function(i){
    rad1=rads[comparisons_to_make[1,i]]
    rad2=rads[comparisons_to_make[2,i]]
    overlaps<-countOverlaps(rad_ranges[[comparisons_to_make[1,i]]],rad_ranges[[comparisons_to_make[2,i]]],
                            type="any",
                            ignore.strand=T,minoverlap = 2L)
    overlaps[overlaps > 1]<-1
    out_winds<-outliers_adj[outliers_adj$radiation==rad1,]
    out_winds<-out_winds[which(overlaps!=0),]
    
    return(out_winds)
  })))
  
  # Final vals for permutations
  pairwise_overlap<-nrow(overlap_res)
  outlier_N<-sapply(rads,function(rad){return(nrow(outliers_adj[outliers_adj$radiation==rad,]))})
  
  # Permute
  iterations<-1000
  permute_out<-data.frame(rbindlist(mclapply(1:iterations,function(iter){
    
    # Draw each radiations windows
    wind_list<-lapply(rads,function(rad){
      winds<-read.table(paste0("data/outlier_beds/",rad,"_50k_windows_SNPcount.bed"))
      winds<-paste0(winds$V1,":",winds$V2,"-",winds$V3)
      rand_winds<-sample(winds,outlier_N[rad])
    })
    
    # Get overlap
    perm_overlap<-calculate.overlap(wind_list)
    names(perm_overlap)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    # Turn calculate.overlap into human readable output
    perm_overlap_counts<-list()
    for (i in seq(1,15,by=1)){
      perm_overlap_counts[[i]]<-length((as.vector(perm_overlap[[i]])))
    }
    perm_overlap_counts<-as.data.table(perm_overlap_counts)
    colnames(perm_overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
    
    return(perm_overlap_counts)
    # End of permute
  },mc.cores=detectCores()-1)))
  
  # Sum up the overlaps
  permute_out$overlapping<-rowSums(permute_out[,1:11])
  
  # Calculate pval
  var_p<-nrow(permute_out[permute_out$overlapping >= pairwise_overlap,])/iterations
  
  # Results tab
  out<-data.frame(variable=var,
                  pairwise_overlap=pairwise_overlap,
                  mean_permute=mean(permute_out$overlapping),
                  quantile=quantile(permute_out$overlapping,0.95),
                  p_val=var_p)
  
  return(out)
  # End of var
})))

# FDR
var_out$fdr<-p.adjust(var_out$p_val,method="fdr")

# Save to outputs
write.table(var_out,
            "outputs/Adjacent_window_permutations_results.txt",
            row.names = F,quote=F,sep="\t")




