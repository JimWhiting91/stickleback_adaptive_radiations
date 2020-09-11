#####################################################################################
# This script looks at shared SNPs within windows and assesses whether overlapping outlier windows are enriched for shared SNPs...
#####################################################################################

# We need these...
lib<-c("multcomp","data.table","VennDiagram","ggplot2","parallel","dplyr","vcfR","UpSetR")
lapply(lib,library,character.only=TRUE)

# Functions
human_readable_overlap<-function(overlap=NULL){
  perm_overlap<-overlap
names(perm_overlap)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")

# Turn calculate.overlap into human readable output
perm_overlap_counts<-list()
for (i in seq(1,15,by=1)){
  perm_overlap_counts[[i]]<-length((as.vector(perm_overlap[[i]])))
}
perm_overlap_counts<-as.data.table(perm_overlap_counts)
colnames(perm_overlap_counts)<-c("ABIN","ABI","ABN","AIN","BIN","AB","AI","AN","BI","BN","IN","A","B","I","N")
return(list(perm_overlap,perm_overlap_counts))
}

# Window size
window_size<-50000

# Now read in lists of all SNPs per radiation
rads<-c("Alaska","BC","Iceland","Scotland")
rad_shorts<-c("Al","Bc","Ic","Sc")
rad_snps<-mclapply(rads,function(rad){
  
  # Read in the snps
  tmp<-read.table(paste0("data/raw_bayenv_outputs/",rad,"_CHRBP"))
  
  # Add ID
  snp_id<-paste0(tmp$V2,"_",tmp$V3)
  return(snp_id)
},mc.cores=4)

# Now compare lists of SNPs for overlap among all
snp_ids<-rad_snps
overlap<-calculate.overlap(snp_ids)
overlap2<-human_readable_overlap(overlap)

# Counts
overlap2[[2]]

# Plot
names(snp_ids)<-rads
upset(fromList(snp_ids),keep.order = F)

# Plot to pdf
pdf("figs/FigureSX_Shared_Polymorphisms_within_radiations.pdf",width=12,height=5)
upset(fromList(snp_ids),keep.order = F)
dev.off()

###### By Windows
all_snps<-data.frame(rbindlist(lapply(rads,function(rad){
  out <- data.frame(snp=snp_ids[[rad]])
  out$snp <- gsub("scaffold_","scaffold-",out$snp)
  out <- tidyr::separate(out,snp,"_",into=c("chr","pos"))
  out$rad <- rad
  return(out)
})))
all_snps$pos<-as.integer(all_snps$pos)

# Make windows
chrs<-unique(all_snps$chr)
windows<-data.frame(rbindlist(lapply(chrs,function(chr){
  winds<-seq(0,max(all_snps[as.character(all_snps$chr) == chr,"pos"]),window_size)
  winds2<-winds+window_size
  bed<-data.frame(chr=chr,
                  start=winds,
                  end=winds2)
})))

# Run through the windows and do the overlap of SNPs
snp_sharing_by_window<-data.frame(rbindlist(mclapply(1:nrow(windows),function(x){
  print(x)
  # Get window
  window<-windows[x,]
  
  # Fetch snps
  tmp_snps<-lapply(rads,function(rad){
    tmp<-all_snps[all_snps$rad==rad,]
    
    tmp2<-tmp[tmp$chr==window$chr[1] &
                tmp$pos > window$start[1],]
    tmp2 <- tmp2[tmp2$pos < window$end[1],]
    
    if(nrow(tmp2) > 0){
    tmp2$snp_id<-paste0(tmp2$chr,"_",tmp2$pos)
    return(tmp2$snp_id)
    } else {
      return(paste0(rad,"_empty"))  
    }
})
  names(tmp_snps)<-c("A","B","I","N")
  overlap<-calculate.overlap(tmp_snps)
  overlap2<-data.frame(human_readable_overlap(overlap)[[2]])
  
  # For each calculate the expected and return fold-enrichment...
  unique_snps <- length(unique(unlist(tmp_snps)))
  overlap_count <- sum(overlap2[1:11])
  
  alaska_count <- length(tmp_snps[[1]])
  bc_count <- length(tmp_snps[[2]])
  iceland_count <- length(tmp_snps[[3]])
  scotland_count <- length(tmp_snps[[4]])

  alaska_odds <- alaska_count/unique_snps
  bc_odds <- bc_count/unique_snps
  iceland_odds <- iceland_count/unique_snps
  scotland_odds <- scotland_count/unique_snps
  
  # Probabilty of "any" overlap
  #expected = (1-((1-alaska_odds) * (1-bc_odds) * (1-iceland_odds) * (1-scotland_odds)))*unique_snps
  
  # Fold enrichment...
  #enrich = (overlap_count+1)/expected
  
  overlap2$ABIN <- overlap2$ABIN / (alaska_odds * bc_odds * iceland_odds * scotland_odds * unique_snps)

  overlap2$ABI <- overlap2$ABI / (alaska_odds * bc_odds * iceland_odds * (1-scotland_odds) * unique_snps)
  overlap2$ABN <- overlap2$ABN / (alaska_odds * bc_odds * (1-iceland_odds) * scotland_odds * unique_snps)
  overlap2$AIN <- overlap2$AIN / (alaska_odds * (1-bc_odds) * iceland_odds * scotland_odds * unique_snps)
  overlap2$BIN <- overlap2$BIN / ((1-alaska_odds) * bc_odds * iceland_odds * scotland_odds * unique_snps)

  overlap2$AB <- overlap2$AB / (alaska_odds * bc_odds * (1-iceland_odds) * (1-scotland_odds) * unique_snps)
  overlap2$AI <- overlap2$AI / (alaska_odds * (1-bc_odds) * iceland_odds * (1-scotland_odds) * unique_snps)
  overlap2$AN <- overlap2$AN / (alaska_odds * (1-bc_odds) * (1-iceland_odds) * scotland_odds * unique_snps)
  overlap2$BI <- overlap2$BI / ((1-alaska_odds) * bc_odds * iceland_odds * (1-scotland_odds) * unique_snps)
  overlap2$BN <- overlap2$BN / ((1-alaska_odds) * bc_odds * (1-iceland_odds) * scotland_odds * unique_snps)
  overlap2$IN <- overlap2$IN / ((1-alaska_odds) * (1-bc_odds) * iceland_odds * scotland_odds * unique_snps)
  
  # Replace NAs with 0s...
  overlap2[1,is.na(overlap2)] <- 0
  
  #out<-cbind(window,overlap2[,1:11])
  out <- data.frame(chr=windows[x,1],
                    start=as.integer(windows[x,2]),
                    end=as.integer(windows[x,3]),
                    overlap2[,1:11])
return(out) 

},mc.cores=detectCores()-1)))


##################################################
# Are overlapping outliers enriched for shared SNPs?
# Set up vectors
rads<-c("Alaska","BC","Iceland","Scotland")
variables_vector<-c("Ca","Gyro","Na","pH","Schisto","Zn","Lake_Area",
                    "Shape_PC1","Shape_PC2","Shape_PC3",
                    "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                    "Gill_Raker_L","Gill_Raker_N")

# Read in all windows and outliers
windows<-data.frame(rbindlist(lapply(1:4, function(x){
  
  # Get windows
  windows<-read.table(paste0("data/outlier_beds/",rads[x],"_50k_windows_SNPcount.bed"))
  windows<-windows[,-4]
  colnames(windows)<-c("chr","BP1","BP2")
  windows$window_id<-paste0(windows$chr,":",windows$BP1,"-",windows$BP2)
  return(windows)
})))

# Remove dups
windows<-windows[!(duplicated(windows$window_id)),]

# Get outliers
outliers<-lapply(1:length(variables_vector),function(y){
  tmp<-read.table(paste0("outputs/50k/50k_",variables_vector[y],"_top_outliers_by_rad.txt"),header=T)
  outliers<-paste0(tmp$chr,":",tmp$BP1,"-",tmp$BP2)
  overlap<-outliers[duplicated(outliers)]
  return(list(outliers,overlap))
})

# Extract
outliers2<-unique(unlist(lapply(outliers,function(input){return(input[[1]])})))
overlap2<-unique(unlist(lapply(outliers,function(input){return(input[[2]])})))

# Mark outliers
windows$outliers<-"Neutral"
windows[windows$window_id %in% outliers2,"outliers"]<-"Outlier"
windows[windows$window_id %in% overlap2,"outliers"]<-"Overlap"

# Get overlapping SNP counts
snp_sharing_by_window$window_id<-paste0(snp_sharing_by_window$chr,":",as.integer(snp_sharing_by_window$start),"-",as.integer(snp_sharing_by_window$end))
windows$enrich<-sapply(1:nrow(windows),function(x){

   tmp <- snp_sharing_by_window[snp_sharing_by_window$window_id == windows$window_id[x],]
   if(nrow(tmp) > 0){
    tmp2 <- tmp[,4:14]
    tmp3 <- tmp2[tmp2 != 0]
   enrich<-as.numeric(sum(log10(tmp3)))
   return(enrich)
   } else {
     return(NA)
   }
})

# Now summarise by category
na.omit(windows) %>%
  group_by(outliers) %>%
  dplyr::summarise(enrich_mean=mean(enrich),
            median=median(enrich),
            sd=sd(enrich),
            N=n())

# Assess significance
windows$outliers_F<-factor(windows$outliers,levels = unique(windows$outliers))
clean_windows<-na.omit(windows)
clean_windows$log_enrich <- log10(clean_windows$enrich+1)

# Visualise
ggplot(clean_windows,aes(x=outliers_F,y=enrich))+geom_boxplot()

glm1<-glm(enrich~outliers_F,data=clean_windows,family="gaussian")
drop1(glm1,test="Chisq")
summary(glht(glm1, mcp(outliers_F="Tukey")))

ggplot(clean_windows,aes(x=outliers,y=log(log_enrich)))+
 # geom_violin()+
  geom_boxplot()

# ##################################################
# # Are overlapping outliers enriched for shared SNPs in specific variables
# rads<-c("Alaska","BC","Iceland","Scotland")
# rad_letters<-c("A","B","I","N")
# variables_vector<-c("Ca","Gyro","Na","pH","Schisto","Zn","Lake_Area",
#                     "Shape_PC1","Shape_PC2","Shape_PC3",
#                     "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
#                     "Gill_Raker_L","Gill_Raker_N")
# 
# # Read in all windows and outliers
# windows<-data.frame(rbindlist(lapply(1:4, function(x){
#   
#   # Get windows
#   windows<-read.table(paste0("data/outlier_beds/",rads[x],"_50k_windows_SNPcount.bed"))
#   windows<-windows[,-4]
#   colnames(windows)<-c("chr","BP1","BP2")
#   windows$window_id<-paste0(windows$chr,":",windows$BP1,"-",windows$BP2)
#   return(windows)
# })))
# 
# # Remove dups
# windows<-windows[!(duplicated(windows$window_id)),]
# 
# # Get outliers
# shared_snps_outliers<-lapply(1:length(variables_vector),function(y){
#   tmp<-read.table(paste0("outputs/50k/50k_",variables_vector[y],"_top_outliers_by_rad.txt"),header=T)
#   # outliers<-paste0(tmp$chr,":",tmp$BP1,"-",tmp$BP2)
#   # overlap<-outliers[duplicated(outliers)]
#   # return(list(outliers,overlap))
#   
#   # Separate for list of vectors per rad
#   outlier_vecs<-lapply(rads,function(rad){return(tmp[tmp$radiation == rad,"window_id"])})
#   names(outlier_vecs)<-rads
#   
#   # Separate into outliers and overlapping
#   outliers<-as.character(unlist(outlier_vecs))
#   overlapping<-outliers[duplicated(outliers)]
#   outliers<-outliers[!(outliers %in% overlapping)]
#   
#   # Now separate out into 
#   neutral_shared<-snp_sharing_by_window[!(snp_sharing_by_window$window_id %in% c(overlapping,outliers)),]
#   outlier_shared<-snp_sharing_by_window[snp_sharing_by_window$window_id %in% outliers,]
#   overlapping_shared<-snp_sharing_by_window[snp_sharing_by_window$window_id %in% overlapping,]
#   
#   # Calculate mean and sd
#   out<-data.frame(window=c("Neutral","Outlier","Overlapping"),
#                   avg_shared=c(mean(rowSums(neutral_shared[,4:14])),
#                                mean(rowSums(outlier_shared[,4:14])),
#                                mean(rowSums(overlapping_shared[,4:14]))),
#                   N_windows=c(length(rowSums(neutral_shared[,4:14])),
#                               length(rowSums(outlier_shared[,4:14])),
#                               length(rowSums(overlapping_shared[,4:14]))),
#                   var=variables_vector[y])
#   
#   # Compare with wilcox?
#   if(nrow(overlapping_shared) > 0){
#   test<-wilcox.test(rowSums(outlier_shared[,4:14]),
#                     rowSums(overlapping_shared[,4:14]))
#   test.out<-data.frame(var=variables_vector[y],
#                        stat=test$statistic,
#                        p_val=test$p.value)
#   } else {
#     test.out<-data.frame(var=variables_vector[y],
#                          stat=NA,
#                          p_val=NA)
#   }
# 
# return(list(out,test.out))
# })
# 
# # Extract
# shared_snps_by_var<-data.frame(rbindlist(lapply(shared_snps_outliers,function(x){return(x[[1]])})))
# shared_snps_tests<-data.frame(rbindlist(lapply(shared_snps_outliers,function(x){return(x[[2]])})))
# 
# 
# # Assess
# ggplot(shared_snps_outliers,aes(x=shared_sites,y=overlapping))+
#   geom_point()

