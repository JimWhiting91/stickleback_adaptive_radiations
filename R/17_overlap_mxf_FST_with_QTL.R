#################################################################################
# This script compares conventional Fst outliers for Marine x Freshwater comparisons with QTLs

# We need these
lib<-c("PopGenome","ggplot2","data.table","GenomicRanges","dplyr","tidyverse")
lapply(lib,library,character.only=T)

# Fetch the Fst results
rads<-c("Alaska","BC","Iceland","Scotland")
dd<-data.frame(rbindlist(lapply(1:4,function(x){
  tmp<-fread(paste0("outputs/",rads[x],"_MxF_Fst_50kb_windows.txt"))
  return(tmp)
})))

# Return only the outliers
pops<-unique(dd$pop)
outliers<-data.frame(rbindlist(lapply(pops,function(x){
  tmp<-dd[dd$pop == x,]
  out<-tmp[tmp$FST > quantile(tmp$FST, 0.95),]
  
  # We also need to breakup the window ID
  out<-separate(out,window_id,sep=c(":"),into=c("chr2","window"))
  out<-separate(out,window,sep=c("-"),into=c("start","end"))
  return(out)    
})))

# Get the QTLs
qtls<-read.csv("data/marine_fresh_windows_v2.csv")
colnames(qtls)<-c("chr","start","end","source","type")
qtls$chr<-gsub("chr","group",qtls$chr)

# Tidy up bad QTLs
for(i in 1:nrow(qtls)){
  if(qtls$start[i] > qtls$end[i]){
    tmp1<-qtls$start[i] 
    tmp2<-qtls$end[i] 
    qtls$start[i]<-tmp2
    qtls$end[i]<-tmp1
  }
}

qtls$qtl_size<-qtls$end-qtls$start
qtls<-qtls[qtls$qtl_size < 5000000,]
qtl_range<-makeGRangesFromDataFrame(qtls,
                                    ignore.strand = T,
                                    strand.field = "chr",
                                    start.field = "start",
                                    end.field = "end")

# # Make ranges for comparison
# rad_ranges<-lapply(rads,function(rad){
#   tmp<-outliers[outliers$radiation==rad,]
#   ranges<-makeGRangesFromDataFrame(tmp,
#                                    ignore.strand = T,
#                                    strand.field = "chr",
#                                    start.field = "start",
#                                    end.field = "end")
# })
# 
# # Compare each with qtls
# range_overlaps<-lapply(1:4,function(x){
# overlaps<-countOverlaps(qtl_range,rad_ranges[[x]],
#                         type="any",
#                         ignore.strand=T,minoverlap = 2L)
# 
# # For each, retain only those with at least 5 overlaps
# to_keep<-qtls_mxf[which(overlaps > 5),]
# to_keep$overlapping_N<-overlaps[which(overlaps > 5)]
# return(to_keep)
# })
# 
# # How many of each are covered?
# for(i in 1:4){
# print(nrow(range_overlaps[[i]])/nrow(outliers[outliers$radiation==rads[i],]))
# }

########## Repeat for MxF associated #############
# Fetch data
rad_shorts<-c("Al","Bc","Ic","Sc")

MxF_associated<-data.frame(rbindlist(lapply(1:4,function(x){
  tmp<-read.table(paste0("data/outlier_beds/MxF/",rad_shorts[x],"_bayenv_MxF_windows_50000_OUTLIERS.txt"),header=T)
  tmp$rad<-rads[x]
  return(tmp)
})))

MxF_associated<-MxF_associated[MxF_associated$chr %in% paste0("group",as.roman(1:21)),]

# Make ranges for comparison
rad_ranges<-lapply(rads,function(rad){
  tmp<-MxF_associated[MxF_associated$rad==rad,]
  ranges<-makeGRangesFromDataFrame(tmp,
                                   ignore.strand = T,
                                   strand.field = "chr",
                                   start.field = "BP1",
                                   end.field = "BP2")
})

# Compare each with qtls
range_overlaps<-lapply(1:4,function(x){
  tmp_outliers<-MxF_associated[MxF_associated$rad == rads[x],]
  overlaps<-countOverlaps(rad_ranges[[x]],qtl_range,
                          type="any",
                          ignore.strand=T,minoverlap = 2L)
  
  # For each, retain only those with at least 5 overlaps
  to_keep<-tmp_outliers[which(overlaps > 0),]
  to_keep$overlapping_N<-overlaps[which(overlaps > 0)]
  return(to_keep)
})

# How many of each are covered?
for(i in 1:4){
  print(nrow(range_overlaps[[i]])/nrow(MxF_associated[MxF_associated$rad==rads[i],]))
}

