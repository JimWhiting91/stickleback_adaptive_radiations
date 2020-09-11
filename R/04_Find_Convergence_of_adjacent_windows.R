################################################################
# This script reads back in the outlier regions and looks for convergence of any adjacent windows
#################################################################

# Load packages
lib<-as.vector(c("doParallel","dplyr","lsmeans","ggsignif","MASS","data.table","VennDiagram","ggplot2","parallel","effects","lme4","lmerTest"))
lapply(lib,library,character.only=TRUE)

# Source functions
source("R/Genomic_Parallelism_Rcode_functions.R")

# Run over all the window sizes
winds<-c("50k","75k","100k","200k")
wind_size<-c(50000,75000,100000,200000)

rads<-c("Alaska","BC","Iceland","Nuist")

# Run over variables
variables_vector<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn","Lake_Area",
                              "Shape_PC1","Shape_PC2","Shape_PC3",
                              "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                              "Gill_Raker_L","Gill_Raker_N"))

# Make a "types" vector to explore previous evidences
qtl_types<-c(rep("Env",7),
             rep("Shape",3),
             rep("Armour",7),
             rep("Gills",2))

# Labels
labels<-c("All",
          "Alaska & BC & Iceland",
          "Alaska & BC & Scotland",
          "Alaska & Iceland & Scotland",
          "BC & Iceland & Scotland",
          "Alaska & BC",
          "Alaska & Iceland",
          "Alaska & Scotland",
          "BC & Iceland",
          "BC & Scotland",
          "Iceland & Scotland")

# Windows
winds_out<-lapply(1:length(winds),function(x){
  
# Variable
vars_out<-data.frame(rbindlist(lapply(1:length(variables_vector),function(y){
  
# Read in the list of outliers
dd<-read.table(paste0("outputs/",winds[x],"/",winds[x],"_",variables_vector[y],"_top_outliers_by_rad.txt"),header=T)

# Get list of outliers
outliers<-dd[order(dd[,1],dd[,2]),c("radiation","window_id")]
outliers$radiation<-as.character(outliers$radiation)
outliers[outliers$radiation == "Nuist","radiation"]<-"Scotland"


# Concat to only retain unique
outliers2<-data.frame(outliers %>%
                   group_by(window_id) %>%
                   summarise(radiation=paste(unique(radiation),collapse = ' & ')))

# Fix
outliers2$chr=unlist(tstrsplit(outliers2$window_id,":",keep=1))
BPs<-unlist(tstrsplit(outliers2$window_id,":",keep=2))
outliers2$BP1<-as.integer(unlist(tstrsplit(BPs,"-",keep=1)))

# Sort
outliers2<-outliers2[order(outliers2[,3],outliers2[,4]),]

# Group adjacents
outliers2$adj_window_id<-condense_windows(outliers2$window_id)

# Fix the radiations
outliers2_2<-data.frame(rbindlist(lapply(1:length(outliers2$adj_window_id),function(x){
  tmp<-outliers2[outliers2$adj_window_id == outliers2$adj_window_id[x],]
  rads<-unique(unlist(strsplit(tmp$radiation," & ")))
  rads<-rads[order(rads)]
  
  out<-data.frame(window_id = tmp$adj_window_id[1],
                  rads = paste(rads,collapse = ' & '))
  return(out)
})))

# Remove dupes
outliers3<-outliers2_2[!duplicated(outliers2_2),]

# Remove singletons
'%!in%' <- function(x,y)!('%in%'(x,y))
outliers3<-outliers3[outliers3$rads %!in% c("Alaska","BC","Iceland","Scotland"),]

# Get outlier window measures
tmp_chr<-strsplit(as.character(outliers3$window_id),":")
chrs<-unlist(lapply(tmp_chr,function(j){j[1]}))
bps<-unlist(lapply(tmp_chr,function(j){j[2]}))
tmp_bp<-strsplit(bps,"-")
bp1<-unlist(lapply(tmp_bp,function(j){j[1]}))
bp2<-unlist(lapply(tmp_bp,function(j){j[2]}))
outlier_windows<-data.frame(chr=chrs,
                            BP1=bp1,
                            BP2=bp2)


# Compare to bed file of previous Marine x Freshwater outliers
MxF<-read.csv("data/marine_fresh_windows_v2.csv")
MxF$group<-gsub("chr","group",MxF$group)

# Only use the right type of qtls
if(qtl_types[y] != "Env"){
  MxF<-MxF[MxF$Type == qtl_types[y],]
}

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

outliers3$Previous_Studies<-outlier_windows$MxF

# Export
dir.create(paste0("outputs/",winds[x],"/08/"),showWarnings = F)
write.table(outliers3,
            paste0("outputs/",winds[x],"/08/",winds[x],"_",variables_vector[y],"_Adjacent_Convergence.txt"),
            row.names = F, quote = F,sep="\t")

outliers3$var<-variables_vector[y]
return(outliers3)
})))

return(vars_out)
})

write.csv(winds_out[[1]],
          "tables/TableSX_with_QTL_comparisons_50kb.csv",
          row.names = F)
