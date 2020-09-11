######################################################
# Take raw BF beds, and intersect with SNP Count beds

# Load packages
lib<-as.vector(c("data.table","PopGenome","dplyr","data.table","ggplot2","parallel"))
lapply(lib,library,character.only=TRUE)

# Define the grouping
grouping_vector<-c("50k","75k","100k","200k","cM")
# Define the variables
variables_vector<-c("Ca","Gyro","Na","pH","Schisto","Zn","Lake_Area",
                              "Shape_PC1","Shape_PC2","Shape_PC3",
                              "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                              "Gill_Raker_L","Gill_Raker_N")
# Analysis to be looped over Radiations
radiations_vector<-c("Alaska","BC","Iceland","Scotland")

# Run for each radiation, for each group, for each variable
mclapply(radiations_vector,function(rad){
  
  # Do the grouping variable
  lapply(grouping_vector,function(group){
    
    # Read in snp_count
    snp_count<-read.table(paste0("data/outlier_beds/",rad,"_",group,"_windows_SNPcount.bed"))
    
    # Do the variable
    lapply(variables_vector,function(var){
      
      # Read in outlier beds
      outliers<-read.table(paste0("outputs/",rad,"_",var,"_cleaned_outliers.bed"))
      
      # Intersect
      outlier_winds<-data.frame(rbindlist(lapply(1:nrow(snp_count),function(x){
        
        # Subset
        tmp<-outliers[outliers$V1 == snp_count$V1[x] &
                        outliers$V2 >= snp_count$V2[x] &
                        outliers$V2 <= snp_count$V3[x],]
        
        if(nrow(tmp) == 0){
          return(cbind(snp_count[x,1:3],data.frame(V4=NA)))
        } else {
          return(cbind(snp_count[x,1:3],data.frame(V4=nrow(tmp))))
        }
        
      })))
      
      # Clean
      outlier_winds<-na.omit(outlier_winds)
      
      # Write
      write.table(outlier_winds,
                  paste0("data/outlier_beds/",group,"/",group,"_windows_",rad,"_",var,"_cleaned_outliers_v2.bed"),
                         row.names=F,quote = F,col.names = F,sep="\t")
      
    })
  })
  
})