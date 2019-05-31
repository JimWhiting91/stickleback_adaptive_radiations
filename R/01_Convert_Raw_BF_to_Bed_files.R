####
# This script takes the sorted, averaged, otherwise RAW Bayenv output and converts it into bed_files that can be used downstream
# We need an output bed file for each variable, for each radiation
####

#pH [2:4]
#Zn [5:7]
#Na [8:10]
#Ca [11:13]
#Gyro [14:16]
#Schisto [17:19]
#SHAPE.IND.PC1 [20:22]
#SHAPE.IND.PC2 [23:25]
#SHAPE.IND.PC3 [26:28]
#ARMOUR.IND.NOPLATE.PC1 [29:31]
#ARMOUR.IND.NOPLATE.PC2 [32:34]
#ARMOUR.IND.WITHPLATE.PC1 [35:37]
#ARMOUR.IND.WITHPLATE.PC2 [38:40]
#ARMOUR.IND.WITHPLATE.PC3 [41:43]
#ARMOUR.IND.WITHPLATE.PC4 [44:46]
#Avg.resDS1 [47:49]
#Avg.resDS2 [50:52]
#Avg.resPS [53:55]
#Avg.resLP [56:58]
#Avg.resHP [59:61]
#Avg.resBAP [62:64]
#Avg.n.plate [65:67]
#St.Length [68:70]
#Gill.Raker.N [71:73]
#Resid.Raker.L [74:76]

# First we create a vector of variable names
variable_names<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn",
                            "Shape_PC1","Shape_PC2","Shape_PC3",
                            "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                            "Gill_Raker_L","Gill_Raker_N"))
# Now we list the corresponding column identifiers
columns_ID_vector<-as.vector(c(11,14,8,2,17,5,
                               20,23,26,
                               47,50,53,56,59,62,65,
                               71,74))

# Also a vector of Radiation Labels
radiation_names<-as.vector(c("Alaska","BC","Iceland","Scotland"))

# The same analysis will be performed over each radiation, so write a function
radiation_cleaner<-function(x){
  # Read in the data
  dd<-read.table(paste0("data/raw_bayenv_outputs/",radiation_names[x],"_EnvPheno_sorted_AVG_v2.txt"),header=F)
  
  # Read in the SNP IDs
  dd_SNP<-read.table(paste0("data/raw_bayenv_outputs/",radiation_names[x],"_CHRBP"),header=F)
  
  # Remove second column because its all 0s
  dd<-dd[,-2]
  
  # For each variable, we create a bed file with chrom, bp start, bp end, Index, BF, Spearman, Pearson
  for (i in 1:length(variable_names)){
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
  
  # May as well also make the outlier bed file whilst we're here
    # Outliers defined as logBF > 1.5 and top 5% of ABSOLUTE (direction of effect doesn't matter) spearman quantile
    dd_var_outliers<-dd_var[log10(dd_var$BF) > 1.5 & abs(dd_var$Spear) > quantile(abs(dd_var$Spear),probs = 0.95),]
  
  # Write both files to output
    write.table(dd_var,
                paste0("outputs/",radiation_names[x],"_",variable_names[i],"_cleaned.bed"),
                col.names = F,
                row.names = F,
                quote = F)
    write.table(dd_var_outliers,
                paste0("outputs/",radiation_names[x],"_",variable_names[i],"_cleaned_outliers.bed"),
                col.names = F,
                row.names = F,
                quote = F)
  }
}
lapply(1:length(radiation_names),radiation_cleaner)
