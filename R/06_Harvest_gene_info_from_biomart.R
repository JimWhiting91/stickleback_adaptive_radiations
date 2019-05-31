#########################################################################################################
# 06 - Harvester. This script imports parallel associated windows and harvests information on genes contained and their GO if available
#########################################################################################################


lib<-as.vector(c("data.table","patchwork","VennDiagram","ggplot2","parallel","biomaRt"))
lapply(lib,library,character.only=TRUE)

# Define the grouping
grouping_vector<-as.vector(c("50k","75k","100k","200k","genes"))
# Define the variables
variables_vector<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn",
                              "Shape_PC1","Shape_PC2","Shape_PC3",
                              "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                              "Gill_Raker_L","Gill_Raker_N"))
env_pheno_vec<-c(rep("Env",6),
                 rep("Shape",3),
                 rep("Armour",7),
                 rep("Gill",2))

#set up biomart
listMarts()
stickleback=useMart("ensembl",dataset="gaculeatus_gene_ensembl")
listDatasets(stickleback)

# ------
# This script will be repeated over windows function
# ------

window_function<-lapply(1:length(grouping_vector),function(x){
  # Read all of the parallel windows in a list
  parallel_windows<-lapply(1:length(variables_vector),function(y){
    tmp_in<-read.table(paste0("outputs/",grouping_vector[x],"/08/",grouping_vector[x],"_",variables_vector[y],"_Adjacent_Convergence.txt"),header=T,sep = "\t")
    
    # Print variable
    print(grouping_vector[x])
    print(variables_vector[y])
    
    # Get Chr etc
    tmp_in$chr=unlist(tstrsplit(tmp_in$window_id,":",keep=1))
    BPs<-unlist(tstrsplit(tmp_in$window_id,":",keep=2))
    tmp_in$BP1<-as.integer(unlist(tstrsplit(BPs,"-",keep=1)))
    tmp_in$BP2<-as.integer(unlist(tstrsplit(BPs,"-",keep=2)))
  
    
    
    # Now for each row of parallel_tmp, use biomart to find genes and GO
    parallel_biomart<-data.frame(rbindlist(lapply(1:length(tmp_in$window_id),function(z){
      # Get variables
      chr<-tmp_in$chr[z]
      BP1<-tmp_in$BP1[z]
      BP2<-tmp_in$BP2[z]
        
      # Call Biomart
      #Get location information from biomaRt
      tmp_biomart<-getBM(attributes = c("ensembl_gene_id","external_gene_name","go_id","name_1006"),
                        filters= c("chromosome_name","start","end"),
                        values=list(chr,BP1,BP2),
                        mart=stickleback)
      
      # Only return if there are genes
      if(length(tmp_biomart$ensembl_gene_id) > 0){
      
      # Add other columns
      tmp_biomart$Window<-tmp_in$window_id[z]
      tmp_biomart$Radiations<-tmp_in$rads[z]
      tmp_biomart$Variable<-variables_vector[y]
      
      # Return
      return(tmp_biomart)
      } 
      })))
    
    # Set names
    colnames(parallel_biomart)[1:4]<-c("Ensembl ID",
                                       "Common Name",
                                       "GO Term",
                                       "GO Name")
    
    # Replace empty space with NA
    parallel_biomart[parallel_biomart==""]<-NA
    
    # Write to an output
    dir.create(paste0("outputs/",grouping_vector[x],"/06"))
    write.table(parallel_biomart,
                paste0("outputs/",grouping_vector[x],"/06/06_",grouping_vector[x],"_",variables_vector[y],"_parallel_window_genes_biomart_out.txt"),
                row.names = F,sep="\t",quote = F)
    
  })
})
