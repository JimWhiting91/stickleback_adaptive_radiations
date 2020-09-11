# Process Lake Data
library(popbio)

rads<-c("Alaska","BC","Iceland","Nuist")

# Read and sort
data<-lapply(rads,function(x){
  tmp<-lapply(1:10,function(y){
    run<-read.table(paste0("data/raw_bayenv_outputs/lake_area/",x,"_bayenv_lake_tests_",y,".bf"))
    run<-run[order(run$V1),]
    return(as.matrix(run[,2:4]))
  })
  
  # Average the list
  tmp<-popbio::mean.list(tmp)
  
  # Write the output
  write.table(tmp,
              paste0("data/raw_bayenv_outputs/lake_area/",x,"_bayenv_lake_tests_AVG_sorted.txt"), 
              quote=F,row.names=F,sep = "\t",col.names=F)
})
