####### Script for performing marine - freshwater Fst comparisons

# Where are we?
setwd("~/Nottingham/NatEvol_Ecol_Genomic_Parallelism//")

lib<-c("PopGenome","ggplot2","data.table")
lapply(lib,library,character.only=T)

# Set marines
marines<-c("MUD1","LICA","NYPS","OBSM")

# Read in the data
vcfs<-c("Al","Bc","Ic","Sc")
radiations<-c("Alaska","BC","Iceland","Scotland")

# Run in parallel over each radiation
rad_out<-mclapply(1:length(vcfs),function(x){
  
  SNPs<-read.table(paste0("data/",vcfs[x],".SNPs"))[,1:2]
  scafs<-unique(SNPs[,1])
  popmap<-read.table(paste0("data/",vcfs[x],".popmap2"))
  freshwaters<-unique(popmap[!(popmap$V2 %in% marines),"V2"])
  
  # Get comparisons ready
  comparisons<-lapply(1:length(freshwaters),function(pop){
    pop1<-popmap[popmap$V2 == freshwaters[pop],"V1"]
    pop2<-popmap[popmap$V2 == marines[x],"V1"]
    return(list(pop1,pop2))
  })
  
  # Loop over scafs within each pop
  scafs_out<-data.frame(rbindlist(lapply(1:length(scafs),function(i){
  print(scafs[i])  
    # Get chr
    tid<-as.character(scafs[i])
    
    # Get max SNP
    final<-max(SNPs[SNPs$V1 == tid,"V2"])
    
    # Read in the SNPs
    vcf_in<-readVCF(paste0("data/",vcfs[x],".vcf.gz"),tid=tid,frompos=1,topos=final,numcols = 1000,include.unknown = T)
    
    # Perform Fst calculation for each comparison
    comps_out<-data.frame(rbindlist(lapply(1:length(comparisons),function(comp){
      
      # Set the comparison
      vcf2<-set.populations(vcf_in,do.call("list",comparisons[[comp]]),diploid=TRUE)
      vcf2<-set.outgroup(vcf2,FALSE)
      
      # Do window transform to 50kb
      #split into windows
      
      # Check first whether region is large enough
      split_region<-as.integer(strsplit(vcf2@region.names," - ")[[1]][2])
      if(split_region > 50000){
      win_SNP<-sliding.window.transform(vcf2,width=50000,jump=50000,type=2)
      } else {
      win_SNP<-vcf2
      }
      #do pop stats
      win_SNP <-F_ST.stats(win_SNP,mode = "nucleotide")
      
      win_SNP <-neutrality.stats(win_SNP,FAST=FALSE)

      #Get centre of the window
      genome.pos <- sapply(win_SNP@region.names, function(x){
        split <- as.integer(strsplit(x," ")[[1]][c(1,3)])
        split[1]<-split[1]-2
        split[2]<-split[2]-1
        
        # Catch bad windows
        if(split[1] < 50000){split[1] <- 0}
        if(split[2] < 50000){split[2] <- 50000}
        
        val   <- paste0(tid,":",as.integer(split[1]),"-",as.integer(split[2]))
        return(val)
      })
      
      #output results matrix
      FST_mat<-data.frame(chr = tid,
                     FST=win_SNP@nucleotide.F_ST,
                     win_SNP@n.segregating.sites,
                     win_SNP@nuc.diversity.within/50000,
                     genome.pos)
      colnames(FST_mat)<-c("chr","FST","fresh_sites","marine_sites","fresh_pi","marine_pi","window_id")
      FST_mat$pop<-freshwaters[comp]
      FST_mat$radiation<-radiations[x]
      
      # Correct low Fst
      FST_mat[FST_mat$FST < 0, "FST"]<-0
      
      # Return
      return(na.omit(FST_mat))
    })))
   
  return(comps_out)   
  })))

write.table(scafs_out,
            paste0("outputs/",radiations[x],"_MxF_Fst_50kb_windows.txt"),quote = F,row.names = F,sep="\t")
},mc.cores=4)
