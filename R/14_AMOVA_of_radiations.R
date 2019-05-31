#################
# This script performs AMOVA over all freshwater pops nested by radiation
#################

lib<-as.vector(c("vcfR","poppr","data.table","VennDiagram","ggplot2","parallel","dplyr"))
lapply(lib,library,character.only=TRUE)

# Read in the freshwater VCF
vcf_in<-read.vcfR("data/MxF_VCFs/All_fresh_p4r05_FST.vcf")
popmap<-read.table("data/MxF_VCFs/All_fresh_popmap",header=F)
colnames(popmap)<-c("Ind","Pop")

# Get pops into rad order
radiations<-list(as.character(unique(popmap$Pop)[1:19]),
                 as.character(unique(popmap$Pop)[20:37]),
                 as.character(unique(popmap$Pop)[38:55]),
                 as.character(unique(popmap$Pop)[56:73]))
rads<-c("Alaska","BC","Iceland","Scotland")
radmap<-data.frame(rbindlist(lapply(1:4,function(x){return(data.frame(pop=radiations[[x]],
                                                 rad=rads[x]))})))
radmap$rad<-as.character(radmap$rad)
radmap$pop<-as.character(radmap$pop)

# Add to popmap
for(i in 1:nrow(popmap)){
  popmap$Rad[i]<-radmap[as.character(radmap$pop) == as.character(popmap$Pop[i]),"rad"]
}

# Convert to genlight
fresh<-vcfR2genlight(vcf_in)

# Filter popmap
popmap<-popmap[popmap$Ind %in% fresh$ind.names,]

popmap$Continent<-"North America"
popmap[popmap$Rad %in% c("Iceland","Scotland"),"Continent"]<-"Europe"

# Make strata
strata(fresh)<-popmap[,c("Continent","Rad","Pop")]

# Perform amova
amova_res<-poppr.amova(fresh, hier = ~Continent/Rad/Pop)

# Export results to a file
capture.output(amova_res$results, file = "tables/AMOVA_results.txt", append = TRUE)
capture.output(amova_res$componentsofcovariance, file = "tables/AMOVA_covaraince_components.txt", append = TRUE)
