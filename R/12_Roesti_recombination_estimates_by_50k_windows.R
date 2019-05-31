###################################################################################################
# This script asks how recombination rate varies between windows, outliers, and convergent outliers 
###################################################################################################

# Read in res
lib<-as.vector(c("dplyr","data.table","ggplot2","parallel","ggpubr","ggsignif"))
lapply(lib,library,character.only=TRUE)

#
# Run over 50kb results
#

# Set up vectors
rads<-c("Alaska","BC","Iceland","Scotland")
variables_vector<-c("Ca","Gyro","Na","pH","Schisto","Zn",
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

# Read in recomb rates
recomb<-read.table("data/roesti_BROADS1_recomb_rates.txt",header=T)
recomb$chr<-paste0("group",as.roman(recomb$lg))

# Assign rec rate as weighted mean
windows$recomb<-NA
for(i in 1:nrow(windows)){
  recomb_sub<-recomb[recomb$chr == windows$chr[i],]
  recomb_tmp<-recomb_sub[recomb_sub$pos2 >= windows$BP1[i] &
                           recomb_sub$pos1 <= windows$BP2[i],]
  if(length(recomb_tmp[,1]) > 0){
    if(recomb_tmp$pos1[1] < windows$BP1[i]){
      recomb_tmp$pos1[1]<-windows$BP1[i]
    }
    if(recomb_tmp$pos2[nrow(recomb_tmp)] > windows$BP2[i]){
      recomb_tmp$pos2[nrow(recomb_tmp)]<-windows$BP2[i]
    }
    diffs<-recomb_tmp$pos2-recomb_tmp$pos1
    props<-sapply(diffs,function(x){return(x/sum(diffs))})
    windows$recomb[i]<-weighted.mean(recomb_tmp$rate,props)
  }
}

# Get rid of windows with no data
recomb_windows<-na.omit(windows)
recomb_windows$recomb<-as.numeric(recomb_windows$recomb)
recomb_windows$outliers<-as.factor(recomb_windows$outliers)

# Summarise
recomb_windows %>%
  group_by(outliers) %>%
  summarise(mean=mean(recomb),
            median=median(recomb),
            sd=sd(recomb),
            N=n())

# Plot
rec_fig<-ggplot(recomb_windows,aes(x=outliers,y=log10(recomb)))+
  geom_violin(draw_quantiles = 0.5,fill="gray")+
  geom_signif(comparisons = list(c("Neutral", "Outlier"),
                                 c("Neutral","Overlap"),
                                 c("Outlier","Overlap")),   
              map_signif_level=TRUE,
              step_increase = 0.1)+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  labs(y=expression(Recombination~(Log[10]~cM/Mb)),
       x="")

# Plot to pdf
pdf("figs/Recombination_estimates_of_outliers_overlap.pdf")
rec_fig
dev.off()

# Is there a different in recomb rate between treatment groups?
kruskal.test(recomb ~ outliers, data = recomb_windows)

# between which groups?
pairwise.wilcox.test(recomb_windows$recomb, recomb_windows$outliers,
                     p.adjust.method = "BH")
