##### Read in the FST files, identify outliers and overlap matrix

# Get packages
lib<-c("data.table","ggplot2","dplyr","parallel","viridis")
lapply(lib,library,character.only=T)

# Read in the outputs
radiations<-c("Alaska","BC","Iceland","Scotland")
dd<-data.frame(rbindlist(mclapply(radiations,function(x){
  
  tmp<-data.frame(fread(paste0("outputs/",x,"_MxF_Fst_50kb_windows.txt")))
  return(tmp)
},mc.cores = 4)))

# Take top 5% of windows and overlap
fresh_pops<-unique(dd$pop)
outlier_list<-lapply(fresh_pops,function(x){
  
  tmp<-dd[dd$pop == x,]
  outliers<-tmp[tmp$FST > quantile(tmp$FST, 0.95),"window_id"]
  return(outliers)
})

# Populate matrix
comparisons_to_make<-combn(1:length(outlier_list),2)
overlap_mat<-matrix(ncol=length(outlier_list),nrow=length(outlier_list))
colnames(overlap_mat)<-fresh_pops
rownames(overlap_mat)<-fresh_pops

# Fill the matrix
for(i in 1:ncol(comparisons_to_make)){
A<-comparisons_to_make[1,i]
B<-comparisons_to_make[2,i]

overlap<-length(Reduce(intersect,list(outlier_list[[A]],outlier_list[[B]])))
if(length(overlap)==0){
  overlap_mat[A,B]<-0
} else {
overlap_prop<-overlap/min(c(length(outlier_list[[A]]),length(outlier_list[[B]])))
overlap_mat[A,B]<-overlap_prop
overlap_mat[B,A]<-overlap_prop
}
}

# Fill the centre line
for(i in 1:length(fresh_pops)){
  overlap_mat[i,i]<-NA
}

# Plot
plot_dd<-melt(overlap_mat)
overlaps_fig<-ggplot(plot_dd,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  scale_fill_viridis()+
  theme(axis.text.x=element_text(angle=90,hjust=1),
        axis.title=element_blank())+
  labs(fill=expression(MxF~F[ST]~Overlap))

# Plot
pdf("figs/MxF_Fst_overlap_proportions.pdf",width=12,height=11)
overlaps_fig
dev.off()

# Export the matrix
write.table(overlap_mat,
            paste0("outputs/MxF_FST_overlap_props.txt"),quote = F,row.names = F,sep="\t")
