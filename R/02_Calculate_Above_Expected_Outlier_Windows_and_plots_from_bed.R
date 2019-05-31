#################
# This script takes the bed files that group outlier SNPs into windows/genes and performs binomial expectation outlier identification
# to find outlier windows/genes

# Outliers were defined as log(BF) > 1.5 and above 95% quantile of spearmans

#################

lib<-as.vector(c("dplyr","data.table","VennDiagram","ggplot2","parallel","ggpubr"))
lapply(lib,library,character.only=TRUE)

# ------
# This script will be repeated over windows and genes so we put in function
# ------

# Define the grouping
grouping_vector<-as.vector(c("50k","75k","100k","200k","cM"))
# Define the variables
variables_vector<-as.vector(c("Ca","Gyro","Na","pH","Schisto","Zn",
                              "Shape_PC1","Shape_PC2","Shape_PC3",
                              "DS1","DS2","PS","LP","HP","BAP", "Plate_N",
                              "Gill_Raker_L","Gill_Raker_N"))
# Analysis to be looped over Radiations
radiations_vector<-as.vector(c("Alaska","BC","Iceland","Scotland"))

# Define the function
binomial_exp_calculations<-function(x){
  
  var_ass_count<-list()
  
  #Loop over variables
  for (j in 1:length(variables_vector)){
    
    # Need radiation slots in lists for saving figures and outputs
    manhattan_list<-list()
    top_outliers_list<-list()
    scatter_list<-list()
    scatter_list2<-list()
    
    rad_ass_count<-list()
    
    for (i in 1:length(radiations_vector)){
      
      print(grouping_vector[x])
      print(radiations_vector[i])
      print(variables_vector[j])
      
      #Read in SNP count file
      dd_SNP_count<-read.table(paste0("data/outlier_beds/",radiations_vector[i],"_",grouping_vector[x],"_windows_SNPcount.bed"),header=F,sep='\t',fill = TRUE)
      dd_SNP_count$window_id<-paste0(dd_SNP_count$V1,":",dd_SNP_count$V2,"-",dd_SNP_count$V3)
      
      #Read in data for outlier SNPs
      dd_outliers<-read.table(paste0("data/outlier_beds/",grouping_vector[x],"/",grouping_vector[x],"_windows_",radiations_vector[i],"_",variables_vector[j],"_cleaned_outliers.bed_sorted.bed"),
                              header=F,sep='\t',fill = T)

      dd_outliers$window_id<-paste0(dd_outliers$V1,":",dd_outliers$V2,"-",dd_outliers$V3)
      
      #Combine SNP counts for windows back on to cutoffs
      SNP_counts_outliers<-dd_SNP_count[dd_SNP_count$window_id %in% dd_outliers$window_id , ]
      dd_outliers$SNP_N<-SNP_counts_outliers$V4
      
      dd_outliers$SNP_N<-as.integer(dd_outliers$SNP_N)
      
      #Tidy
      colnames(dd_outliers)[1:4]<-c("chr","BP1","BP2","outlier_N")
      
      #Calculate cutoff
      SNP_sim<-as.data.frame(seq(1,max(as.integer(dd_outliers$SNP_N)),by=1))
      colnames(SNP_sim)<-"SNP_N"
      
      # We need to know what the probability for expected number of SNPs is which = N of outlier SNPs over total SNPs
      p<-sum(dd_outliers$outlier_N)/sum(dd_SNP_count$V4)
      
      # Here we calculate binomial expectation for 0.99 quantile
      for (k in 1:max(dd_outliers$SNP_N)){
        SNP_sim$exp[k]<-qbinom (0.99, k, p)
      }
      
      #Calculate outliers
      dd_outliers$exp<-qbinom (0.99, dd_outliers$SNP_N, p)
      dd_outliers$index<-seq(1:length(dd_outliers$chr))
      
      # Write intermediate output to file
      top_outliers<-subset(dd_outliers,outlier_N > exp)
      top_outliers_index<-top_outliers$index
      
      # Save information about associated SNP count and associated window count
      snp.window.count.dd<-data.frame(ass.snp.count=sum(dd_outliers$outlier_N),
                                      ass.window.count=length(top_outliers_index),
                                      window_size=grouping_vector[x],
                                      radiation=radiations_vector[i],
                                      variable=variables_vector[j])
      rad_ass_count[[i]]<-snp.window.count.dd
      
      # We also want to save the top_outliers into single variable files
      top_outliers$variable<-variables_vector[j]
      # Add radiation label
      top_outliers$radiation<-radiations_vector[i]
      # Save to list
      top_outliers_list[[i]]<-as.data.table(top_outliers)
      #Highlight in dataset
      dd_outliers$top<-rep("No",length(dd_outliers$chr))
      dd_outliers$top[top_outliers_index]<-c("Yes")
      
      
      #Plot basic first
      g_outlier<-ggplot(dd_outliers,aes(x=SNP_N,y=outlier_N))+
        geom_point(aes(color=top))+
        geom_line(data=SNP_sim,aes(y=exp),color=c("blue1"),alpha=0.5)+
        xlab(paste0("SNP count for ",grouping_vector[x]))+
        ylab(paste0("Outlier count for ",grouping_vector[x]))+
        scale_color_manual(values=c("black","red2"))+
        ggtitle(paste0(radiations_vector[i]," ",variables_vector[j]))+
        theme_bw()+
        theme(axis.title.x = element_text(face="bold", size=14),
              axis.title.y = element_text(face="bold", size=14),
              legend.position = "none",
              title = element_text(face="bold",size=16))
      
      # #Make a label data frame to add on later (because for some reason it won't do it properly)
      # label_dd<-matrix(nrow = 1,ncol = 3)
      # label_dd[1,1]<-paste0(radiations_vector[i]," ",variables_vector[j])
      # label_dd[1,3]<-max(dd_outliers$outlier_N-0.5)
      # label_dd[1,2]<-max(dd_outliers$SNP_N/4)
      # colnames(label_dd)<-c("label","x","y")
      # label_dd<-as.data.frame(label_dd)
      # label_dd$x<-as.numeric(as.character(label_dd$x))
      # label_dd$y<-as.numeric(as.character(label_dd$y))
      # 
      scatter_list[[i]]<-g_outlier
      # scatter_list2[[i]]<-label_dd
      
      # Now we need to plot the manhattans, so combine dd_SNP_count with dd_outliers and match up columns
      dd_SNP_count$SNP_N<-dd_SNP_count$V4
      dd_SNP_count$V4<-rep(0,length(dd_SNP_count$V4))
      dd_SNP_count$exp<-qbinom(0.99,dd_SNP_count$SNP_N,p)
      dd_SNP_count$index<-seq(1:length(dd_SNP_count$V1))
      dd_SNP_count$top<-rep("No",length(dd_SNP_count$V1))
      # Match colnames
      colnames(dd_SNP_count)<-colnames(dd_outliers)
      
      # Now we remove from dd_SNP_count the outlier rows
      dd_SNP_count_no_outliers<-dd_SNP_count[!(dd_SNP_count$window_id %in% dd_outliers$window_id),]
      
      # Now combine the removed total file with the outlier rows so that outlier_N column is appropriate
      dd_manhat<-rbind(dd_SNP_count_no_outliers,dd_outliers)
      
      # Re-sort the data 
      dd_manhat<-dd_manhat[order(dd_manhat$chr,dd_manhat$BP1),]
      
      # Calculate residual against expected as outlier_N - Exp
      dd_manhat$resid<-dd_manhat$outlier_N-dd_manhat$exp
      
      # Provide the x-axis index value
      dd_manhat$window_index<-seq(1:length(dd_manhat$chr))
      
      # Need to define odd and even chromosomes for colour
      odd_chr<-c("groupI","groupIII","groupV","groupVII","groupIX","groupXI","groupXIII","groupXV","groupXVII","groupXXI")
      even_chr<-c("groupII","groupIV","groupVI","groupVIII","groupX","groupXII","groupXIV","groupXVI","groupXVIII","groupXX")
      dd_manhat$chr_colour<-rep("None",length(dd_manhat$chr))
      dd_manhat[dd_manhat$chr %in% odd_chr,"chr_colour"]<-"ODD"
      dd_manhat[dd_manhat$chr %in% even_chr,"chr_colour"]<-"EVEN"
      
      # Assign for outside combination with other radiations
      dd_manhat$radiation<-rep(radiations_vector[i],length(dd_manhat$chr))
      manhattan_list[[i]]<-data.table(dd_manhat)
    }
    #rbindlist
    top_outliers<-as.data.frame(rbindlist(top_outliers_list))
    
    # Save ass counts
    rad_ass_count_out<-data.frame(rbindlist(rad_ass_count))
    var_ass_count[[j]]<-rad_ass_count_out
    
    # Save to output
    write.table(top_outliers,
                paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_",variables_vector[j],"_top_outliers_by_rad.txt"),
                row.names = F,
                quote = F)
    
    # Also need to calculate.overlap values to look for convergence
    dd_overlap<-calculate.overlap(list(top_outliers_list[[1]]$window_id,
                                       top_outliers_list[[2]]$window_id,
                                       top_outliers_list[[3]]$window_id,
                                       top_outliers_list[[4]]$window_id))
    
    # Turn calculate.overlap into human readable output with column headers
    dd_overlap_counts<-list()
    for (i in seq(1,15,by=1)){
      dd_overlap_counts[[i]]<-length((as.vector(dd_overlap[[i]])))
    }
    dd_overlap_counts<-data.frame(dd_overlap_counts)
    colnames(dd_overlap_counts)<-c("All_4",
                                   "Alaska_&_BC_&_Iceland","Alaska_&_BC_&_Scotland","Alaska_&_Iceland_&_Scotland","BC_&_Iceland_&_Scotland",
                                   "Alaska_&_BC","Alaska_&_Iceland","Alaska_&_Scotland","BC_&_Iceland","BC_&_Scotland","Iceland_&_Scotland",
                                   "Alaska","BC","Iceland","Scotland")
    
    # Write to an output
    write.table(dd_overlap_counts,
                paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_",variables_vector[j],"Venn_overlaps.txt"),
                row.names = F,
                quote = F,
                sep="\t")
    
    # Stack scatters
    p1<-scatter_list[[1]]
    p2<-scatter_list[[2]]
    p3<-scatter_list[[3]]
    p4<-scatter_list[[4]]
    
    pdf(paste0("figs/02/",grouping_vector[x],"/02_",grouping_vector[x],"_",variables_vector[j],"_binom_exp_outlier_scatters.pdf"))
    print(ggarrange(p1,p2,p3,p4,ncol=2,nrow=2))
    dev.off()
    
    
    # Back to manhattans, combine each of the saved outputs and make figure using facet_wraps
    dd_manhat_all<-data.frame(rbindlist(manhattan_list))
    
    # Adjust genomic location based on size of chr
    chr<-as.vector(c("groupI",
                     "groupII",
                     "groupIII",
                     "groupIV",
                     "groupV",
                     "groupVI",
                     "groupVII",
                     "groupVIII",
                     "groupIX",
                     "groupX",
                     "groupXI",
                     "groupXII",
                     "groupXIII",
                     "groupXIV",
                     "groupXV",
                     "groupXVI",
                     "groupXVII",
                     "groupXVIII",
                     "groupXIX",
                     "groupXX",
                     "groupXXI"))
    location<-as.numeric(c(0,
                           28185914,
                           51481566,
                           68280072,
                           100913020,
                           113164417,
                           130248092,
                           158185535,
                           177554239,
                           197803718,
                           213461158,
                           230167210,
                           248568277,
                           268651407,
                           283897868,
                           300096632,
                           318212420,
                           332815561,
                           349098277,
                           369338937,
                           389071008))
    chr_locations<-data.frame(cbind(chr,location))
    chr_locations$location<-as.numeric(as.character(chr_locations$location))
    
    dd_manhat_all$wind_start<-seq(1:length(dd_manhat_all$chr))
    dd_manhat_all$wind_end<-seq(1:length(dd_manhat_all$chr))
    
    #First remove scaffolds
    dd_manhat_all<-dd_manhat_all[dd_manhat_all$chr %in% chr,]
    
    for (i in 1:length(dd_manhat_all$chr)){
      chr_info<-(chr_locations[chr_locations$chr %in% dd_manhat_all$chr[i],2])
      dd_manhat_all$wind_start[i]<-dd_manhat_all$BP1[i]+chr_info
      dd_manhat_all$wind_end[i]<-dd_manhat_all$BP2[i]+chr_info
    }
    
    # Extract duplicated convergent windows into new dataframe to plot over
    dd_manhat_top<-dd_manhat_all[dd_manhat_all$top == "Yes",]
    dd_manhat_top<-dd_manhat_top[duplicated(dd_manhat_top$window_id),]
    dd_manhat_top<-dd_manhat_all[dd_manhat_all$window_id %in% dd_manhat_top$window_id,]
    dd_manhat_top<-dd_manhat_top[dd_manhat_top$top == "Yes",]
    
    # Count how many times each appears to determine level of convergence
    tryCatch(for (i in 1:length(dd_manhat_top$chr)){
      dd_manhat_top$convergence_N[i]<-as.character(sum(dd_manhat_top$window_id == dd_manhat_top$window_id[i]))
    }, error=function(e) NULL)
    
    # Set plotting levels
    dd_manhat_top$plot_y[dd_manhat_top$radiation=="Alaska"]<-as.numeric(max(dd_manhat_all$resid[dd_manhat_all$radiation=="Alaska"]))
    dd_manhat_top$plot_y[dd_manhat_top$radiation=="BC"]<-as.numeric(max(dd_manhat_all$resid[dd_manhat_all$radiation=="BC"]))
    dd_manhat_top$plot_y[dd_manhat_top$radiation=="Iceland"]<-as.numeric(max(dd_manhat_all$resid[dd_manhat_all$radiation=="Iceland"]))
    dd_manhat_top$plot_y[dd_manhat_top$radiation=="Scotland"]<-as.numeric(max(dd_manhat_all$resid[dd_manhat_all$radiation=="Scotland"]))
    
    #Add x axis to plot
    chrlabs<-c(14092957,
               39833740,
               59880819,
               84596546,
               107038718.5,
               121706254.5,
               144216813.5,
               167869887,
               187678978.5,
               205632438,
               221814184,
               239367743.5,
               258609842,
               276274637.5,
               291997250,
               309154526,
               325513990.5,
               340956919,
               379204972.5,
               399000000)
    chrNam<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","20","21")
    
    # Plot manhattans in ggplot
    manhattan_p<-ggplot(dd_manhat_all,aes(x=wind_start,y=resid,colour=chr_colour))+
      geom_segment(aes(x = wind_start, y = 0, xend = wind_end, yend = resid))+
      facet_wrap(~radiation,ncol=1,scales= "free_y",strip.position="right")+
      theme_bw()+
      theme(axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16,angle=90),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            strip.text = element_text(face="bold", size =14),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position = "none",
            title=element_text(size=16,face="bold"))+
      scale_colour_manual(breaks = c("ODD","EVEN"),
                          values = c("black","grey45"))+
      ggtitle(variables_vector[j])+
      ylab("Residual associated SNPs over expected")+
      xlab("Genome Position")+
      ylim(-1,NA)+
      scale_x_continuous(name="Chromosome",breaks=chrlabs,labels=chrNam) +
      geom_point(data=dd_manhat_top[dd_manhat_top$convergence_N == "2",],inherit.aes = F,aes(x=wind_start,y=plot_y+0.5),colour="red2",shape=16,size=2)+
      geom_point(data=dd_manhat_top[dd_manhat_top$convergence_N == "3",],inherit.aes = F,aes(x=wind_start,y=plot_y+1.5),colour="forestgreen",shape=17,size=3)+
      geom_point(data=dd_manhat_top[dd_manhat_top$convergence_N == "4",],inherit.aes = F,aes(x=wind_start,y=plot_y+2.5),colour="gold2",shape=18,size=4)
    
    pdf(paste0("figs/02/",grouping_vector[x],"/02_",grouping_vector[x],"_",variables_vector[j],"_manhattans.pdf"))
    print(manhattan_p)
    dev.off()
    
  }
  
  var_ass_count_out<-data.frame(rbindlist(var_ass_count))
  return(var_ass_count_out)
}

#Run through trycatch because errors out if there are no convergent windows to plot in manhattans
total_ass_count<-lapply(1:length(grouping_vector),binomial_exp_calculations)
total_ass_count_out<-data.frame(rbindlist(total_ass_count))
# Save total ass count file
write.table(total_ass_count_out,
            "tables/Counts_Associated_SNP_Window_All_vars_All_rads.txt",
            sep="\t",row.names = F,quote = F)


##### Histograms of count data results #####
# For each window size, produce side by side comparisons of associated SNP counts and associated window counts
hist_graphs<-mclapply(1:length(grouping_vector),function(x){
  # Subset
  tmp_counts<-total_ass_count_out[total_ass_count_out$window_size == grouping_vector[x],]
  
  # Make facet histogram of SNP count
  g1<-ggplot(tmp_counts,aes(ass.snp.count))+
    geom_histogram()+
    facet_wrap(~radiation,ncol = 1)+
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.title.y = element_text(face="bold", size=16,angle=90),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          strip.text = element_text(face="bold", size =14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "none")+
    xlab("Associated SNP Count")
  
  # Make facet histogram of Window count
  g2<-ggplot(tmp_counts,aes(ass.window.count))+
    geom_histogram()+
    facet_wrap(~radiation,ncol = 1)+
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.title.y = element_text(face="bold", size=16,angle=90),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          strip.text = element_text(face="bold", size =14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "none")+
    xlab(paste0("Associated Window Count (",grouping_vector[x],")"))
  
  # Plot side by side
  hists<-ggarrange(g1,g2,ncol = 2,labels = "AUTO")
  return(hists)
})

# Write all to pdf
pdf("figs/02/Histograms_SNP_Window_Counts_all.pdf",height=12,width = 9)
for (i in 1:length(hist_graphs)){
  print(hist_graphs[[i]])
}
dev.off()

##### Are there differences in the numbers of called SNPs or windows across variables? #####
# First do SNP test
SNP_tmp<-total_ass_count_out[total_ass_count_out$window_size == grouping_vector[x],]
# Model
snp.glm<-glm(ass.snp.count~variable+radiation,SNP_tmp,family="gaussian")
# Check residuals
print(shapiro.test(snp.glm$residuals)) # OK
print(anova(snp.glm,test="F"))

for (i in 1:length(grouping_vector)){
  tmp_counts<-total_ass_count_out[total_ass_count_out$window_size == grouping_vector[i],]
  # Model
 window.glm<-glm(ass.window.count~variable+radiation,tmp_counts,family="gaussian")
  # Check residuals
  print(shapiro.test(window.glm$residuals)) # OK
  print(anova(window.glm,test="F"))
}


##### Stitch together outputs into human readable format #####

for (x in 1:length(grouping_vector)){
  data_list<-list()
  for (j in 1:length(variables_vector)){
    dd<-read.table(paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_",variables_vector[j],"Venn_overlaps.txt"),header = T,sep="\t")
    data_list[[j]]<-as.data.table(dd)
  }
  all_Venn<-as.data.frame(rbindlist(data_list))
  rownames(all_Venn)<-variables_vector
  
  summed_outliers<-as.vector(colSums(all_Venn[,1:15]))
  all_Venn[length(variables_vector)+1,]<-summed_outliers
  rownames(all_Venn)[length(rownames(all_Venn))]<-"TOTAL"
  colnames(all_Venn)<-c("All_4",
                        "Alaska_&_BC_&_Iceland","Alaska_&_BC_&_Scotland","Alaska_&_Iceland_&_Scotland","BC_&_Iceland_&_Scotland",
                        "Alaska_&_BC","Alaska_&_Iceland","Alaska_&_Scotland","BC_&_Iceland","BC_&_Scotland","Iceland_&_Scotland",
                        "Alaska","BC","Iceland","Scotland")
  write.table(all_Venn,
              paste0("outputs/",grouping_vector[x],"/",grouping_vector[x],"_Venn_overlaps.txt"),
              quote = F)
}

##### How do N SNPs and Residual scores of overlaps compare with general outliers #####
# For each window size, examine outliers
hist_graphs<-mclapply(1:length(grouping_vector),function(x){
  
  # We will pool all the outliers together across variables for a single analysis (still grouped by window size)
  N_Outputs<-list.files(paste0("outputs/",grouping_vector[x]))
  
  # Get the outlier files
  to_keep<-N_Outputs[grep("top_outliers_by_rad",N_Outputs)]
  
  # Read in
  out_dd<-data.frame(rbindlist(mclapply(1:length(to_keep),function(y){
    tmp_in<-fread(paste0("outputs/",grouping_vector[x],"/",to_keep[y]),fill=T)
    return(tmp_in)
  }),fill = T))
  
  # Residual calc
  out_dd$resid<-out_dd$outlier_N-out_dd$exp
  
  # Asign 'Overlap' or 'Unique'
  out_dd$overlap<-"Unique"
out_dd<-data.frame(rbindlist(lapply(1:length(variables_vector),function(y){
  
  # Subset by variable
  tmp<-out_dd[out_dd$variable == variables_vector[y],]
  
  # Work backwards
  dups<-tmp$window_id[duplicated(tmp$window_id)]
  tmp[tmp$window_id %in% dups,"overlap"]<-"Overlap"
  return(tmp)
})))  

# Plot histograms for N SNPs and Resid
SNP_N_hist<-ggplot(data=out_dd,aes(x=SNP_N,fill=overlap))+
  geom_histogram(bins=20,col="black")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text = element_text(size=20))+
  facet_wrap(~overlap,scales="free_y")+
  labs(y="Count",x=expression(SNP~N~window^-1))

resid_hist<-ggplot(data=out_dd,aes(x=resid,fill=overlap))+
  geom_histogram(bins=20,col="black")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text = element_text(size=20))+
  facet_wrap(~overlap,scales="free_y")+
  labs(y="Count",x="Residual outlier SNP above expected")

# Analyse with Poisson GLMs
# SNP N
out_dd %>% group_by(overlap) %>% summarise(mean=mean(SNP_N),sd=sd(SNP_N))
glm1<-glm(SNP_N~overlap,data=out_dd,family="poisson")
drop1(glm1,test = "Chisq")

# Residual
out_dd %>% group_by(overlap) %>% summarise(mean=mean(resid),sd=sd(resid))
glm2<-glm(resid~overlap,data=out_dd,family="poisson")
drop1(glm2,test = "Chisq")

return(list(SNP_N_hist,resid_hist))
})

# Plot histograms away
pdf("figs/Histograms_overlap_VS_unique_windows.pdf")
for(i in 1:length(hist_graphs)){
  print(ggarrange(hist_graphs[[i]][[1]],
                  hist_graphs[[i]][[2]],ncol=1,nrow=2,labels = "AUTO"))
}
dev.off()
