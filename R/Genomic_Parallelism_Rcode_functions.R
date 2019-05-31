## Function to examine overlap of windows 

condense_windows<-function(windows){
  if(length(windows)==0){
    return(NA)
  } else {
    # Make data.frame
    tmp<-data.frame(outlier_start=NA,
                    outlier_end=NA,
                    chr=unlist(tstrsplit(windows,":",keep=1)))
    
    BPs<-unlist(tstrsplit(windows,":",keep=2))
    tmp$BP1<-as.integer(unlist(tstrsplit(BPs,"-",keep=1)))
    tmp$BP2<-as.integer(unlist(tstrsplit(BPs,"-",keep=2)))
    # Find starts of windows
    tmp$outlier_start[1]<-tmp$BP1[1]
    tmp$outlier_end[nrow(tmp)]<-tmp$BP2[nrow(tmp)]
    # Loop over to fill
    for(i in 2:(nrow(tmp))){
      if(tmp$BP1[i] == (tmp$BP2[i-1])) {
      } else {
        tmp$outlier_start[i]<-tmp$BP1[i]
      }
    }
    for(i in 1:(nrow(tmp)-1)){
      if(tmp$BP2[i] == (tmp$BP1[i+1])) {
      } else {
        tmp$outlier_end[i]<-tmp$BP2[i]
      }
    }
    
    # Now loop over to get rid of NAs
    for(i in 1:nrow(tmp)){
      if(is.na(tmp$outlier_start[i])==TRUE){
        tmp$outlier_start[i] <- tmp$outlier_start[i-1] 
      }
    }
    for(i in nrow(tmp):1){
    if(is.na(tmp$outlier_end[i])==TRUE){
      tmp$outlier_end[i] <- tmp$outlier_end[i+1] 
    }
    }
    
    # Remove the NAs
    out<-paste0(tmp$chr,":",tmp$outlier_start,"-",tmp$outlier_end)
    return(out)
  }
}

# Input of function needs to be list of vectors
overlap_detect<-function(outliers=NULL,extra=NULL,window_size=NULL){
  
  # Set up parallel
  library(parallel)
  library(doParallel)
  library(tidyr)
  library(VennDiagram)
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  outlier_list<-lapply(1:length(outliers),function(x){
    
    # Expand out each vector
    tidy_tmp<-tidyr::separate(as.data.frame(outliers[[x]]),col="outliers[[x]]",c("chr","BP"),":") %>%
      tidyr::separate(col=BP,into=c("BP1","BP2"),"-")
    
    # Add/Subtract Extra for each
    tidy_tmp$BP1<-as.integer(tidy_tmp$BP1)-as.numeric(extra)
    tidy_tmp$BP2<-as.integer(tidy_tmp$BP2)+as.numeric(extra)
    
    # Return to list
    return(tidy_tmp)
  })
  
  # Perform all possible levels of comparison
  output_list<-list()
  for (i in length(outliers):2){
    compare_sets<-as.data.frame(combn(seq(1:length(outliers)),m=i))
    
    output_list2<-list()
    # Split compare_pairs into individual comparison columns
    for (j in 1:length(compare_sets)){
      comp_tmp<-outlier_list[[(compare_sets[1,j])]]
      
      # Assess against each of the other data frames
      for (k in 2:length(compare_sets$V1)){
        comp_tmp2<-outlier_list[[(compare_sets[k,j])]]
        
        # Check for convergence along each row
        overlap_vec<-unlist(foreach(l=1:length(comp_tmp$chr)) %dopar% {
          length(rownames(unique(rbind(comp_tmp2[comp_tmp2$chr == comp_tmp[l,1] & comp_tmp2$BP1 >= comp_tmp[l,2] & comp_tmp2$BP1 <= comp_tmp[l,3],],
                                       comp_tmp2[comp_tmp2$chr == comp_tmp[l,1] & comp_tmp2$BP2 >= comp_tmp[l,2] & comp_tmp2$BP2 <= comp_tmp[l,3],]))))
        })
        # Add overlap_vec back onto comp_tmp
        comp_tmp<-cbind(comp_tmp,overlap_vec)
      }
      # Remove non-convergence
      comp_tmp_out<-comp_tmp[(apply(comp_tmp, 1, function(row) all(row !=0 ))),]
      # Save output
      output_list2[[j]]<-paste0(comp_tmp_out$chr,":",comp_tmp_out$BP1,"-",comp_tmp_out$BP2)
    }
    output_list[[i-1]]<-output_list2
  }
  # Remove duplicates from the output list and return
  un <- unlist(unlist(rev(output_list),recursive = F))
  res <- Map(`[`, unlist(rev(output_list),recursive = F), relist(!duplicated(un), skeleton = unlist(rev(output_list),recursive = F)))
  
  # Tidy output
  res2<-lapply(1:length(res),function(x){
    
    # Expand out each vector
    tidy_tmp<-tidyr::separate(as.data.frame(res[[x]]),col="res[[x]]",c("chr","BP"),":") %>%
      tidyr::separate(col=BP,into=c("BP1","BP2"),"-")
    
    # Add/Subtract Extra for each
    tidy_tmp$BP1<-as.integer(tidy_tmp$BP1)+as.integer(extra)
    tidy_tmp$BP2<-as.integer(tidy_tmp$BP2)-as.integer(extra)
    
    # Remove NAs
    tidy_tmp<-na.omit(tidy_tmp)
    
    if(length(tidy_tmp$chr) != 0)
    {
    # Vectorise
    tidy_vec<-paste0(tidy_tmp$chr,":",tidy_tmp$BP1,"-",tidy_tmp$BP2)
    
    # Return to list
    return(tidy_vec)
    } else {return(NA)}
  })
  
  # Find values that don't need modifying
  non_overlap<-calculate.overlap(outliers)
  res3<-lapply(1:length(res2),function(i){
    edit_vec<-setdiff(res2[[i]],non_overlap[[i]])
    keep_vec<-setdiff(res2[[i]],edit_vec)
    
    # How much to add on?
    if ( (((extra/window_size)/0.5)*0.5) < round((extra/window_size)/0.5)*0.5 ) {
      to_add<-(2*(round((extra/window_size)/0.5)*0.5))
    } else {
      to_add<-(2*(round((extra/window_size)/0.5)*0.5))+1
    }
    
    # Add it on
    edit2<-lapply(1:length(list(edit_vec)),function(x){
      
      # Expand out each vector
      tidy_tmp<-tidyr::separate(as.data.frame(edit_vec),col="edit_vec",c("chr","BP"),":") %>%
        tidyr::separate(col=BP,into=c("BP1","BP2"),"-")
      
      # Add/Subtract Extra for each
      tidy_tmp$BP1<-as.integer(tidy_tmp$BP1)-(to_add*window_size)
      tidy_tmp$BP2<-as.integer(tidy_tmp$BP2)+(to_add*window_size)
      tidy_tmp<-na.omit(tidy_tmp)
      
      # Vectorise
      if(length(tidy_tmp$chr) != 0)
      {
      tidy_vec<-paste0(tidy_tmp$chr,":",as.integer(tidy_tmp$BP1),"-",as.integer(tidy_tmp$BP2))
      
      # Return to list
      return(tidy_vec)
      } else {
        return(NA)
      }
    })
    
    return(sort(unlist(list(keep_vec,edit2))))
  })
  return(res3)
}

plot_cor_mat<-function(input_matrix,cor=FALSE){
  if(cor==TRUE){
    cormat <- round(cor((input_matrix)),2)
  } else {
    cormat <- round((input_matrix),2)
  }
  
  library(reshape2)
  melted_cormat <- melt(cormat)
  library(ggplot2)
  ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  upper_tri <- get_upper_tri(cormat)
  upper_tri
  
  # Melt the correlation matrix
  library(reshape2)
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Heatmap
  library(ggplot2)
  ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  # Reorder the correlation matrix
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  # Print the heatmap
  #print(ggheatmap)
  
  ggheatmap + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      axis.text = element_text(size=16),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
}

reorder_cormat <- function(cormat){
  # Mirror matrix
  cormat[lower.tri(cormat)]<-t(cormat)[lower.tri(cormat)]
  dd <- dist(cormat)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


# Function for making vectors
calculate_change<-function(data=NULL,var=NULL,rad=NULL){
  ddtmp<-data
  # Don't scale as scaling is already performed
  tmp_scale<-data.frame(val=ddtmp[,var],
                        RAD=ddtmp$RADIATION)
  vec_in<-tmp_scale[tmp_scale$RAD == rad,"val"]
  
  # Set up, data is scaled so vectors are comparable
  tmp<-data.frame(value=sort(vec_in),
                  rank=seq(1:length(vec_in)))
  # Model
  mod<-lm(tmp$value~tmp$rank)
  
  # Subtract top from bottom
  change<-max(predict(mod))-min(predict(mod))
  
  # Return
  return(change)
}

# Plot PCA
plot_pca<-function(PCA_out,PCA_scores){
  tmp<-data.frame(RADS=dd$RADIATION,
                  PC1=PCA_scores[,1],
                  PC2=PCA_scores[,2])
  
  var1<-round(summary(PCA_out)[]$importance[2,1]*100,1)
  var2<-round(summary(PCA_out)[]$importance[2,2]*100,1)
  
g1<-ggplot(tmp,aes(x=PC1,y=PC2,colour=RADS,shape=RADS))+
  geom_point(size=3)+
  stat_ellipse(level=0.95,show.legend = F)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(size=16),
        legend.text = element_text(size=16),
        legend.position = "right")+
  labs(colour="Radiation",shape="Radiation",x=paste0("PC1 (",var1," %)"),y=paste0("PC2 (",var2," %)"))+
         scale_colour_manual(values = c("red2","blue2","gold2","forestgreen"))+
  scale_shape_manual(values=c(15:18))
return(g1)
}

# Plot PCA
plot_scaled<-function(PC1,PC2,scale1,scale2){
  tmp<-data.frame(RADS=dd$RADIATION,
                  PC1=PC1,
                  PC2=PC2)
  
  g1<-ggplot(tmp,aes(x=PC1,y=PC2,colour=RADS,shape=RADS))+
    geom_point(size=3)+
    stat_ellipse(level=0.95,show.legend = F)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=15),
          axis.text = element_text(size=16),
          legend.text = element_text(size=16),
          legend.position = "right",
          legend.title = element_blank())+
    labs(colour="Radiation",shape="Radiation",x=paste0("PC1 (",scale1," %)"),y=paste0("PC2 (",scale2," %)"))+
    scale_colour_manual(values = c("red2","blue2","gold2","forestgreen"))+
    scale_shape_manual(values=c(15:18))
  return(g1)
}


##### Modified version of Stuart's f.theta R code, there's calculates angles based on correlations but we've not enough data ######
f.theta.jim<-function(table.of.diff.in.means, unique.id, select.col){
  library(Morpho)
  angle.matrix.radians <- matrix(nrow = length(unique.id), ncol = length(unique.id))
  for(i in 1:length(unique.id)){
    for(j in 1:length(unique.id)){
      angle.matrix.radians[i,j] <- round(angle.calc(t(table.of.diff.in.means[i,select.col]), t(table.of.diff.in.means[j,select.col])), 3) #
    }
  }
  rownames(angle.matrix.radians) <- unique.id
  colnames(angle.matrix.radians) <- unique.id
  angle.matrix.degrees <- round(angle.matrix.radians*(180/pi), 3)
  angle.output <- list(angle.matrix.radians, angle.matrix.degrees)
  names(angle.output) <- c("theta.radians", "theta.degrees")
  return(angle.output)
}
