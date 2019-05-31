##### This script simulates vectors of dummy variables to examine expected vectors of random variables ######
# Stuart Functions
source("R/Stuart_Nature_Rcode_functions.R")
source("R/Genomic_Parallelism_Rcode_functions.R")

# Load packages
lib<-as.vector(c("pbapply","ggpubr","ggrepel","data.table","VennDiagram","ggplot2","parallel","dplyr"))
lapply(lib,library,character.only=TRUE)

# Read in data
dd<-read.csv("data/All_Lakes_Env_and_Pheno_v2.csv")
colnames(dd)[1]<-"RADIATION"
envs<-c("Ca","Gyro","Na","pH","Schisto","Zn")

# Read in armour_data
armour_vec<-c("Avg.resDS1","Avg.resDS2","Avg.resPS",
              "Avg.resLP","Avg.resHP","Avg.resBAP","Avg.n.plate")

# Read in shape
shape_dd<-read.table("data/phenotype_data/body_shape_residuals_separate_radiations_and_marine.txt",header=T)
shape_dd<-shape_dd[shape_dd$system == "freshwater",]
shape_vec<-colnames(shape_dd)[c(7:32)]

# Redesign
shape_avg<-data.frame(shape_dd %>% group_by(LAKE) %>% summarise_all(funs(mean)))
shape_avg<-shape_avg[order(dd$LAKE),]

# Add to dd
shape_mat<-matrix(nrow = nrow(dd),ncol=length(shape_vec))
for(i in 1:nrow(dd)){
  for(j in 1:length(shape_vec)){
    shape_mat[i,j]<-shape_avg[as.character(shape_avg$LAKE)==as.character(dd$LAKE[i]),shape_vec[j]]
}
}
dd_all<-cbind(dd,data.frame(shape_mat))
colnames(dd_all)[44:69]<-shape_vec

##### Set up each analysis as a list ######
group_list<-list(envs,shape_vec,armour_vec)
group_vec<-c("ENV","SHAPE","ARMOUR")

###### What are our null expectations of vectors? ######

# To do this, we simulate random variables and calculate the vectors through them
rads<-c("ALASKA","BC","ICELAND","SCOTLAND")
rad_vec<-rads
iterations<-1000

# Check if file exists
if(file.exists(paste0("outputs/05_Null_vector_simulations_",iterations,".txt"))=="FALSE"){
simmed_vectors<-data.frame(rbindlist(lapply(1:iterations,function(iter){
  set.seed(iter)
  change_out<-lapply(1:length(group_list),function(x){
    
    # Get variables
    tmp_dd<-dd_all[,c(group_list[[x]],"RADIATION")]
    
    # Now sim distributions
    sims<-lapply(1:length(group_list[[x]]),function(y){
      out<-unlist(lapply(1:length(rads),function(RAD){
        mean_tmp<-mean(tmp_dd[tmp_dd$RADIATION == rads[RAD],group_list[[x]][y]])
        sd_tmp<-sd(tmp_dd[tmp_dd$RADIATION == rads[RAD],group_list[[x]][y]])
        simmed<-rnorm(nrow(tmp_dd[tmp_dd$RADIATION == rads[RAD],]),
                      mean_tmp,sd_tmp)
        return(simmed)
      }))
      return(out)
    })
    
    # Set up sim data
    dd_sim1<-data.frame(do.call("cbind",sims))
    colnames(dd_sim1)<-group_list[[x]]
    dd_sim1$RADIATION<-c(rep("ALASKA",19),
                        rep("BC",18),
                        rep("ICELAND",18),
                        rep("SCOTLAND",18))
    
    # Now we condense down to a PCA
    if(group_vec[x] == "SHAPE"){
    sim_pca<-prcomp(dd_sim1[,group_list[[x]]],center=T,scale=F)
    } else {
      sim_pca<-prcomp(dd_sim1[,group_list[[x]]],center=T,scale=T)
    }
    
    # Get PC proportions
    PC_props1<-(sim_pca[[1]]^2)/sum(sim_pca[[1]]^2)
    PC_props<-PC_props1[PC_props1 > 0.1]
    if(length(PC_props) < 3){
      PC_props<-PC_props1[1:3]
    }
    dd_sim<-data.frame(RADIATION=dd_sim1$RADIATION,
                       sim_pca$x[,1:length(PC_props)])
    colnames(dd_sim)[2:(length(PC_props)+1)]<-paste0("PC",seq(1,length(PC_props),1))
    
    # Scale PCs
    for(i in 1:length(PC_props)){
      dd_sim[,i+1]<-dd_sim[,i+1]*(PC_props[i]*100)
    }
    
    # Set up tmp_var
    tmp_var<-paste0("PC",seq(1,length(PC_props),1))
    
    # Across tmp_var calculate vectors for all radiations
    out<-data.frame(rbindlist(lapply(1:length(tmp_var),function(y){
      A_tmp<-data.frame(vec=calculate_change(dd_sim,tmp_var[y],"ALASKA"),
                        rad="ALASKA",
                        var=tmp_var[y])
      B_tmp<-data.frame(vec=calculate_change(dd_sim,tmp_var[y],"BC"),
                        rad="BC",
                        var=tmp_var[y])
      I_tmp<-data.frame(vec=calculate_change(dd_sim,tmp_var[y],"ICELAND"),
                        rad="ICELAND",
                        var=tmp_var[y])
      N_tmp<-data.frame(vec=calculate_change(dd_sim,tmp_var[y],"SCOTLAND"),
                        rad="SCOTLAND",
                        var=tmp_var[y])
      return(rbind(A_tmp,B_tmp,I_tmp,N_tmp))
    })))
    
    # Add label
    out$group<-group_vec[x]
    
    # Return
    out
  })
  
  # Calculate angles
  
  angles_out<-lapply(1:length(group_list),function(x){
    
    # List out vectors
    vec_dd<-t(data.frame(lapply(1:length(rad_vec),function(y){
      return(change_out[[x]][change_out[[x]]$rad==rad_vec[y],"vec"])
    })))
    
    avg_vec_diffs<-data.frame(vec_dd)
    colnames(avg_vec_diffs)<-paste0("PC",seq(1,ncol(vec_dd),1))
    
    # Stuart calc functions
    # Use Stuart et al function to calculate theta and L
    rad.theta <- f.theta.jim(avg_vec_diffs, unique.id = rad_vec, select.col =1:length(avg_vec_diffs))
    rad.theta.radians <- rad.theta[[1]]
    rad.theta.degrees <- rad.theta[[2]]
    
    rad.theta.degrees_vec<-rad.theta.degrees[upper.tri(rad.theta.degrees, diag = FALSE)]
    
    rad.deltaL <- f.deltaL(avg_vec_diffs, unique.id = rad_vec, select.col =1:length(avg_vec_diffs[1,]))
    rad.theta.deltaL_vec<-abs(rad.deltaL[[1]][upper.tri(rad.deltaL[[1]], diag = FALSE)])
    
    # Combine the two
    out_mat<-matrix(nrow = 2,ncol = length(combn(rad_vec,2)[1,]))
    out_mat[1,]<-rad.theta.degrees_vec
    out_mat[2,]<-rad.theta.deltaL_vec
    out_dd<-data.frame(out_mat,
                       measure=c("Theta_Deg","deltaL"))
    
    # Rename columns
    colnames(out_dd)[1:6]<-c("Alaska_BC",
                             "Alaska_Iceland",
                             "BC_Iceland",
                             "Alaska_Scotland",
                             "BC_Scotland",
                             "Iceland_Scotland")
    
    # Reorder columns alphabetically
    out_dd<-out_dd[,order(colnames(out_dd))]
    
    # Give rad Info
    out_dd$group<-group_vec[x]
    
    # Also get length of vectors
    out_dd_L<-rad.deltaL[[2]]
    colnames(out_dd_L)<-c("Vector_L","Radiation")
    out_dd_L$group<-group_vec[x]
    
    # Return from CHR
    return(list(out_dd,out_dd_L,rad.theta.radians,rad.deltaL[[1]]))
  })
  
  # Extract the outputs
  # out_dd
  out_dd<-data.frame(rbindlist(lapply(1:length(angles_out),function(x){
    return(angles_out[[x]][[1]])})))
  
  # Modify out_dd so that we only have angles less than 90
  for(i in 1:6){
    tmp<-data.frame(new=180-out_dd[out_dd$measure == "Theta_Deg",i],
                    old=out_dd[out_dd$measure == "Theta_Deg",i])
    out_dd[out_dd$measure == "Theta_Deg",i]<-do.call(pmin,tmp)
  }
  
  # Transform out_dd and export for modelling
  out_dd_2<-data.frame(t(out_dd))
  colnames(out_dd_2)<-c("ENV_theta","ENV_deltaL",
                        "SHAPE_theta","SHAPE_deltaL",
                        "ARMOUR_theta","ARMOUR_deltaL")
  out_dd_2<-out_dd_2[1:6,]
  out_dd_2$rads<-rownames(out_dd_2)
  
  return(out_dd_2)
})))
write.table(simmed_vectors,
            paste0("outputs/05_Null_vector_simulations_",iterations,".txt"),
            row.names = F,quote = F,sep="\t")
}

# Read in results if needed
simmed_vectors<-read.table(paste0("outputs/05_Null_vector_simulations_",iterations,".txt"),
                           header=T)

# Read back in observed data and estimate p-vals
out_dd_2<-read.table("outputs/05_Observed_theta_delta_vectors.txt",header=T)
out_dd_2<-out_dd_2[c(1:6,9)]

# Get quantiles
comps<-unique(simmed_vectors$rads)
p_cutoffs<-data.frame(rbindlist(lapply(1:6,function(x){
  tmp<-out_dd_2[,c(7,x)]
  for(i in 1:nrow(tmp)){
    tmp$p_val[i]<-(length(simmed_vectors[simmed_vectors$rads == tmp$rads[i] &
                                           as.numeric(as.character(simmed_vectors[,x])) < as.numeric(as.character(tmp[i,2])),1]))/iterations
  }
  colnames(tmp)<-c("Rads","Measure")
  tmp$Variable<-colnames(out_dd_2)[x]
  return(tmp)
})))

#########################################################################
# We run Gill rakers separately, without PCs
# To do this, we simulate random variables and calculate the vectors through them
rads<-c("ALASKA","BC","ICELAND","SCOTLAND")
rad_vec<-rads
iterations<-1000

gill_vec<-c("Gill.Raker.N","Resid.Raker.L")

# Check if file exists
if(file.exists(paste0("outputs/05_Null_vector_simulations_GILL_",iterations,".txt"))=="FALSE"){
  simmed_vectors<-data.frame(rbindlist(lapply(1:iterations,function(iter){
    set.seed(iter)
      
      # Get variables
      tmp_dd<-dd_all[,c(gill_vec,"RADIATION")]
      # Scale
      tmp_dd$Gill.Raker.N<-scale(tmp_dd$Gill.Raker.N)
      tmp_dd$Resid.Raker.L<-scale(tmp_dd$Resid.Raker.L)
      
      # Now sim distributions
      sims<-lapply(1:length(gill_vec),function(y){
        out<-unlist(lapply(1:length(rads),function(RAD){
          mean_tmp<-mean(tmp_dd[tmp_dd$RADIATION == rads[RAD],gill_vec[y]])
          sd_tmp<-sd(tmp_dd[tmp_dd$RADIATION == rads[RAD],gill_vec[y]])
          simmed<-rnorm(nrow(tmp_dd[tmp_dd$RADIATION == rads[RAD],]),
                        mean_tmp,sd_tmp)
          return(simmed)
        }))
        return(out)
      })
      
      # Set up sim data
      dd_sim1<-data.frame(do.call("cbind",sims))
      colnames(dd_sim1)<-gill_vec
      dd_sim1$RADIATION<-c(rep("ALASKA",19),
                           rep("BC",18),
                           rep("ICELAND",18),
                           rep("SCOTLAND",18))
      
      # Across tmp_var calculate vectors for all radiations
      out<-data.frame(rbindlist(lapply(1:length(gill_vec),function(y){
        A_tmp<-data.frame(vec=calculate_change(dd_sim1,gill_vec[y],"ALASKA"),
                          rad="ALASKA",
                          var=gill_vec[y])
        B_tmp<-data.frame(vec=calculate_change(dd_sim1,gill_vec[y],"BC"),
                          rad="BC",
                          var=gill_vec[y])
        I_tmp<-data.frame(vec=calculate_change(dd_sim1,gill_vec[y],"ICELAND"),
                          rad="ICELAND",
                          var=gill_vec[y])
        N_tmp<-data.frame(vec=calculate_change(dd_sim1,gill_vec[y],"SCOTLAND"),
                          rad="SCOTLAND",
                          var=gill_vec[y])
        return(rbind(A_tmp,B_tmp,I_tmp,N_tmp))
      })))
      
      # Add label
      out$group<-"GILL"
    
    # Calculate angles
      # List out vectors
      vec_dd<-t(data.frame(lapply(1:length(rad_vec),function(y){
        return(out[out$rad==rad_vec[y],"vec"])
      })))
      
      avg_vec_diffs<-data.frame(vec_dd)
      colnames(avg_vec_diffs)<-paste0("PC",seq(1,ncol(vec_dd),1))
      
      # Stuart calc functions
      # Use Stuart et al function to calculate theta and L
      rad.theta <- f.theta.jim(avg_vec_diffs, unique.id = rad_vec, select.col =1:length(avg_vec_diffs))
      rad.theta.radians <- rad.theta[[1]]
      rad.theta.degrees <- rad.theta[[2]]
      
      rad.theta.degrees_vec<-rad.theta.degrees[upper.tri(rad.theta.degrees, diag = FALSE)]
      
      rad.deltaL <- f.deltaL(avg_vec_diffs, unique.id = rad_vec, select.col =1:length(avg_vec_diffs[1,]))
      rad.theta.deltaL_vec<-abs(rad.deltaL[[1]][upper.tri(rad.deltaL[[1]], diag = FALSE)])
      
      # Combine the two
      out_mat<-matrix(nrow = 2,ncol = length(combn(rad_vec,2)[1,]))
      out_mat[1,]<-rad.theta.degrees_vec
      out_mat[2,]<-rad.theta.deltaL_vec
      out_dd<-data.frame(out_mat,
                         measure=c("Theta_Deg","deltaL"))
      
      # Rename columns
      colnames(out_dd)[1:6]<-c("Alaska_BC",
                               "Alaska_Iceland",
                               "BC_Iceland",
                               "Alaska_Scotland",
                               "BC_Scotland",
                               "Iceland_Scotland")
      
      # Reorder columns alphabetically
      out_dd<-out_dd[,order(colnames(out_dd))]
      
      # Give rad Info
      out_dd$group<-group_vec[x]
      
      # Also get length of vectors
      out_dd_L<-rad.deltaL[[2]]
      colnames(out_dd_L)<-c("Vector_L","Radiation")
      out_dd_L$group<-group_vec[x]
      
    # Modify out_dd so that we only have angles less than 90
    for(i in 1:6){
      tmp<-data.frame(new=180-out_dd[out_dd$measure == "Theta_Deg",i],
                      old=out_dd[out_dd$measure == "Theta_Deg",i])
      out_dd[out_dd$measure == "Theta_Deg",i]<-do.call(pmin,tmp)
    }
    
    # Transform out_dd and export for modelling
    out_dd_2<-data.frame(t(out_dd[,1:6]))
    colnames(out_dd_2)<-c("GILL_theta","GILL_deltaL")
    out_dd_2<-out_dd_2[1:6,]
    out_dd_2$rads<-rownames(out_dd_2)
    
    return(out_dd_2)
  })))
  write.table(simmed_vectors,
              paste0("outputs/05_Null_vector_simulations_GILL_",iterations,".txt"),
              row.names = F,quote = F,sep="\t")
}

# Read in results if needed
simmed_vectors<-read.table(paste0("outputs/05_Null_vector_simulations_GILL_",iterations,".txt"),
                           header=T)

# Read back in observed data and estimate p-vals
out_dd_2<-read.table("outputs/05_Observed_theta_delta_vectors.txt",header=T)
out_dd_2<-out_dd_2[c(7:9)]

# Get quantiles
comps<-unique(simmed_vectors$rads)
p_cutoffs<-data.frame(rbindlist(lapply(1:2,function(x){
  tmp<-out_dd_2[,c(3,x)]
  for(i in 1:nrow(tmp)){
    tmp$p_val[i]<-(length(simmed_vectors[simmed_vectors$rads == tmp$rads[i] &
                                           as.numeric(as.character(simmed_vectors[,x])) < as.numeric(as.character(tmp[i,2])),1]))/iterations
  }
  colnames(tmp)<-c("Rads","Measure")
  tmp$Variable<-colnames(out_dd_2)[x]
  return(tmp)
})))



