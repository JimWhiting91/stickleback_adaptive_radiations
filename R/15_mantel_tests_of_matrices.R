####################################
# Analysis of how distance predicts parallelism - Matrix analyses

# Get packages
lib<-c("data.table","ggplot2","dplyr","MASS","effects","vegan","ade4","ape")
lapply(lib,library,character.only=T)

# Read in matrices
MxF_mat<-as.matrix(read.table("outputs/MxF_FST_overlap_props.txt",header=T))
env_mat<-as.matrix(read.table("outputs/Freshwater_env_distance_mat.txt",header=T))
pheno_mat<-as.matrix(read.table("outputs/Freshwater_pheno_distance_mat.txt",header=T))
parallel_mat<-as.matrix(read.table("outputs/Freshwater_parallel_distance_mat.txt",header=T))

# Get FST genetic distance mat
FST_mat<-read.csv("outputs/All_Radiations_All_Pops_Fst.csv")
FST_mat<-as.matrix(FST_mat[,2:ncol(FST_mat)])

# Change TROU to TROUT
cols<-colnames(FST_mat)
cols[cols == "TROS"]<-"TROF"
colnames(FST_mat)<-cols
rownames(FST_mat)<-cols

cols<-colnames(MxF_mat)
cols[cols == "TROS"]<-"TROF"
colnames(MxF_mat)<-cols
rownames(MxF_mat)<-cols

# Get order
pop_order<-as.character(colnames(env_mat))

# Add rownmaes
rownames(env_mat)<-pop_order
rownames(pheno_mat)<-pop_order
rownames(parallel_mat)<-pop_order

# Reorder
FST_mat<-FST_mat[pop_order,pop_order]
pheno_mat<-pheno_mat[pop_order,pop_order]
parallel_mat<-parallel_mat[pop_order,pop_order]
MxF_mat<-MxF_mat[pop_order,pop_order]

# Perform mantel tests
# First test, does environment explain marine x freshwater FST overlap
mantel(-1*env_mat, MxF_mat, method="spearman", permutations=9999,parallel=4)
mantel(-1*pheno_mat, MxF_mat, method="spearman", permutations=9999,parallel=4)
mantel(-1*FST_mat, MxF_mat, method="spearman", permutations=9999,parallel=4)

# Now do partial mantel tests whilst correcting for genetic distance
mantel.partial(-1*env_mat, MxF_mat, -1*FST_mat, method = "spearman", permutations = 9999,parallel = 4)
mantel.partial(-1*pheno_mat, MxF_mat, -1*FST_mat, method = "spearman", permutations = 9999,parallel = 4)

# 
# # Alternative mantel test method
# mantel.rtest(env_mat, MxF_mat, nrepet = 9999)
# 
# # Alteranative tests
# mantel.test(env_mat, MxF_mat, method="pearson", permutations=9999,graph=TRUE)
# mantel(pheno_mat, MxF_mat, method="pearson", permutations=9999,parallel=4)
# mantel(FST_mat, MxF_mat, method="pearson", permutations=9999,parallel=4)



