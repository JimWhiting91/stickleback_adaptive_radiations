######################################################
# Make Phylogeny for all marine and freshwater fish
lib <- c("adegenet","vcfR","ape","ggplot2","ggtree","phytools","parallel","poppr","phyloch")
lapply(lib,library,character.only=T)

######################################################################
# Read in VCF
vcf_in <- read.vcfR("outputs/all_marine_fresh_stacks_output_p5r05/all_marine_fresh_stacks_output_p5r05_linkage_filtered_nosex.vcf")
gen_dd <- vcfR2genind(vcf_in)

######################################################################
# Build tree and bootstrap
pop_tree <- aboot(gen_dd, tree = "bionj", distance = nei.dist, sample = 6, showtree = T, cutoff = 0, quiet = F,threads = 6,root = FALSE)
saveRDS(pop_tree,"all_marine_fresh_phylo_bootstrapped_tree.rds")
pop_tree <- readRDS("outputs/all_marine_fresh_phylo_bootstrapped_tree.rds")

# Visualise and check
pdf("test_tree.pdf",width=20,height=40)
plot.phylo(pop_tree, cex = 0.8, font = 2, adj = 0)
nodelabels(pop_tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
dev.off()

######################################################################
# Collapse low support nodes...
tre <- pop_tree
temp <- tre
myBoots <- pop_tree$node.label

N <- length(tre$tip.label)
toCollapse <- match(which(myBoots<80)+N, temp$edge[,2])
temp$edge.length[toCollapse] <- 0
tre3 <- di2multi(temp, tol=0.00001)
write.tree(tre3, file = "all_fresh_marine_phylo_bootstrapped_tree")

# Visualise tree
pdf("figs/test_tree_bootstrapped.pdf",width=30,height=30)
plot(tre3,type="unrooted")
dev.off()

######################################################################
# Plot tree
tree_to_plot <- tre3

# Plot simple to get nodes
pdf("figs/test_tree_nodes.pdf",width=100,height=100)
plot(tree_to_plot,"unrooted")
dev.off()


# Remove low quality individuals
to_drop <- c("CORC150970","HOGG150597","LYND151111","BRAN150613","KENN150671","KENN150672")
for(ind in to_drop){
tree_to_plot <- drop.tip(tree_to_plot,ind)
}

# Colour tips
popmap <- read.table("data/MxF_VCFs/all_fresh_marine_phylo_p5.popmap")
popmap <- popmap[popmap$V1 %in% tree_to_plot$tip.label,]

# Find Clades...
groups_to_plot <- c("ALASKA","BC","ICE","SCOT")
groups_list <- lapply(groups_to_plot,function(x){return(popmap[popmap$V2 == x,"V1"])})
tree_cols <- edge.color(tree_to_plot,groups_list,col=c("red2","blue2","gold2","forestgreen"),bgcol="black")

# Add tip cols...
marine_tips <- c("MUD","LICA","NYPS","OBSM")
marine_list <- lapply(marine_tips,function(x){return(popmap[grep(x,popmap$V1),"V1"])})
tip_cols <- tip.color(tree_to_plot,marine_list,col=c("red2","blue2","gold2","forestgreen"),bgcol="white")

# Plot in base R
pdf("figs/Figure1_all_marine_fresh_tree.pdf",width=20,height=20)
plot(tree_to_plot,"unrooted",
     edge.color=tree_cols,
     show.tip.label=FALSE)
tiplabels(tip=which(tip_cols != "white"),pch=19,col=tip_cols[tip_cols != "white"],cex=2)
dev.off()

pdf("figs/Figure1_all_marine_fresh_tree_with_tips.pdf",width=100,height=100)
plot(tree_to_plot,"unrooted",
     edge.color=tree_cols,
     show.tip.label=TRUE)
tiplabels(tip=which(tip_cols != "white"),pch=19,col=tip_cols[tip_cols != "white"],cex=2)
dev.off()

