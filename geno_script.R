# SSR GENOTYPING PRACTICAL – SESSION 3
  
  ## SMC, Clustering, Phylogeny, Ordination & Mantel Test
  
  ## **SECTION 0: Package Installation and Loading**
  
  
# 0.1 Install required packages (run once)
install.packages(c(
  "poppr", "adegenet", "ape", "vegan", "ggplot2",
  "dendextend", "cluster", "pvclust", "NbClust",
  "factoextra", "rlang", "ggpubr", "purrr"
))

# 0.2 Load libraries
library(poppr)
library(adegenet)
library(ape)
library(vegan)
library(ggplot2)
library(dendextend)
library(cluster)
library(pvclust)
library(NbClust)
library(factoextra)

## **SECTION 1: Import SSR Genotype Data**

# 1.1 Read SSR genotype data
ssr <- read.csv("H_SSR.csv", row.names = 1, na.strings = c("", "NA"))

# 1.2 Inspect data
head(ssr)

## **SECTION 2: Convert SSR Data to Genetic Object**
  
  
# 2.1 Convert to genind object (diploid, codominant SSRs)
gen <- df2genind(
  ssr,
  ploidy = 2,
  sep = "_",
  NA.char = NA,
  type = "codom"
)

# 2.2 Summary of genind object
summary(gen)

## **SECTION 3: Locus-Level Summary Statistics**
  
  
# 3.1 Capture summary text
gen_summary_text <- capture.output(summary(gen))

# 3.2 Extract number of alleles per locus
num_alleles_line <- grep("Number of alleles per locus", gen_summary_text, value = TRUE)
num_alleles_values <- as.numeric(unlist(strsplit(gsub(".*: ", "", num_alleles_line), " ")))

# 3.3 Extract observed heterozygosity
obs_het_line <- grep("Observed heterozygosity", gen_summary_text, value = TRUE)
obs_het_values <- as.numeric(unlist(strsplit(gsub(".*: ", "", obs_het_line), " ")))

# 3.4 Extract expected heterozygosity
exp_het_line <- grep("Expected heterozygosity", gen_summary_text, value = TRUE)
exp_het_values <- as.numeric(unlist(strsplit(gsub(".*: ", "", exp_het_line), " ")))

# 3.5 Compile locus summary table
locus_summary <- data.frame(
  Locus = locNames(gen),
  NumAlleles = num_alleles_values,
  ObsHet = obs_het_values,
  ExpHet = exp_het_values
)

# 3.6 Export locus summary
write.csv(locus_summary, "SSR_Locus_Summary_exact.csv", row.names = FALSE)

## **SECTION 4: Filtering Monomorphic Loci**
  
  
# 4.1 Identify loci with no allelic variation
locus_var <- sapply(locNames(gen), function(locus) {
  length(unique(alleles(gen, locus)))
})

monomorphic_loci <- names(locus_var[locus_var <= 1])

# 4.2 Remove monomorphic loci
gen_var <- if (length(monomorphic_loci) > 0) {
  gen[, -which(locNames(gen) %in% monomorphic_loci)]
} else {
  gen
}

## **SECTION 5: Simple Matching Coefficient (SMC)**
  
  
# 5.1 Convert genind to allele matrix
allele_matrix <- tab(gen_var, NA.method = "mean")
n_ind <- nInd(gen_var)

# 5.2 Compute SMC similarity matrix
smc <- matrix(
  0,
  nrow = n_ind,
  ncol = n_ind,
  dimnames = list(indNames(gen_var), indNames(gen_var))
)

for (i in 1:n_ind) {
  for (j in i:n_ind) {
    match_count <- sum(allele_matrix[i,] == allele_matrix[j,], na.rm = TRUE)
    total_count <- sum(!is.na(allele_matrix[i,]) & !is.na(allele_matrix[j,]))
    smc[i, j] <- match_count / total_count
    smc[j, i] <- smc[i, j]
  }
}

diag(smc) <- 1

# 5.3 Convert similarity to distance
smc_dist <- as.dist(1 - smc)

# 5.4 Export matrices
write.csv(as.matrix(smc), "SSR_SMC_similarity.csv", row.names = TRUE)

  
## **SECTION 6: Hierarchical Clustering (UPGMA)**
  
  
# 6.1 Perform UPGMA clustering
hclust_obj <- hclust(smc_dist, method = "average")

## **SECTION 7: Optimal Number of Clusters**
  
### **7.1 Elbow Method**
  
  
dist_mat <- as.matrix(smc_dist)

wss <- sapply(2:10, function(k) {
  clust <- cutree(hclust_obj, k)
  sum(sapply(unique(clust), function(i) {
    members <- which(clust == i)
    if (length(members) == 1) return(0)
    sum(dist_mat[members, members]^2) / 2
  }))
})

png("Elbow_plot.png", 800, 600)
plot(2:10, wss, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Within-cluster sum of squares",
     main = "Elbow Method")
dev.off()

### **7.2 Silhouette Method**
sil_width <- sapply(2:10, function(k) {
  clust <- cutree(hclust_obj, k)
  ss <- silhouette(clust, smc_dist)
  mean(ss[,3])
})

png("Silhouette_plot.png", 800, 600)
plot(2:10, sil_width, type = "b", pch = 19,
     xlab = "k",
     ylab = "Average silhouette width",
     main = "Silhouette Method")
dev.off()

  
### **7.3 Gap Statistic**
set.seed(123)
gap_stat <- clusGap(as.matrix(dist_mat), FUN = hcut, K.max = 10, B = 100)

png("Gap_statistic_plot.png", 800, 600)
plot(gap_stat, main = "Gap Statistic")
dev.off()

### **7.4 Final Cluster Assignment**
  
K <- 4 ## replace with optimal
clusters <- cutree(hclust_obj, k = K)

write.csv(
  data.frame(Sample = names(clusters), Cluster = clusters),
  "SSR_clusters_optimal.csv",
  row.names = FALSE
)

## **SECTION 8: Phylogenetic Analysis**
  
## **8.1 UPGMA with Bootstrap**
  
  
set.seed(123)
pv <- pvclust(
  smc,
  method.hclust = "average",
  method.dist = "euclidean",
  nboot = 1000
)

plot(pv)
pvrect(pv, alpha = 0.95)

write.tree(as.phylo(pv$hclust), "SSR_UPGMA_1000boot.nwk")

### **8.2 Neighbor-Joining Tree**
  
  
nj_tree <- nj(smc_dist)

boot_nj <- boot.phylo(
  phy = nj_tree,
  x = smc,
  FUN = function(xx) nj(dist(xx)),
  B = 1000
)

plot(nj_tree)
nodelabels(boot_nj, frame = "n", cex = 0.7)

write.tree(nj_tree, "SSR_NJ_1000boot.nwk")

## **SECTION 9: Principal Coordinate Analysis (PCoA)**
  
  
# 9.1 Run PCoA
pcoa <- cmdscale(smc_dist, eig = TRUE, k = 10)

# 9.2 Variance explained
var_exp <- 100 * pcoa$eig[1:10] / sum(pcoa$eig)
png("PCoA_scree_plot.png", width = 800, height = 600)
barplot(
  var_exp[1:10],
  names.arg = paste0("PC", 1:10),
  ylim = c(0, 10),                      # <-- FIXED Y-AXIS
  ylab = "Percentage of Variance Explained",
  xlab = "PCoA Axes",
  main = "PCoA Scree Plot (First 10 Axes)",
  col = "grey70"
)
dev.off()

dev.off()
# 9.3 Export scree data
write.csv(
  data.frame(PC = paste0("PC", 1:10), VariancePercent = var_exp),
  "PCoA_ScreePlot_First10PCs.csv",
  row.names = FALSE
)

# 9.4 Export individual coordinates
pc_coords <- as.data.frame(pcoa$points[,1:10])
pc_coords$Sample <- rownames(pc_coords)
write.csv(pc_coords, "PCoA_IndividualCoordinates_First10PCs.csv", row.names = FALSE)

## **SECTION 10: PCoA Visualization**
  
### **10.1 Cluster-based PCoA**
png("PCoA_clusters.png", 800, 600)
plot(
  pcoa$points[,1], pcoa$points[,2],
  col = as.factor(clusters), pch = 19,
  xlab = paste0("PCoA1 (", round(var_exp[1], 2), "%)"),
  ylab = paste0("PCoA2 (", round(var_exp[2], 2), "%)"),
  main = "PCoA by Genetic Clusters"
)
legend("topright",
       legend = levels(as.factor(clusters)),
       col = 1:length(unique(clusters)),
       pch = 19,
       title = "Cluster")
dev.off()


### **10.2 Population-based PCoA**
pop_info <- read.csv("H_Pop.csv", row.names = 1)
pop <- as.factor(pop_info[rownames(smc),1])

png("PCoA_population.png", 800, 600)
plot(
  pcoa$points[,1], pcoa$points[,2],
  col = pop, pch = 19,
  xlab = paste0("PCoA1 (", round(var_exp[1], 2), "%)"),
  ylab = paste0("PCoA2 (", round(var_exp[2], 2), "%)"),
  main = "PCoA by Population"
)
legend("topright",
       legend = levels(pop),
       col = 1:length(levels(pop)),
       pch = 19,
       title = "Population")
dev.off()


## **SECTION 11: Mantel Test (Genotype–Phenotype Association)**

# 11.1 Read phenotype data
pheno <- read.csv("H_Pheno.csv", row.names = 1)

# 11.2 Compute distance matrices
pheno_dist <- dist(scale(pheno), method = "euclidean")
smc_dist <- as.dist(1 - smc)

# 11.3 Mantel test
mantel_test <- mantel(pheno_dist, smc_dist, permutations = 999)
print(mantel_test)

# 11.4 Mantel scatter plot
png("Mantel_plot.png", 800, 600)
plot(as.vector(smc_dist), as.vector(pheno_dist),
     pch = 19, col = "steelblue",
     xlab = "Genetic Distance (SSR)",
     ylab = "Phenotypic Distance")
abline(lm(as.vector(pheno_dist) ~ as.vector(smc_dist)),
       col = "red", lty = 2)
dev.off()


# **FINAL OUTPUTS (NUMBERED CLEARLY)**
  
#1. **SSR_Locus_Summary_exact.csv** – locus diversity statistics
#2. **SSR_SMC_similarity.csv** – genetic similarity matrix
#3. **SSR_SMC_distance.csv** – genetic distance matrix
#4. **SSR_clusters_optimal.csv** – cluster membership
#5. **Silhouette_plot.png** – silhouette diagnostics
#6. **Gap_statistic_plot.png** – gap statistic
#7. **SSR_UPGMA_1000boot.nwk** – UPGMA tree
#8. **SSR_NJ_1000boot.nwk** – NJ tree
#9. **PCoA_ScreePlot_First10PCs.csv** – variance explained
#10. **PCoA_IndividualCoordinates_First10PCs.csv** – PCoA coordinates
#11. **PCoA_clusters.png** – PCoA by clusters
#12. **PCoA_population.png** – PCoA by population
#13. **Mantel_plot.png** – genotype–phenotype relationship

