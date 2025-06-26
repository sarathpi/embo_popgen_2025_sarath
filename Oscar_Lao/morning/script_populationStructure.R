library(smacof);
library("genio");
library("BEDMatrix");
library("vegan");
library(LEA);
library(ggplot2)

rm(list=ls());

#########################################################################
# 1. Admixture model
# Goal: to understand from a coalescence point of view what admixture is
# Exercise:
# Create function called model.admixed.pop with a fastSimcoal2 model that, 
# given two populations A and B, create a new Admix pop by admixture with 
# proportions pa and 1-pa.

folder.fastSimcoal2 <- "C:\\Users\\u9424\\OneDrive - Universitat Pompeu Fabra\\Grants\\2025\\EMBO Naples June\\Second day\\";

setwd(folder.fastSimcoal2);

model.admixed.pop <- function(sample_1, sample_2, sample_admix, effective_population_size_1, effective_population_size_2, effective_population_size_admixed, time_admixture, percentage_from_population_1, time_split, effective_population_size_ancestral, number_of_blocks)
{
  lines <- c(
    "//Number of population samples (demes)",
    "3",
    "//Population effective sizes (number of genes)",
    effective_population_size_1,
    effective_population_size_2,
    effective_population_size_admixed,
    "//Sample sizes",
    sample_1,
    sample_admix,
    sample_2,
    "//Growth rates\t: negative growth implies population expansion",
    "0",
    "0",
    "0",
    "//Number of migration matrices : 0 implies no migration between demes",
    "0",
    "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
    "3 historical event",
    paste(time_admixture,"1 0",percentage_from_population_1,"1 0 0",sep = " "),
    paste(time_admixture,"1 2 1 1 0 0",sep = " "),
    paste(time_split,"0 2 1",effective_population_size_ancestral/effective_population_size_2, "0 0",sep = " "),
    "//Number of independent loci [chromosome]",
    paste(number_of_blocks,"1",sep=" ")
  )
  
  # Repeating the block 22 times
  block <- c(
    "//Per chromosome: Number of linkage blocks",
    "1",
    "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
    "DNA 1000000 1.0E-8 1.855284327902964E-8 0.0"
  )
  
  # Append number_of_blocks blocks
  for (i in 1:number_of_blocks) {
    lines <- c(lines, block)
  }
  
  # Write to file
  writeLines(lines, "DemographicModelSplitR.par")
  
  exe <- ".\\fsc28.exe"
  args <- c("-i", ".\\DemographicModelSplitR.par", "-x", "-s0", "-d", "-n", "1", "-q", "-G")
  # Execute
  system2(exe, args = args)
  
  data.t <- read.table(file=paste(folder.fastSimcoal2,"\\DemographicModelSplitR\\DemographicModelSplitR_1_1.gen", sep=""), header = T);
  
  # First four columns are snp info
  # haplotype matrix. Rows are haplotypes, columns are positions
  
  H <- t(as.matrix(data.t[,5:ncol(data.t)]));
  # First 100 haplotypes is pop 1, next pop 2
  rownames(H) <- c(rep("A",sample_1),rep("Admix",sample_admix),rep("B",sample_2));
  return(list(position = data.t[,1:4],haplotype_matrix = H));
}
#########################################################################
# 2. Project in two dimensions the relationships of a dataset generated with the model.admixed.pop
# Goal: to understand how MDS works in practice
# Exercise:
# Create a dataset of 100 chromosomes in each population. The effective population size of pop1 is 1000, the same for pop 2 and pop 3.
# the time of admixture is set to 10 generations and the time of split between the two parental populations to 10000 generations ago.
# the effective population size ancestral is 2000.
# the percentage of admixture from population 1 is 0.3.
# run the model. Generate a distance matrix and run a classical multidimensional scaling (cmdscale).
# Plot the first two dimensions.

sample_1 <- 100;
sample_2 <-100;
sample_admix <- 100;
effective_population_size_1 <- 1000;
effective_population_size_2 <- 1000;
effective_population_size_admixed <- 1000;
time_admixture <- 10;
time_split <- 10000;
effective_population_size_ancestral <- 2000;
number_of_blocks <- 10;
percentage_from_population_1 <- 0.3;
sim <- model.admixed.pop(sample_1, sample_2, sample_admix, effective_population_size_1, effective_population_size_2, effective_population_size_admixed, time_admixture, percentage_from_population_1, time_split, effective_population_size_ancestral, number_of_blocks)
# first 5 columns are information of the SNPs
haplotypes <- sim$haplotype_matrix;
# Genotypes are computed by collapsing two contiguous haplotypes
genotypes <- haplotypes[seq(from=1,to=nrow(haplotypes),by=2),] + haplotypes[seq(from=2,to=nrow(haplotypes),by=2),];
# Provide labels
pop <- rownames(haplotypes)[seq(from=1,to=nrow(haplotypes),by=2)];
rownames(genotypes) <- pop;
# Compute the standard deviation of each marker to decide if it is worth being included
sd.snps <- apply(genotypes,2,sd);
consider.snps <- which(sd.snps > 0.1);
genotypes.clean <- genotypes[,consider.snps];
# scale each marker by substracting the mean and dividing by the standard deviation
genotypes.clean.scaled <- scale(genotypes.clean);
# compute the euclidean distance between each two individuals
dmatrix <- dist(genotypes.clean.scaled);
# conduct the classical multidimensional scaling
mds.result <- cmdscale(dmatrix, k = (nrow(as.matrix(dmatrix))-1),add=T,eig=T);
points.mds <- mds.result$points[,1:2];
rownames(points.mds) <- pop;
plot(mds.result$eig); 
# plot the result (nicer with ggplot2, of course)
plot(mds.result$points[,1:2], xlab = paste("DIM 1 (",round(100*mds.result$eig[1]/sum(mds.result$eig),2) ," %)"), ylab = paste("DIM 2 (",round(100*mds.result$eig[2]/sum(mds.result$eig),2) ," %)"))
points(mds.result$points[pop=="A",1], mds.result$points[pop=="A",2], pch = 19, col ="red");
points(mds.result$points[pop=="B",1], mds.result$points[pop=="B",2], pch = 19, col ="blue");

#########################################################################
# 3. Basic admixture estimator
# Goal: to understand how the dimensions of the MDS relate to admixture parameters
# Exercise:
# Create a function that takes the first dimension of the MDS eigenvector output.
# Rescales the values of the Admix individuals to (x-meanPop1)/(meanPop2-meanPop1)

extract.admixture <- function(points.mds)
{
  popA <- which(rownames(points.mds)=="A");
  popB <- which(rownames(points.mds)=="B");
  popAd <- which(rownames(points.mds)=="Admix");
  meanA <- mean(points.mds[popA,1]);
  meanB <- mean(points.mds[popB,1]);
  return(1-(points.mds[popAd,1]-meanA)/(meanB-meanA));
}

#########################################################################
# 4. Admixture proportions and position in the first dimension of the scaled MDS
# Goal: to understand how the dimensions of the MDS relate to admixture parameters
# Exercise:
# Generate 1000 samples, changing at random the percentage of admixture from population 1

sample_1 <- 100;
sample_2 <-100;
sample_admix <- 100;
effective_population_size_1 <- 1000;
effective_population_size_2 <- 1000;
effective_population_size_admixed <- 1000;
time_admixture <- 10;
time_split <- 10000;
effective_population_size_ancestral <- 2000;
number_of_blocks <- 10;

percentage.admixture <- rep(NA,100);
percentage.admixture.estimated <- rep(NA,100);
for(rep in 1:100)
{
  percentage_from_population_1 <- runif(1,0,1);
  percentage.admixture[rep] <- percentage_from_population_1;
  sim <- model.admixed.pop(sample_1, sample_2, sample_admix,effective_population_size_1, effective_population_size_2, effective_population_size_admixed, time_admixture, percentage_from_population_1, time_split, effective_population_size_ancestral, number_of_blocks)
  # first 5 columns are information of the SNPs
  haplotypes <- sim$haplotype_matrix;
  # Genotypes are computed by collapsing two contiguous haplotypes
  genotypes <- haplotypes[seq(from=1,to=nrow(haplotypes),by=2),] + haplotypes[seq(from=2,to=nrow(haplotypes),by=2),];
  # Provide labels
  pop <- rownames(haplotypes)[seq(from=1,to=nrow(haplotypes),by=2)];
  rownames(genotypes) <- pop;
  # Compute the standard deviation of each marker to decide if it is worth being included
  sd.snps <- apply(genotypes,2,sd);
  consider.snps <- which(sd.snps > 0.1);
  genotypes.clean <- genotypes[,consider.snps];
  # scale each marker by substracting the mean and dividing by the standard deviation
  genotypes.clean.scaled <- scale(genotypes.clean);
  # compute the euclidean distance between each two individuals
  dmatrix <- dist(genotypes.clean.scaled);
  # conduct the classical multidimensional scaling
  mds.result <- cmdscale(dmatrix, k = (nrow(as.matrix(dmatrix))-1),add=T,eig=T);
  points.mds <- mds.result$points[,1:2];
  rownames(points.mds) <- pop;
  percentage.admixture.estimated[rep] <- median(extract.admixture(points.mds))  
}

plot(percentage.admixture, percentage.admixture.estimated);
summary(lm(percentage.admixture.estimated ~ percentage.admixture));

# What happens if the time of admixture is much older?

#########################################################################
# 5. Variance in the admixture proportions and time of admixture
# Goal: to understand how the dimensions of the MDS relate to time of admixture parameters
# Exercise:
# Generate 100 samples, with a fixed percentage of admixture from population 1 of 0.5
# change the time of admixture between 1 generation to 100 generations
# estimate the variance in the percentage of estimated admixture
# plot the time  of admixture vs the variance in the percentage of estimated admixture

effective_population_size_1 <- 1000;
effective_population_size_2 <- 1000;
effective_population_size_admixed <- 1000;
time_split <- 10000;
effective_population_size_ancestral <- 2000;
number_of_blocks <- 10;
percentage_from_population_1 <- 0.5;

time.admixture <- rep(NA,100);
sd.percentage.admixture.estimated <- rep(NA,100);
for(rep in 1:100)
{
  time_admixture <- round(runif(1,1,5000));
  time.admixture[rep] <- time_admixture;
  sim <- model.admixed.pop(sample_1, sample_2, sample_admix,effective_population_size_1, effective_population_size_2, effective_population_size_admixed, time_admixture, percentage_from_population_1, time_split, effective_population_size_ancestral, number_of_blocks)
  # first 5 columns are information of the SNPs
  haplotypes <- sim$haplotype_matrix;
  # Genotypes are computed by collapsing two contiguous haplotypes
  genotypes <- haplotypes[seq(from=1,to=nrow(haplotypes),by=2),] + haplotypes[seq(from=2,to=nrow(haplotypes),by=2),];
  # Provide labels
  pop <- rownames(haplotypes)[seq(from=1,to=nrow(haplotypes),by=2)];
  rownames(genotypes) <- pop;
  # Compute the standard deviation of each marker to decide if it is worth being included
  sd.snps <- apply(genotypes,2,sd);
  consider.snps <- which(sd.snps > 0.1);
  genotypes.clean <- genotypes[,consider.snps];
  # scale each marker by substracting the mean and dividing by the standard deviation
  genotypes.clean.scaled <- scale(genotypes.clean);
  # compute the euclidean distance between each two individuals
  dmatrix <- dist(genotypes.clean.scaled);
  # conduct the classical multidimensional scaling
  mds.result <- cmdscale(dmatrix, k = (nrow(as.matrix(dmatrix))-1),add=T,eig=T);
  points.mds <- mds.result$points[,1:2];
  rownames(points.mds) <- pop;
  sd.percentage.admixture.estimated[rep] <- sd(extract.admixture(points.mds))  
}

plot(time.admixture, sd.percentage.admixture.estimated);

#########################################################################
# 6. Grid of 2D stepping stone pops
# Goal: to understand from a coalescence point of view what admixture is in spatial models
# Exercise:
# Create function called model.admixed.pop with a fastSimcoal2 model that 
# models a two dimensional grid of populations, connected to the geographically closest
# elements.

model.2D.grid.pop <- function(sample_sizes, nes, migration_rate, number_of_blocks)
{
  if(length(sample_sizes)!=length(nes))
  {
    print("The lenght of the vector of sample_size must be the same as Ne!");
    print(c(length(sample_sizes),length(nes)));
    return(NA);
  }
  # size of grid
  n <- length(sample_sizes)^0.5
  
  if(round(n)-n!=0)
  {
    print("Length of populations must be pow of an integer");
    return(NA);
  }
  total_pops <- length(sample_sizes);
  grow_rates <- rep(0, total_pops);
  migration_grid <- matrix(nrow=total_pops, ncol = total_pops, 0);
  # Isolation by distance
  for(po in 1:total_pops)
  {
    pop <- po -1;
    row_pop1 <- floor(pop/n);
    col_pop1 <- pop%%n;
    for(po2 in 1:total_pops)
    {
      pop2 <- po2 - 1;
      row_pop2 <- floor(pop2/n);
      col_pop2 <- pop2%%n;      
      d <- (row_pop1-row_pop2)^2 + (col_pop1-col_pop2)^2;
      if(d!=0 && d <= 2)
      {
        # receives in forward from r, c isolation migrants
        migration_grid[po,po2] <- 1;                
      }      
    }
  }
  # scale by column
  for(c in 1:total_pops)
  {
    migration_grid[,c] <- migration_rate*migration_grid[,c];
  }
  
  mmatrix <- apply(migration_grid, 1, function(row) paste(row, collapse = " "))
  
  lines <- c(
    "//Number of population samples (demes)",
    total_pops,
    "//Population effective sizes (number of genes)",
    nes,
    "//Sample sizes",
    sample_sizes,
    "//Growth rates\t: negative growth implies population expansion",
    grow_rates,
    "//Number of migration matrices : 0 implies no migration between demes",
    "1",
    "//migration matrix",
    mmatrix,
    "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
    "0 historical event",
    "//Number of independent loci [chromosome]",
    paste(number_of_blocks,"1",sep=" ")
  )
  
  # Repeating the block 22 times
  block <- c(
    "//Per chromosome: Number of linkage blocks",
    "1",
    "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
    "DNA 1000000 1.0E-8 1.855284327902964E-8 0.0"
  )
  
  # Append number_of_blocks blocks
  for (i in 1:number_of_blocks) {
    lines <- c(lines, block)
  }
  
  # Write to file
  writeLines(lines, "DemographicModelSplitR.par")
  
  exe <- ".\\fsc28.exe"
  args <- c("-i", ".\\DemographicModelSplitR.par", "-x", "-s0", "-d", "-n", "1", "-q", "-G")
  # Execute
  system2(exe, args = args)
  
  data.t <- read.table(file=paste(folder.fastSimcoal2,"\\DemographicModelSplitR\\DemographicModelSplitR_1_1.gen", sep=""), header = T);
  
  # First four columns are snp info
  # haplotype matrix. Rows are haplotypes, columns are positions
  
  H <- t(as.matrix(data.t[,5:ncol(data.t)]));
  # First 100 haplotypes is pop 1, next pop 2
  lab <- c();
  for(s in 1:total_pops)
  {
    lab <- c(lab,c(rep(s,sample_sizes[s])));
  }
  rownames(H) <- lab
  return(list(position = data.t[,1:4],haplotype_matrix = H));
}

#########################################################################
# 7. MDS output from a 2D stepping stone model
# Goal: to understand how an MDS works in a 2D stepping stone
# Exercise:
# Generate 16 populations (a grid of 4 by 4 populations), each with Ne = 500 chromosomes
# sample 20 chromosomes by population. Run MDS and plot the result with ggplot.

sample_sizes <- rep(20,16);
Ne <- rep(500, 16);
migration_rate <- 10^-5;
number_of_blocks <- 10;
sim <- model.2D.grid.pop(sample_sizes, Ne, migration_rate, number_of_blocks)
# first 5 columns are information of the SNPs
haplotypes <- sim$haplotype_matrix;
# Genotypes are computed by collapsing two contiguous haplotypes
genotypes <- haplotypes[seq(from=1,to=nrow(haplotypes),by=2),] + haplotypes[seq(from=2,to=nrow(haplotypes),by=2),];
# Provide labels
pop <- rownames(haplotypes)[seq(from=1,to=nrow(haplotypes),by=2)];
rownames(genotypes) <- pop;
# Compute the standard deviation of each marker to decide if it is worth being included
sd.snps <- apply(genotypes,2,sd);
consider.snps <- which(sd.snps > 0.1);
genotypes.clean <- genotypes[,consider.snps];
# scale each marker by substracting the mean and dividing by the standard deviation
genotypes.clean.scaled <- scale(genotypes.clean);
# compute the euclidean distance between each two individuals
dmatrix <- dist(genotypes.clean.scaled);
# conduct the classical multidimensional scaling
mds.result <- cmdscale(dmatrix, k = (nrow(as.matrix(dmatrix))-1),add=T,eig=T);
points.mds <- mds.result$points[,1:2];
rownames(points.mds) <- pop;

ma <- data.frame(model = pop, x = points.mds[,1], y = points.mds[,2]);

ggplot(ma, aes(x=x, y=y, color=model)) +
  geom_point() +
  labs(
    x = paste("DIM 1 (",round(100*mds.result$eig[1]/sum(mds.result$eig),2) ," %)"),
    y = paste("DIM 2 (",round(100*mds.result$eig[2]/sum(mds.result$eig),2) ," %)")
  )

plot(mds.result$eig); 
# plot the result (nicer with ggplot2, of course)
plot(mds.result$points[,1:2], xlab = paste("DIM 1 (",round(100*mds.result$eig[1]/sum(mds.result$eig),2) ," %)"), ylab = paste("DIM 2 (",round(100*mds.result$eig[2]/sum(mds.result$eig),2) ," %)"))

#########################################################################
# 7. Biased samples
# Goal: to understand how the MDS can be biased by different sample sizes
# Exercise:
# Take the distance matrix. Remove the first 6 rows and columns. Run the MDS with the new matrix.

dmatrix.lite <- as.matrix(dmatrix)[-c(1,6),-c(1,6)];
mds.result <- cmdscale(dmatrix.lite, k = (nrow(dmatrix.lite)-1),add=T,eig=T);
points.mds <- mds.result$points[,1:2];
rownames(points.mds) <- pop[-c(1,6)];
plot(mds.result$eig); 

ma <- data.frame(model = pop[-c(1,6)], x = points.mds[,1], y = points.mds[,2]);
ggplot(ma, aes(x=x, y=-y, color=model)) +
  geom_point() +
  labs(
    x = paste("DIM 1 (",round(100*mds.result$eig[1]/sum(mds.result$eig),2) ," %)"),
    y = paste("DIM 2 (",round(100*mds.result$eig[2]/sum(mds.result$eig),2) ," %)")
  )

# This behaviour occurs also if we have samples that are more related than expected!

#########################################################################
# 8. Making a plink file
# Goal: to understand how to manage plink files
# Exercise:
# Create a dataset with one of the models and store it in plink

effective_population_size_1 <- 1000;
effective_population_size_2 <- 1000;
effective_population_size_admixed <- 1000;
time_split <- 10000;
effective_population_size_ancestral <- 2000;
number_of_blocks <- 10;
percentage_from_population_1 <- 0.5;
time.admixture <- rep(NA,100);
sd.percentage.admixture.estimated <- rep(NA,100);
time_admixture <- 20;
sim <- model.admixed.pop(sample_1, sample_2, sample_admix,effective_population_size_1, effective_population_size_2, effective_population_size_admixed, time_admixture, percentage_from_population_1, time_split, effective_population_size_ancestral, number_of_blocks)
# first 5 columns are information of the SNPs
haplotypes <- sim$haplotype_matrix;
# Genotypes are computed by collapsing two contiguous haplotypes
genotypes <- haplotypes[seq(from=1,to=nrow(haplotypes),by=2),] + haplotypes[seq(from=2,to=nrow(haplotypes),by=2),];
# Provide labels
pop <- rownames(haplotypes)[seq(from=1,to=nrow(haplotypes),by=2)];
rownames(genotypes) <- pop;
# Compute the standard deviation of each marker to decide if it is worth being included
sd.snps <- apply(genotypes,2,sd);
consider.snps <- which(sd.snps > 0.1);
genotypes.clean <- genotypes[,consider.snps];

# store the genotype matrix in Plink format
bim.data <- data.frame(sim$position[consider.snps,1], paste("rs",sim$position[consider.snps,2],sep=""), rep(0,length(consider.snps)), sim$position[consider.snps,2], sim$position[consider.snps,3],sim$position[consider.snps,4]);
colnames(bim.data) <- c("chr", "id", "posg", "pos", "alt", "ref");
fam.data <- data.frame(pop,paste(pop,seq(from=1,to=length(pop),by=1),sep="_"),rep(0,length(pop)),rep(0,length(pop)),rep(0,length(pop)),rep(-9,length(pop)));
colnames(fam.data) <- c("fam", "id", "pat", "mat", "sex", "pheno");
bed.matrix <- t(genotypes.clean);
colnames(bed.matrix) <- fam.data$id;
write_plink("modelAdmixture_1_1", bed.matrix, bim = bim.data, fam = fam.data);

# Lets check that the dataset was properly written
kk <- BEDMatrix("modelAdmixture_1_1")
rm(kk);

# You can use the bed file to run ADMIXTURE in linux
# Two files are produced. Q -> The admixture proportions of each individual. P -> the frequencies in the ancestral 

#########################################################################
# 10. Estimating ancestry proportions with sNMF
# Goal: to estimate ancestry proportions using one algorithm for global ancestry estimation 
# Exercise:
# use sNMF, which has been coded in R, to estimate ancestry proportions. BEWARE: Results can be different
# from the ones from ADMIXTURE as the algorithm for ancestry optimization and the assumptions of the unsupervised model
# ARE DIFFERENT

# Step 1
# create the geno data
#genotypes.clean has n individuals as rows and m SNVs as columns
# we use the clean, because we will have to take into account 
write.geno(genotypes.clean, "genotypes.geno")
project = NULL;
project = snmf("genotypes.geno",
               K = 1:5,
               entropy = TRUE,
               repetitions = 2,
               project = "new");

# How good is each model (at each K)?
plot(project, col = "blue", pch = 19, cex = 1.2)
# Extract Q matrix (ancestry coefficients)
Q <- as.matrix(Q(project, K = 2, run = best))
# Create an index that orders individuals by population label
ordered_indices <- order(pop)
# Reorder Q matrix and population labels
Q_ordered <- Q[ordered_indices, ]
pop_ordered <- pop[ordered_indices]
# Define colors
my.colors <- c("tomato", "gold")

# Plot the barplot manually using barplot(), not barchart()
barplot(t(Q_ordered),
        col = my.colors,
        border = NA,
        space = 0,
        xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = "Ancestry matrix")

# Add population labels on the x-axis
axis(1, at = seq(0.5, nrow(Q_ordered) - 0.5, by = 1),
     labels = pop_ordered,
     las = 2,
     cex.axis = 0.4)

read.table(file=paste(project@runs[best][[1]]@directory,project@runs[best][[1]]@G.output.file, sep=""),header=F)

###########################################################################
# 11. Explore and Convert Genotype Data
# Goal: Load a genotype matrix and inspect its structure.
# Exercise:
#  Load a PLINK .bed/.bim/.fam) using SNPRelate.
# Convert it into a GDS (Genomic Data Structure) format.
# Extract a genotype matrix for further analysis.

library(SNPRelate)

# Convert VCF to GDS (use your own file path)
snpgdsVCF2GDS("your_data.vcf", "your_data.gds", method = "biallelic.only")
genofile <- snpgdsOpen("your_data.gds")

# Inspect sample and SNP data
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))

#########################################################################
# 12. Principal Component Analysis (PCA)
# Goal: Visualize substructure based on genetic variation.
# Exercise:
# Perform PCA on the genotype data.
# Plot the first two principal components.

pca <- snpgdsPCA(genofile, autosome.only = TRUE)
pc.percent <- pca$varprop * 100

# Tidy data for plotting
pca_df <- data.frame(
  sample.id = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2]
)

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA of Genotype Data")

#########################################################################
# 13. Cluster Individuals Based on Genetic Similarity
# Goal: Use clustering methods to define subgroups.
# Exercise:
# Run a clustering algorithm (e.g., k-means or DAPC from adegenet).
# Visualize population clusters.

library(adegenet)

# Extract genotypes into a genlight object
genlight_obj <- snpgds2genlight(genofile)

# Run DAPC
grp <- find.clusters(genlight_obj, max.n.clust = 10)
dapc_res <- dapc(genlight_obj, grp$grp)

scatter(dapc_res, scree.da = TRUE)

#########################################################################
# 14. Cluster Individuals after feature extraction
# Goal: Use clustering methods to define subgroups.
# Exercise:
# Run mclust and dbscan to visualize population clusters in the PCA.

library(mclust);

clus <- Mclust(smac$conf[,1:5], G= 1:10)
plot(clus)

library("dbscan");
cl <- hdbscan(smac$conf[,1:2], minPts = 5)
plot(smac$conf[,1:2], col=cl$cluster+1, pch=20)

#########################################################################
# 15. Pairwise Genetic Distances
# Goal: Compute and visualize genetic distances between individuals.
# Exercise:
# Calculate IBS (identity-by-state) distance matrix.
# Visualize with a heatmap.

ibs <- snpgdsIBS(genofile, num.thread=2)
dist_matrix <- 1 - ibs$ibs
heatmap(as.matrix(dist_matrix), Rowv=NA, Colv=NA, main="IBS Distance")

#########################################################################
# 16. Working with other types of distances
# Goal: To see how to project in a low dimensional space genetic distance matrices that are
# non-Euclidean.
# Exercise:
# Use SMACOF ordinal to reduce the dimensionality of IBD matrix.

rm(list= ls());
data.fibla <- read.table(file = "Location_lat_lon_population395.txt", header = T);

ibd.matrix <- matrix(nrow = nrow(data.fibla), ncol = nrow(data.fibla),0);
rownames(ibd.matrix) <- data.fibla$DNA_code;
colnames(ibd.matrix) <- data.fibla$DNA_code;

for(chr in 1:22)
{
  data.hbd <- read.table(gzfile(paste("out.chr", chr, ".hbd.gz", sep="")), header =F);
  for(e in 1:nrow(data.hbd))
  {
    name <- strsplit(data.hbd[e,1],"_")[[1]][2]
    r <- which(data.fibla$DNA_code==name);
    ibd.matrix[r,r] <- ibd.matrix[r,r] + data.hbd[e,8];
  }
  data.ibd <- read.table(gzfile(paste("out.chr", chr, ".ibd.gz", sep="")), header =F);
  for(e in 1:nrow(data.ibd))
  {
    name_one <- strsplit(data.ibd[e,1],"_")[[1]][2];
    name_two <- strsplit(data.ibd[e,3],"_")[[1]][2];
    r_one <- which(data.fibla$DNA_code==name_one);
    r_two <- which(data.fibla$DNA_code==name_two);
    ibd.matrix[r_one,r_two] <- ibd.matrix[r_one,r_two] + data.ibd[e,8];
    ibd.matrix[r_two,r_one] <- ibd.matrix[r_one,r_two];
  }  
}

dis.ibd.matri<- (max(ibd.matrix)- ibd.matrix)/(max(ibd.matrix)-min(ibd.matrix));

diag(dis.ibd.matri) <- 0;

# Classical multidimensional scaling
cmd <- (cmdscale(dis.ibd.matri))

# non-metric multidimensional scaling with SMACOF
smac <- mds(dis.ibd.matri, 14, type = "ordinal");

# compare both matrices
res.prot <- protest(smac$conf[,1:2], cbind(data.fibla$Lon_X, data.fibla$Lat_Y), permutations = 9999);

res.prot.cmd <- protest(cmd, cbind(data.fibla$Lon_X, data.fibla$Lat_Y), permutations = 9999);

#########################################################################
# 17. Local ancestry
# Goal: To see that the genome can be modeled as a mosaic of ancestries
# Exercise:
# Build a file for RFMix from simulated data (which we know it is already PHASED).
# This is the optimal scenario as the data is for sure phased
# In real life. 1) Phase your data (ShapeIt, Beagle), generate the files
# vcf, reference ancestry, genetic map
# run rfmix with phasing errors option

rm(list = ls());

#############################
folder.fastSimcoal2 <- "C:\\Users\\u9424\\OneDrive - Universitat Pompeu Fabra\\Grants\\2025\\EMBO Naples June\\Second day\\";

setwd(folder.fastSimcoal2);

model.admixed.pop <- function(sample_1, sample_2, sample_admix, effective_population_size_1, effective_population_size_2, effective_population_size_admixed, time_admixture, percentage_from_population_1, time_split, effective_population_size_ancestral, number_of_blocks, length_block)
{
  lines <- c(
    "//Number of population samples (demes)",
    "3",
    "//Population effective sizes (number of genes)",
    effective_population_size_1,
    effective_population_size_2,
    effective_population_size_admixed,
    "//Sample sizes",
    sample_1,
    sample_admix,
    sample_2,
    "//Growth rates\t: negative growth implies population expansion",
    "0",
    "0",
    "0",
    "//Number of migration matrices : 0 implies no migration between demes",
    "0",
    "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
    "3 historical event",
    paste(time_admixture,"1 0",percentage_from_population_1,"1 0 0",sep = " "),
    paste(time_admixture,"1 2 1 1 0 0",sep = " "),
    paste(time_split,"0 2 1",effective_population_size_ancestral/effective_population_size_2, "0 0",sep = " "),
    "//Number of independent loci [chromosome]",
    paste(number_of_blocks,"1",sep=" ")
  )
  
  # Repeating the block 22 times
  block <- c(
    "//Per chromosome: Number of linkage blocks",
    "1",
    "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
    paste("DNA", length_block, "1.0E-8 1.855284327902964E-8 0.0", sep = " ")
  )
  
  # Append number_of_blocks blocks
  for (i in 1:number_of_blocks) {
    lines <- c(lines, block)
  }
  
  # Write to file
  writeLines(lines, "DemographicModelSplitR.par")
  
  exe <- ".\\fsc28.exe"
  args <- c("-i", ".\\DemographicModelSplitR.par", "-x", "-s0", "-d", "-n", "1", "-q", "-G")
  # Execute
  system2(exe, args = args)
  
  data.t <- read.table(file=paste(folder.fastSimcoal2,"\\DemographicModelSplitR\\DemographicModelSplitR_1_1.gen", sep=""), header = T);
  
  # First four columns are snp info
  # haplotype matrix. Rows are haplotypes, columns are positions
  
  H <- t(as.matrix(data.t[,5:ncol(data.t)]));
  # First 100 haplotypes is pop 1, next pop 2
  rownames(H) <- c(rep("A",sample_1),rep("Admix",sample_admix),rep("B",sample_2));
  return(list(position = data.t[,1:4],haplotype_matrix = H));
}

#############################

setwd("C:\\Users\\u9424\\OneDrive - Universitat Pompeu Fabra\\Grants\\2025\\EMBO Naples June\\Second day\\Morning\\Practical\\")
effective_population_size_1 <- 1000;
effective_population_size_2 <- 1000;
effective_population_size_admixed <- 1000;
time_split <- 10000;
effective_population_size_ancestral <- 2000;
number_of_blocks <- 10;
percentage_from_population_1 <- 0.5;
time_admixture <- 20;
sample_1 <- 100;
sample_2 <- 100;
sample_admix <- 100;
sim <- model.admixed.pop(sample_1, sample_2, sample_admix,effective_population_size_1, effective_population_size_2, effective_population_size_admixed, time_admixture, percentage_from_population_1, time_split, effective_population_size_ancestral, 1,20000000)

he <- "##fileformat=VCFv4.2"

vcf_header <- function(samples)
{
  vcf <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sep = "\t");
  samples.m <- paste(samples,collapse = "\t");
  vcf_header= paste(vcf, samples.m, sep = "\t");  
  return(vcf_header);
}

vcf_position <- function(position)
{
  p <- as.character(position);
  m <- paste(p[1], p[2],".",p[3],p[4],".",  "PASS", ".", "GT", sep = "\t", collapse = "");
  return(gsub(" ", "", m))  
}

positions <- apply(sim$position,1,vcf_position);

vcf_haplotype <- function(h)
{
  h1 <- h[seq(from=1,to=length(h),by=2)];
  h2 <- h[seq(from=2,to=length(h),by=2)];
  hh <- paste(h1,h2,sep="|");
  return(paste(hh,collapse="\t"));
}

# region map
genetic_map <- function(position, recombination_rate)
{
  chr = position$Chrom;
  pos = position$Pos;
  genetic_dist <- position$Pos*recombination_rate;
  return(data.frame(chr,pos, genetic_dist));
}

sam <- rownames(sim$haplotype_matrix);
sam <- sam[seq(from=1,to=length(sam),by=2)];
samples <- paste(sam,"_",as.character(seq(from = 1, to= length(sam), by = 1)),sep ="");

samples.admix <- samples[sam=="Admix"];
admix <- sim$haplotype_matrix[rownames(sim$haplotype_matrix)=="Admix",];
v <- vcf_header(samples.admix);

whole.matrix.admix <- c(he,v,paste(positions,apply(admix,2,vcf_haplotype), sep = "\t"));
writeLines(whole.matrix.admix, gzfile("admix.vcf.gz"))

samples.reference <- samples[sam!="Admix"];
reference <- sim$haplotype_matrix[rownames(sim$haplotype_matrix)!="Admix",];
v <- vcf_header(samples.reference);
whole.matrix.reference <- c(he,v,paste(positions,apply(reference,2,vcf_haplotype), sep = "\t"));
writeLines(whole.matrix.reference, gzfile("reference.vcf.gz"))

# Sample map. rfmix is REALLY tricky with spaces. USE TABs!

ancestry <- sam[sam!="Admix"];
sample_map <- paste(samples[sam!="Admix"],ancestry[seq(from=1,to=length(ancestry),by=1)], sep ="\t");
writeLines(c("",sample_map), "sample_map.txt")
write.table(genetic_map(sim$position, 10^-8), "genetic_map_chr.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

exe <- "rfmix"
args <- c("-f", "admix.vcf.gz", "-r", "reference.vcf.gz", "-m", "sample_map.txt", "-g", "genetic_map_chr.txt","-o", "admixed_output", "--chromosome=1","--n-threads=11")
system2(exe, args = args)

# run rfmix
# ./rfmix -f admix.vcf.gz -r reference.vcf.gz -m sample_map.txt -g genetic_map_chr.txt -o admixed_output --chromosome=1 --n-threads=11

#########################################################################
# 18. Local ancestry
# Goal: To see that the genome can be modeled as a mosaic of ancestries
# Exercise:
# Build a file for RFMix from simulated data (which we know it is already PHASED).

#Output
#*.fb.tsv – Forward-backward posterior probabilities (key file)
#*.msp.tsv – Most likely ancestry tracts (compact summary)
#*.alleles – Copy of the input alleles
#*.classes – Copy of the reference panel population labels
#*.log – Log of the run

# Read RFMix msp.tsv
rfmix_data <- read.table("admixed_output.msp.tsv", header = TRUE, sep = "\t", skip = 1, comment.char = "@")

# How does it looks like? Extract for one of the individuals the ancestry A

