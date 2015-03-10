## arguments are effect size and run id
args = commandArgs(trailingOnly = TRUE)
data <- read.table(file="/net/gs/vol1/home/aclark4/population_structure/code/snps_processed.txt",header=T,sep="\t")

library(emma)
nof_individuals <- ncol(data) - 2
nof_snps <- nrow(data)
num_chr = max(data$chr)

## get values needed for EMMA calculations
m <- as.matrix(data[,3:ncol(data)])

## define variables
causal_sites <- rep(0,nof_simulations)
all_pvalues <- matrix(data=NA, nrow=nof_snps, ncol=nof_simulations)
all_phens <- matrix(data=NA, nrow=nof_individuals, ncol=nof_simulations)

## pick a random snp as the causal site
random_snp <- sample(1:nrow(data),1)
causal_sites[i] <- random_snp
print(random_snp)
print(data[random_snp,])
	
## get the indices of the two alleles
data1 <- m[random_snp,]
non_NA <- which(!is.na(data1))
data1 <- data1[non_NA]
alleles <- unique(data1)
a1.index <- which(data1==alleles[1])
a2.index <- which(data1==alleles[2])

simulated_phenotypes <- rep(NA,length(data1))
phenotypes_all_strains <- rep(NA,nof_individuals)
	
## fixed effect size (difference in mean between two groups, standard deviation = 1)
fixed_effect <- as.numeric(args[1])
	
## simulate data for the two alleles
simulated.allele1 <- rnorm(length(a1.index),fixed_effect,1)
simulated.allele2 <- rnorm(length(a2.index),0,1)
simulated_phenotypes[a1.index] <- simulated.allele1
simulated_phenotypes[a2.index] <- simulated.allele2

phenotypes_all_strains[non_NA] = simulated_phenotypes
all_phens[,i] = phenotypes_all_strains

## use EMMA to test for genome-wide association,
## with kinship matrix for each chromosome
K_all <- emma.kinship(m[,non_NA],method="additive",use="all")
m_s = m[,non_NA]
ys <- as.matrix(simulated_phenotypes)
rs <- emma.REML.t(ys,m_s,K_all)
p <- rs$ps
all_pvalues[,i] <- p
        
save(causal_sites,all_pvalues,all_phens,file=paste("/net/gs/vol1/home/aclark4/population_structure/code/output/outfile_effect",args[1], "_", args[2],".Rdata",sep=""))
