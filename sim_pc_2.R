library(data.table)
args = commandArgs(trailingOnly = TRUE)
id = args[1]
new_prefix = args[2]
X = as.matrix(read.table('pc_2_all.txt', header=F))
Y = as.matrix(read.table(paste('output/', new_prefix, '.phen.txt', sep=''), header=F))
genotypes = fread('subset_302_homozygotes_ped_for_R.txt', header = F)[,1:302,with=F]
# add 1 because indexed by 0
marker_inds = sort(read.table('subset_302_homozygotes_tag_inds.txt', header = F)) + 1
n_sites = length(marker_inds)
pvals_coef = rep(-1, n_sites)
## intercept is in file but automatically included in lm so we'll
## replace first column with genotypes
i = 1
for (m in marker_inds) {
    X[,1] = t(genotypes[m])
    m = summary(lm(Y ~ X))
    #pvals_model[i] = pf(m$fstatistic[1], m$fstatistic[2], m$fstatistic[3], lower.tail = FALSE)
    pvals_coef[i] = m$coefficients[2,4]
    i = i + 1
}
write.table(pvals_coef, file=paste('output/', new_prefix,
                            '_no_K_pc_2.assoc.txt', sep = ''),
            row.names = F)
