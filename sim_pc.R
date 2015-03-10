library(data.table)
args = commandArgs(trailingOnly = TRUE)
id = args[1]
pc_type = args[2]
pc = as.matrix(read.table(paste('pc_6_', pc_type, '.txt', sep=''), header=F))
Y = as.matrix(read.table(paste('output/subset_302_homozygotes_',
    id, '.phen.txt', sep=''), header=F))
genotypes = fread('subset_302_homozygotes_ped_for_R.txt', header = F)[,1:302,with=F]
# add 1 because indexed by 0
marker_inds = sort(read.table('subset_302_homozygotes_tag_inds.txt', header = F)) + 1
n_sites = length(marker_inds)
for (n in c(0, 2, 4, 6)) {
    pvals_model = rep(-1, n_sites)
    pvals_coef = rep(-1, n_sites)
    ## intercept is in file but automatically included in lm so we'll
    ## replace first column with genotypes
    X = as.matrix(pc[,1:(n+1)])
    i = 1
    for (m in marker_inds) {
        X[,1] = t(genotypes[m])
        m = summary(lm(Y ~ X))
        pvals_model[i] = pf(m$fstatistic[1], m$fstatistic[2], m$fstatistic[3], lower.tail = FALSE)
        pvals_coef[i] = m$coefficients[2,4]
        i = i + 1
    }
    write.table(pvals_model, file=paste('output/subset_302_homozygotes_', id, '_no_K_pc_',
                           n, '_', pc_type, '_model.assoc.txt', sep = ''),
                row.names = F)

    write.table(pvals_coef, file=paste('output/subset_302_homozygotes_', id, '_no_K_pc_',
                           n, '_', pc_type, '_coef.assoc.txt', sep = ''),
                row.names = F)
}
