f = open('subset_302_homozygotes_pruned.8.Q', 'r')
f_cov = open('admixture_covariates_8_pops.txt', 'w')
line = f.readline()
while line != '':
    line = line.split()
    f_cov.write('1')
    for p in line[:-1]:
        f_cov.write(' ' + p)
    f_cov.write('\n')
    line = f.readline()
f.close()
f_cov.close()
