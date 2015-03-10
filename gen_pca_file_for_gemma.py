# 2, 4, 6 pc files

num_pcs = [0,2,4,6]

f_all = open('prcomp_all_snps_x.txt', 'r')
pcs_all = [[] for i in range(max(num_pcs))]
f_all.readline()
line = f_all.readline()
while line != '':
    line = line.strip().split()
    for i in range(max(num_pcs)):
        pcs_all[i].append(line[i])
    line = f_all.readline()
f_all.close()

for n in num_pcs:
    f_cov_all = open('pc_' + str(n) + '_all.txt', 'w')
    for row in range(len(pcs_all[0])):
        f_cov_all.write('1')
        for j in range(n):
            f_cov_all.write(' ' + pcs_all[j][row])
        f_cov_all.write('\n')
    f_cov_all.close()

f_common = open('prcomp_common_snps_x.txt', 'r')
pcs_common = [[] for i in range(max(num_pcs))]
f_common.readline()
line = f_common.readline()
while line != '':
    line = line.strip().split()
    for i in range(max(num_pcs)):
        pcs_common[i].append(line[i])
    line = f_common.readline()
f_common.close()

for n in num_pcs:
    f_cov_common = open('pc_' + str(n) + '_common.txt', 'w')
    for row in range(len(pcs_common[0])):
        f_cov_common.write('1')
        for j in range(n):
            f_cov_common.write(' ' + pcs_common[j][row])
        f_cov_common.write('\n')
    f_cov_common.close()

