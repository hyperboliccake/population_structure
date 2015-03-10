import sys

freq_cutoff = .05

ped_lines = open('subset_302_homozygotes.ped', 'r').readlines()

num_strains = len(ped_lines)
num_snps = (len(ped_lines[0].split()) - 6)/2
all_alleles = [[] for n in range(num_snps)]
all_alleles_2 = [[] for n in range(num_snps)]
for strain in range(num_strains):
    alleles = ped_lines[strain].split()[6:]
    for i in range(num_snps):
        all_alleles[i].append(alleles[2 * i])
        all_alleles_2[i].append(alleles[2 * i + 1])

keep_inds = []
minor = []
for i in xrange(num_snps):
    alleles = all_alleles[i] 
    if alleles == all_alleles_2[i]:
        alleles_set = list(set(alleles))
        if len(alleles_set) == 2:
            maf = float(alleles.count(alleles_set[0])) / num_strains
            m = alleles_set[0]
            if maf > .5:
                maf = 1 - maf
                m = alleles_set[1]
            if maf >= freq_cutoff:
                keep_inds.append(i)
                minor.append(m)

num_keep = len(keep_inds)

f_out = open('pairwise_ld.txt', 'w')
for i in xrange(num_keep):
    print(i)
    mi = minor[i]
    for j in xrange(i + 1, num_keep):
        print(j)
        mj = minor[j]
        i0 = 0
        j0 = 0
        i0j0 = 0
        for k in range(num_strains):
            if all_alleles[keep_inds[i]][k] == mi:
                i0 += 1
                if all_alleles[keep_inds[j]][k] == mj:
                    j0 += 1
                    i0j0 += 1
            elif all_alleles[keep_inds[j]][k] == mj:
                j0 += 1
        i0 = float(i0)/num_strains
        j0 = float(j0)/num_strains
        i0j0 = float(i0j0)/num_strains
        D = i0j0 - i0 * j0
        try:
            r2 = pow(D, 2) / (i0 * j0 * (1 - i0) * (1 - j0))
        except:
            print i, j
            print i0, j0, i0j0
            print all_alleles[keep_inds[i]]
            print all_alleles[keep_inds[j]]
            sys.stdout.flush()
        f_out.write(str(keep_inds[i]) + ' ' + str(keep_inds[j]) + ' ')
        f_out.write(str(D) + ' ' + str(r2) + ' ')
        f_out.write(str(i0) + ' ' + str(j0) + '\n')
f_out.close()        
                    
