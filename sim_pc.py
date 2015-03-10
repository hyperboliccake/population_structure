import sys
import os
import random
import copy
import scipy.stats
import shutil
from helpers import *

def main():

    orig_prefix = 'subset_302_homozygotes'
    new_prefix = orig_prefix + '_' + sys.argv[1]

    shuffle = False
    if len(sys.argv) > 3 and sys.argv[3] == 'shuffle':
        new_prefix += '_shuffled'
        shuffle = True

    print('choosing causal snp')
    # choose random snp from all snps that meet filter criteria (this is
    # more than just the tag snps)
    f_filtered = open(orig_prefix + '_filtered_inds.txt', 'r')
    snp_inds = [int(i) for i in f_filtered.readline().split()]
    f_filtered.close()
    random_snp_ind = random.choice(snp_inds) # index in original ped/map
    f_causal = open('output/' + new_prefix + '.causal.txt', 'w')
    f_causal.write(str(random_snp_ind) + '\n')
    f_causal.close()

    print('reading in genotypes and assigning phenotypes')
    # assign phenotypes based on chosen snp (alleles aren't converted to
    # 0s and 1s, and rows are strains)
    ped_lines = [line.split() for line in open(orig_prefix + '.ped', 'r').readlines()]
    num_strains = int(len(ped_lines))
    num_snps = int((len(ped_lines) - 6) / 2)
    # these genotypes are from the unfiltered ped file, so the indices in
    # snp_inds will correctly correspond to these rows
    alleles = []
    for strain in range(num_strains):
        alleles.append(ped_lines[strain][6::2][random_snp_ind])
    alleles_set = list(set(alleles))
    major = alleles_set[0]
    minor = alleles_set[1]
    if alleles.count(minor) > alleles.count(major):
        major = alleles_set[1]
        minor = alleles_set[0]
    maf = float(alleles.count(minor)) / len(alleles)
    phens = []

    h2 = float(sys.argv[2])
    fixed_effect = pow(float(h2 * (len(alleles)-1)) / ((h2-1) * len(alleles) * (maf-1) * maf), .5)

    for allele in alleles:
        mu = 0
        if allele == minor:
            # fixed effect size - should this depend on minor
            # allele frequency to keep power constant?
            mu = fixed_effect
        phens.append(random.gauss(mu, 1))

    # shuffle?
    if shuffle:
        random.shuffle(phens)

    # store the phenotypes so we can look at them later
    f_phen = open('output/' + new_prefix + '.phen.txt', 'w')
    for p in phens:
        f_phen.write(str(p) + '\n')
    f_phen.close()
    

    # create fam files
    phenotypes_to_fam(orig_prefix, new_prefix, phens)

    # association tests
    print('running association tests') 

    #####
    # whole-genome kinship matrix plus varying number of principle components
    #####
    print('- whole genome K')

    os.symlink(orig_prefix + '_tag.bed', new_prefix + '.bed')
    os.symlink(orig_prefix + '_tag.bim', new_prefix + '.bim')

    num_pc = [0, 2, 4, 6]
    pc_types = ['all', 'common']

    for pc_type in pc_types:
        for n in num_pc:
            kinship_prefix_d = 'output/' + orig_prefix
            kinship_prefix_u = 'output/' + orig_prefix
            os.system('/net/gs/vol1/home/aclark4/software/gemma' + \
                          ' -bfile ' + new_prefix + \
                          ' -d ' + kinship_prefix_d + '.eigenD.txt' + \
                          ' -u ' + kinship_prefix_u + '.eigenU.txt' + \
                          ' -lmm 4' + \
                          ' -notsnp' + \
                          ' -o ' + new_prefix + \
                          ' -c pc_' + str(n) + '_' + pc_type + '.txt')
            os.system('mv output/' + new_prefix + '.assoc.txt output/' + new_prefix + '_whole_K_pc_' + str(n) + '_' + pc_type + '.assoc.txt')

    os.remove(new_prefix + '.bed')
    os.remove(new_prefix + '.bim')

    # then do all the same for no kinship matrix
    for pc_type in pc_types:
        os.system('Rscript sim_pc.R ' + sys.argv[1] + ' ' + pc_type)

 
    #####
    # t-test (just to confirm this matches kinship 0/no pc case)
    #####
    print('- t-test')

    tag_inds = [int(x) for x in open(orig_prefix + '_tag_inds.txt', 'r').readline().split()]
    tag_inds.sort()
    num_markers = len(tag_inds)

    f_assoc_t = open('output/' + new_prefix + '_t.assoc.txt', 'w')
    for i in tag_inds:
        alleles = []
        for strain in range(num_strains):
            alleles.append(ped_lines[strain][6 + 2 * i])
        a1 = list(set(alleles))[0]
        phens_1 = []
        phens_2 = []
        for phen_ind in range(len(phens)):
            if alleles[phen_ind] == a1:
                phens_1.append(phens[phen_ind])
            else:
                phens_2.append(phens[phen_ind])

        p = scipy.stats.ttest_ind(phens_1, phens_2)[1]
        f_assoc_t.write(str(p) + '\n')
    f_assoc_t.close()
     
    os.remove(new_prefix + '.fam')

main()
