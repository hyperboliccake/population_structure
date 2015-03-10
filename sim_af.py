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
    # basic whole-genome kinship matrix
    #####
    print('- whole genome K')

    os.symlink(orig_prefix + '_tag.bed', new_prefix + '.bed')
    os.symlink(orig_prefix + '_tag.bim', new_prefix + '.bim')
    gemma(new_prefix, 'output/' + orig_prefix, 'output/' + orig_prefix, '_whole_K')
    os.remove(new_prefix + '.bed')
    os.remove(new_prefix + '.bim')
    os.remove('output/' + new_prefix + '_whole_K.log.txt')

    #####
    # kinship matrix by allele frequency bins
    #####
    print('- K by af bins')

    f_assoc_af = open('output/' + new_prefix + '_af.assoc.txt', 'w')

    bed_bim_files = sorted(os.listdir('bed_af'))
    for i in range(len(bed_bim_files)/2):
        orig_prefix_af = bed_bim_files[2 * i][:-4]
        new_prefix_af = orig_prefix_af + '_' + sys.argv[1]

        if shuffle:
            new_prefix_af += '_shuffled'

        os.symlink(new_prefix + '.fam', new_prefix_af + '.fam')
        os.symlink('bed_af/' + orig_prefix_af + '.bed', new_prefix_af + '.bed')
        os.symlink('bed_af/' + orig_prefix_af + '.bim', new_prefix_af + '.bim')
        
        gemma(new_prefix_af, 'output/' + orig_prefix + '_af_eigenD/' + orig_prefix_af, 'output/' + orig_prefix + '_af_eigenU/' + orig_prefix_af)

        f_out_af = open('output/' + new_prefix_af + '.assoc.txt', 'r')
        f_out_af.readline()
        line = f_out_af.readline()
        while line != '':
            f_assoc_af.write(line)
            line = f_out_af.readline()
        f_out_af.close()

        os.remove('output/' + new_prefix_af + '.assoc.txt')
        os.remove('output/' + new_prefix_af + '.log.txt')
        os.remove(new_prefix_af + '.fam')
        os.remove(new_prefix_af + '.bim')
        os.remove(new_prefix_af + '.bed')

    f_assoc_af.close()
    
    os.remove(new_prefix + '.fam')

main()
