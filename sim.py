import sys
import os
import random
import copy
import scipy.stats
import shutil
from helpers import *

def main():

    orig_prefix = 'subset_302_homozygotes'
    new_prefix = orig_prefix + '_' + sys.argv[1] + '_' + sys.argv[2]

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

    """
    #####
    # local kinship matrix (by chromosome)
    #####
    print('- local K by chromosome')

    f_assoc_local_chr = open('output/' + new_prefix + '_local_chr.assoc.txt', 'w')

    for chrm in range(1, 17):
        chrm = str(chrm)
        orig_prefix_chrm = orig_prefix + '_chr' + chrm 
        new_prefix_chrm =  orig_prefix_chrm + '_' + sys.argv[1]
        if shuffle:
            new_prefix_chrm += '_shuffled'
        os.symlink(new_prefix + '.fam', new_prefix_chrm + '.fam')
        os.symlink('bed_chr/' + orig_prefix_chrm + '.bed', new_prefix_chrm + '.bed')
        os.symlink('bed_chr/' + orig_prefix_chrm + '.bim', new_prefix_chrm + '.bim')
        
        gemma(new_prefix_chrm, 'output/' + orig_prefix + '_chr_eigenD/' + orig_prefix_chrm, 'output/' + orig_prefix + '_chr_eigenU/' + orig_prefix_chrm)

        f_out_local = open('output/' + new_prefix_chrm + '.assoc.txt', 'r')
        f_out_local.readline()
        line = f_out_local.readline()
        while line != '':
            f_assoc_local_chr.write(line)
            line = f_out_local.readline()
        f_out_local.close()

        os.remove('output/' + new_prefix_chrm + '.assoc.txt')
        os.remove('output/' + new_prefix_chrm + '.log.txt')
        os.remove(new_prefix_chrm + '.fam')
        os.remove(new_prefix_chrm + '.bim')
        os.remove(new_prefix_chrm + '.bed')

    f_assoc_local_chr.close()

    #####
    # local kinship matrix (excluding current marker's chromosome)
    #####
    print('- local K excluding current chromosome')

    f_assoc_local_exclude_chr = open('output/' + new_prefix + '_local_exclude_chr.assoc.txt', 'w')

    for chrm in range(1, 17):
        chrm = str(chrm)
        orig_prefix_chrm = orig_prefix + '_chr' + chrm 
        new_prefix_chrm =  orig_prefix_chrm + '_' + sys.argv[1]
        if shuffle:
            new_prefix_chrm += '_shuffled'
        os.symlink(new_prefix + '.fam', new_prefix_chrm + '.fam')
        os.symlink('bed_chr/' + orig_prefix_chrm + '.bed', new_prefix_chrm + '.bed')
        os.symlink('bed_chr/' + orig_prefix_chrm + '.bim', new_prefix_chrm + '.bim')
        
        gemma(new_prefix_chrm, 'output/' + orig_prefix + '_exclude_chr_eigenD/' + orig_prefix_chrm, 'output/' + orig_prefix + '_exclude_chr_eigenU/' + orig_prefix_chrm)

        f_out_local = open('output/' + new_prefix_chrm + '.assoc.txt', 'r')
        f_out_local.readline()
        line = f_out_local.readline()
        while line != '':
            f_assoc_local_exclude_chr.write(line)
            line = f_out_local.readline()
        f_out_local.close()

        os.remove('output/' + new_prefix_chrm + '.assoc.txt')
        os.remove('output/' + new_prefix_chrm + '.log.txt')
        os.remove(new_prefix_chrm + '.fam')
        os.remove(new_prefix_chrm + '.bim')
        os.remove(new_prefix_chrm + '.bed')

    f_assoc_local_exclude_chr.close()

    #####
    # local kinship matrix (by linkage group, with windows around all sites in group)
    #####
    print('- local K by linkage group')

    rsids = [int(line.split()[1][2:]) for line in open(orig_prefix + '_tag.map', 'r').readlines()]
    rsids.sort()
    tag_bim_lines = open(orig_prefix + '_tag.bim', 'r').readlines()

    f_assoc_local_ld = open('output/' + new_prefix + '_local_ld.assoc.txt', 'w')

    for i in range(len(rsids)):        
        rsid = 'rs' + str(rsids[i])
        orig_prefix_tag = orig_prefix + '_' + rsid
        new_prefix_tag =  orig_prefix_tag + '_ld_' + sys.argv[1]
        if shuffle:
            new_prefix_tag += '_shuffled'
        os.symlink(new_prefix + '.fam', new_prefix_tag + '.fam')
        os.symlink('bed/' + orig_prefix_tag + '.bed', new_prefix_tag + '.bed')
        f_bim = open(new_prefix_tag + '.bim', 'w')
        f_bim.write(tag_bim_lines[i])
        f_bim.close()
        
        gemma(new_prefix_tag, 'output/' + orig_prefix + '_ld_1000_eigenD/' + orig_prefix_tag, '/net/akey/vol2/aclark4/eigen/' + orig_prefix + '_ld_1000_eigenU/' + orig_prefix_tag)

        f_out_local = open('output/' + new_prefix_tag + '.assoc.txt', 'r')
        f_out_local.readline()
        line = f_out_local.readline()
        while line != '':
            f_assoc_local_ld.write(line)
            line = f_out_local.readline()
        f_out_local.close()

        os.remove('output/' + new_prefix_tag + '.assoc.txt')
        os.remove('output/' + new_prefix_tag + '.log.txt')
        os.remove(new_prefix_tag + '.fam')
        os.remove(new_prefix_tag + '.bim')
        os.remove(new_prefix_tag + '.bed')
        
    f_assoc_local_ld.close()

    #####
    # local kinship matrix by 1000bp window
    #####
    print('- local K by 1000bp window')

    f_assoc_local_window = open('output/' + new_prefix + '_local_window.assoc.txt', 'w')

    for i in range(len(rsids)):        
        rsid = rsids[i]
        orig_prefix_window = orig_prefix + '_' + rsid
        new_prefix_window =  orig_prefix_window + '_window_' + sys.argv[1]
        if shuffle:
            new_prefix_window += '_shuffled'
        os.symlink(new_prefix + '.fam', new_prefix_window + '.fam')
        os.symlink('bed/' + orig_prefix_window + '.bed', new_prefix_window + '.bed')
        f_bim = open(new_prefix_window + '.bim', 'w')
        f_bim.write(tag_bim_lines[i])
        f_bim.close()
        
        gemma(new_prefix_window, 'output/' + orig_prefix + '_1000_eigenD/D/' + orig_prefix_window, 'output/' + orig_prefix + '_1000_eigenU/U/' + orig_prefix_window)

        f_out_local = open('output/' + new_prefix_window + '.assoc.txt', 'r')
        f_out_local.readline()
        line = f_out_local.readline()
        while line != '':
            f_assoc_local_window.write(line)
            line = f_out_local.readline()
        f_out_local.close()

        os.remove('output/' + new_prefix_window + '.assoc.txt')
        os.remove('output/' + new_prefix_window + '.log.txt')
        os.remove(new_prefix_window + '.fam')
        os.remove(new_prefix_window + '.bim')
        os.remove(new_prefix_window + '.bed')
        
    f_assoc_local_window.close()
    """
    #####
    # t-test
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
