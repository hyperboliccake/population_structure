import sys
import os
import random
import copy
import scipy.stats
import shutil
from helpers import *

@profile
def main():

    orig_prefix = 'subset_302_homozygotes'
    new_prefix = orig_prefix + '_' + sys.argv[1]

    shuffle = False
    if len(sys.argv) > 2 and sys.argv[2] == 'shuffle':
        shuffle = True
        new_prefix += '_shuffled'

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
    phens = []

    for allele in alleles:
        mu = 0
        if allele == minor:
            # fixed effect size - should this depend on minor
            # allele frequency to keep power constant?
            mu = 1 
        phens.append(random.gauss(mu, 1))

    # shuffling for permutation tests!
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
    os.symlink(orig_prefix + '_tag.bed', new_prefix + '.bed')
    os.symlink(orig_prefix + '_tag.bim', new_prefix + '.bim')
    gemma(new_prefix, 'output/' + orig_prefix, 'output/' + orig_prefix, '_whole_K')
    os.remove(new_prefix + '.bed')
    os.remove(new_prefix + '.bim')

    #####
    # local kinship matrix (one site at a time)
    #####
    window = 1000
    tag_inds = [int(x) for x in open(orig_prefix + '_tag_inds.txt', 'r').readline().split()]
    tag_inds.sort()
    rsids = [line.split()[1] for line in open(orig_prefix + '.map', 'r').readlines()]

    # automatically sorted by chr, pos
    bim_lines = open(orig_prefix + '.bim', 'r').readlines() 

    f_assoc_local = open('output/' + new_prefix + '_local.assoc.txt', 'w')
    f_assoc_whole = open('output/' + new_prefix + '_whole_K_single.assoc.txt', 'w')
    f_assoc_t = open('output/' + new_prefix + '_t.assoc.txt', 'w')

    eigenD_tar = orig_prefix + '_' + str(window) + '.eigenD.tar'
    eigenU_tar = orig_prefix + '_' + str(window) + '.eigenU.tar'

    for i in tag_inds:
        rsid = rsids[i]
        orig_prefix_single = orig_prefix + '_' + rsid
        new_prefix_single = orig_prefix_single + '_' + sys.argv[1]
        if shuffle:
            new_prefix_single += '_shuffled'

        # bed and bim files don't depend on window so we can just do this
        # once for all methods we're using
        f_bim_single = open(new_prefix_single + '.bim', 'w')
        bim_line = bim_lines[i]
        assert(bim_line.split()[1] == rsid)
        f_bim_single.write(bim_line)
        f_bim_single.close()
        #os.system('ln -s bed/' + orig_prefix_single + '.bed ' + new_prefix_single + '.bed')
        os.symlink('bed/' + orig_prefix_single + '.bed', new_prefix_single + '.bed')

        os.symlink(new_prefix + '.fam', new_prefix_single + '.fam')

        #####
        # single site kinship matrix
        #####

        # extract eigenvector/value files for this site and rename them
        #os.system('unzip -j output/' + orig_prefix + '_' + str(window) + '.eigenD.zip D/' + orig_prefix_single + '.eigenD.txt -d ' + new_prefix)
        #os.system('unzip -j output/' + orig_prefix + '_' + str(window) + '.eigenU.zip U/' + orig_prefix_single + '.eigenU.txt -d ' + new_prefix)
        #os.symlink('output/' + orig_prefix + '_' + str(window) + '_eigenD/D/' + orig_prefix_single + '.eigenD.txt', new_prefix_single + '.eigenD.txt' )

        #gemma(new_prefix_single, new_prefix + '/' + orig_prefix_single)
        gemma(new_prefix_single, 'output/' + orig_prefix + '_' + str(window) + '_eigenD/D/' + orig_prefix_single, 'output/' + orig_prefix + '_' + str(window) + '_eigenU/U/' + orig_prefix_single)

        f_out_local = open('output/' + new_prefix_single + '.assoc.txt', 'r')
        f_out_local.readline()
        results = f_out_local.readline()
        f_assoc_local.write(results)
        f_out_local.close()

        #####
        # whole genome kinship matrix (confirmed that this comes out identical to doing all sites at once)
        #####

        #gemma(new_prefix_single, 'output/' + orig_prefix)

        #f_out_whole = open('output/' + new_prefix_single + '.assoc.txt', 'r')
        #f_out_whole.readline()
        #results = f_out_whole.readline()
        #f_assoc_whole.write(results)
        #f_out_whole.close()

        #####
        # t-test
        #####

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

        # uncorrected t-test, unshuffled
        p = scipy.stats.ttest_ind(phens_1, phens_2)[1]
        f_assoc_t.write(str(p) + '\n')

        # remove files from all these tests
        os.remove('output/' + new_prefix_single + '.assoc.txt')
        os.remove('output/' + new_prefix_single + '.log.txt')
        #os.remove(new_prefix + '/' + orig_prefix_single + '.eigenD.txt')
        #os.remove(new_prefix + '/' + orig_prefix_single + '.eigenU.txt')
        os.remove(new_prefix_single + '.fam')
        os.remove(new_prefix_single + '.bim')
        os.remove(new_prefix_single + '.bed')


    f_assoc_local.close()
    f_assoc_whole.close()
    f_assoc_t.close()
    os.remove(new_prefix + '.fam')
    #os.rmdir(new_prefix)


    #####
    # emma (whole-genome kinship matrix)
    #####

main()
