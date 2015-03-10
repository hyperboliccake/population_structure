import sys
import os
import random
import copy
import scipy.stats
import shutil
from helpers import *

def main():

    freq_group = sys.argv[1]
    num_qtls = int(sys.argv[2])
    h2 = [float(x) for x in sys.argv[3:(3 + num_qtls)]]
    sim_id = sys.argv[3 + num_qtls]
    orig_prefix = 'subset_302_homozygotes'
    new_prefix = orig_prefix + '_' + freq_group + '_' + str(num_qtls)
    for i in range(num_qtls):
        new_prefix += '_' + str(h2[i])
    new_prefix += '_' + sim_id

    print('choosing causal snp')
    # choose random snp from all snps that meet filter criteria (this is
    # more than just the tag snps)

    f_causal = open('output/' + new_prefix + '.causal.txt', 'w')
    f_filtered = open(orig_prefix + '_filtered_inds_' + freq_group + '.txt', 'r')
    snp_inds = f_filtered.readline().split()
    f_filtered.close()
    random_snp_inds = []
    for i in range(num_qtls):
        r = random.choice(snp_inds)
        random_snp_inds.append(int(r)) # index in original ped/map
        f_causal.write(r + ' ')
    f_causal.write('\n')
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
    for i in range(num_qtls):
        alleles_current = []
        for strain in range(num_strains):
            alleles_current.append(ped_lines[strain][6::2][random_snp_inds[i]])
        alleles.append(alleles_current)
    alleles_sets = [list(set(a)) for a in alleles]
    major = [a[0] for a in alleles_sets]
    minor = [a[1] for a in alleles_sets]
    maf = []
    for i in range(num_qtls):
        if alleles[i].count(minor[i]) > alleles.count(major[i]):
            major[i] = alleles_sets[i][1]
            minor[i] = alleles_sets[i][0]
        maf.append(float(alleles[i].count(minor[i])) / len(alleles[i]))

    fixed_effects = [pow(float(h2[i] * (num_strains-1)) / ((h2[i]-1) * num_strains * (maf[i]-1) * maf[i]), .5) for i in range(num_qtls)]

    phens = []
    for s in range(num_strains):
        mu = 0
        for i in range(num_qtls):
            # should we randomly choose between major and minor?
            if alleles[i][s] == minor[i]:
                mu += fixed_effects[i]
        phens.append(random.gauss(mu, 1))

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
    # local kinship matrix by 1000bp window
    #####
    print('- local K by 1000bp window')

    rsids = [int(line.split()[1][2:]) for line in open(orig_prefix + '_tag.map', 'r').readlines()]
    rsids.sort()
    tag_bim_lines = open(orig_prefix + '_tag.bim', 'r').readlines()

    f_assoc_local_window = open('output/' + new_prefix + '_local_window.assoc.txt', 'w')

    for i in range(len(rsids)):
        rsid = 'rs' + str(rsids[i])
        orig_prefix_window = orig_prefix + '_' + rsid
        new_prefix_window = new_prefix + '_' + rsid
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

    #####
    # whole-genome K + admixture
    #####

    print('- whole-genome K with admixture covariates')

    os.symlink(orig_prefix + '_tag.bed', new_prefix + '.bed')
    os.symlink(orig_prefix + '_tag.bim', new_prefix + '.bim')

    kinship_prefix_d = 'output/' + orig_prefix
    kinship_prefix_u = 'output/' + orig_prefix
    os.system('/net/gs/vol1/home/aclark4/software/gemma' + \
                  ' -bfile ' + new_prefix + \
                  ' -d ' + kinship_prefix_d + '.eigenD.txt' + \
                  ' -u ' + kinship_prefix_u + '.eigenU.txt' + \
                  ' -lmm 4' + \
                  ' -notsnp' + \
                  ' -o ' + new_prefix + \
                  ' -c admixture_covariates_8_pops.txt')
    os.system('mv output/' + new_prefix + '.assoc.txt output/' + new_prefix + '_whole_K_admixture.assoc.txt')

    os.remove(new_prefix + '.bed')
    os.remove(new_prefix + '.bim')
    os.remove('output/' + new_prefix + '.log.txt')

    #####
    # admixture
    #####

    print('- admixture covariates only')
    os.system('Rscript sim_admixture.R ' + sim_id + ' ' + new_prefix)

    #####
    # whole-genome K + PCs (first 2)
    #####

    print('- whole-genome K with 2 PC covariates')

    os.symlink(orig_prefix + '_tag.bed', new_prefix + '.bed')
    os.symlink(orig_prefix + '_tag.bim', new_prefix + '.bim')

    kinship_prefix_d = 'output/' + orig_prefix
    kinship_prefix_u = 'output/' + orig_prefix
    os.system('/net/gs/vol1/home/aclark4/software/gemma' + \
                  ' -bfile ' + new_prefix + \
                  ' -d ' + kinship_prefix_d + '.eigenD.txt' + \
                  ' -u ' + kinship_prefix_u + '.eigenU.txt' + \
                  ' -lmm 4' + \
                  ' -notsnp' + \
                  ' -o ' + new_prefix + \
                  ' -c pc_2_all.txt')
    os.system('mv output/' + new_prefix + '.assoc.txt output/' + new_prefix + '_whole_K_pc_2.assoc.txt')
    
    os.remove(new_prefix + '.bed')
    os.remove(new_prefix + '.bim')
    os.remove('output/' + new_prefix + '.log.txt')

    #####
    # PCs (first 2)
    #####

    print('- 2 PC covariates only')
    os.system('Rscript sim_pc_2.R ' + sim_id + ' ' + new_prefix)

    #####
    # K by maf bins
    #####

    print('- K by maf bins')

    f_assoc_af = open('output/' + new_prefix + '_af.assoc.txt', 'w')

    bed_bim_files = sorted(os.listdir('bed_af'))
    for i in range(len(bed_bim_files)/2):
        orig_prefix_af = bed_bim_files[2 * i][:-4]
        new_prefix_af = new_prefix + '_af_bin_' + str(i)

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
