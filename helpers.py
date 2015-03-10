import sys
import os
import math

def extract_ped_indices(ped_lines, map_lines, new_prefix, inds, make_homo=True):
    '''makes a new ped file named new_prefix.ped with only snp indices
    present for each strain in ped_lines; also creates corresponding
    map file and bed/bim files'''
    
    f_ped_subset = open(new_prefix + '.ped', 'w')
    
    for strain in range(len(ped_lines)):
        for info in ped_lines[strain][:6]:
            f_ped_subset.write(info + '  ')
        for i in inds:
            a1 = ped_lines[strain][6 + 2 * i]
            if make_homo:
                f_ped_subset.write(a1 + ' ' + a1 + '  ')
            else:
                a2 = ped_lines[strain][7 + 2 * i]
                f_ped_subset.write(a1 + ' ' + a2 + '  ')
        f_ped_subset.write('\n')

    f_ped_subset.close()

    f_map_subset = open(new_prefix + '.map', 'w')
    for i in inds:
        for info in map_lines[i]:
            f_map_subset.write(info + ' ')
        f_map_subset.write('\n')
    f_map_subset.close()

    os.system('plink --noweb --file ' + new_prefix + ' --make-bed --out ' + new_prefix)

def copy_fam_nonmissing(orig_prefix, new_prefix):
    '''copies fam file to a new file with all missing values (-9)
    replaced by 1'''

    f_fam = open(orig_prefix + '.fam', 'r')
    fam = f_fam.read().replace('-9', '1')
    f_fam.close()
    f_fam_new = open(new_prefix + '.fam', 'w')
    f_fam_new.write(fam)
    f_fam_new.close()

def gemma(new_prefix, kinship_prefix_d, kinship_prefix_u, append = ''):
    '''runs gemma with specified kinship eigen decomposition files'''

    os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + new_prefix + ' -d ' + kinship_prefix_d + '.eigenD.txt' + ' -u ' + kinship_prefix_u + '.eigenU.txt -lmm 4 -notsnp -o ' + new_prefix)

    if append != '':
        os.system('mv output/' + new_prefix + '.assoc.txt output/' + new_prefix + append + '.assoc.txt')
        os.system('mv output/' + new_prefix + '.log.txt output/' + new_prefix + append + '.log.txt')

def filter_ped(prefix, ped_lines, freq_cutoff = .05, \
               freq_upper_cutoff = 1, snp_inds = None, suffix = ''):
    '''creates a file that lists snp indices in ped_lines that are
    biallelic and have MAF above freq_cutoff'''

    num_strains = len(ped_lines)
    num_snps = (len(ped_lines[0]) - 6)/2
    all_alleles = [[] for n in range(num_snps)]
    for strain in range(num_strains):
        alleles = ped_lines[strain][6::2]
        for i in range(num_snps):
            all_alleles[i].append(alleles[i])
    keep_inds = []    
    f_fil = open(prefix + '_filtered_inds' + suffix + '.txt', 'w')
    if snp_inds==None:
        snp_inds = range(num_snps)
    for i in snp_inds:
        alleles = all_alleles[i] 
        alleles_set = list(set(alleles))
        maf = float(alleles.count(alleles_set[0])) / num_strains
        if maf > .5:
            maf = 1 - maf
        if len(alleles_set) == 2 and \
                maf >= freq_cutoff and \
                maf < freq_upper_cutoff:
            keep_inds.append(i)
            f_fil.write(str(i) + ' ')
    print(len(keep_inds))
    f_fil.close()

def phenotypes_to_fam(orig_prefix, new_prefix, phens):
    '''writes phenotypes to a fam file, with strains from
    orig_prefix.fam'''

    f_fam = open(orig_prefix + '.fam', 'r')
    f_fam_new = open(new_prefix + '.fam', 'w')
    line = f_fam.readline()
    i = 0
    while line != '':
        # chop off 1\n and replace with simulated phenotype (we're
        # assuming fam file has strains in same order as file used above)
        assert(line[-3:-1] == ' 1')
        f_fam_new.write(line[:-2] + str(phens[i]) + '\n')
        line = f_fam.readline()
        i += 1
    f_fam.close()
    f_fam_new.close()


def gen_bed_chr(prefix, ped_lines, map_lines, inds):
    '''generates a bed/bim file with only marker indices, one file for
    each chromosome'''
    
    os.system('mkdir bed_chr')
    
    num_markers = len(inds)
    start_ind = 0
    end_ind = 0
    for chrm in range(1, 17):
        chrm = str(chrm)
        prefix_chr = prefix + '_chr' + chrm

        while end_ind < num_markers and map_lines[inds[end_ind]][0] == chrm:
            end_ind += 1

        extract_ped_indices(ped_lines, map_lines, prefix_chr, inds[start_ind:end_ind])

        os.remove(prefix_chr + '.ped')
        os.remove(prefix_chr + '.fam')
        os.remove(prefix_chr + '.map')
        os.remove(prefix_chr + '.nosex')
        os.remove(prefix_chr + '.log')
        os.system('mv ' + prefix_chr + '.bed bed_chr/')
        os.system('mv ' + prefix_chr + '.bim bed_chr/')

        start_ind = end_ind

def gen_bed_af(prefix, ped_lines, map_lines, inds, af, num_bins):
    '''generates a bed/bim file with only marker indices, one file for
    each allele frequency bin'''
    
    os.system('mkdir bed_af')
    
    num_markers = len(inds)

    lower_bound = 1
    for a in af:
        if a != -1 and a < lower_bound:
            lower_bound = a
    bin_size = float(max(af) - lower_bound) / num_bins
    bins = [[] for b in range(num_bins)]
    for i in inds:
        if af[i] >= 0:
            bin_ind = int(math.floor((af[i] - lower_bound) / bin_size))
            if bin_ind == num_bins:
                assert(af[i] == .5)
                bin_ind -= 1
            bins[bin_ind].append(i)
    
    for bin_ind in range(num_bins):
        
        prefix_bin = prefix + '_af_bin_' + str(bin_ind * bin_size + lower_bound)
        print(prefix_bin)
        sys.stdout.flush()

        extract_ped_indices(ped_lines, map_lines, prefix_bin, bins[bin_ind])

        os.remove(prefix_bin + '.ped')
        os.remove(prefix_bin + '.fam')
        os.remove(prefix_bin + '.map')
        os.remove(prefix_bin + '.nosex')
        os.remove(prefix_bin + '.log')
        os.system('mv ' + prefix_bin + '.bed bed_af/')
        os.system('mv ' + prefix_bin + '.bim bed_af/')

def find_minor_allele_freqs(prefix, ped_lines):
    num_sites = (len(ped_lines[0]) - 6) / 2
    num_strains = len(ped_lines)
    f_out = open(prefix + '_maf.txt', 'w')
    for i in range(num_sites):
        a = []
        for line in ped_lines:
            a.append(line[6 + 2 * i])
        a_set = list(set(a))
        if len(a_set) != 2:
            f_out.write('-1 N\n')
        else:
            maf = float(a.count(a_set[0])) / num_strains
            if maf <= .5:
                f_out.write(str(maf) + ' ' + a_set[0] + '\n')
            else:
                f_out.write(str(1 - maf) + ' ' + a_set[1] + '\n')
    f_out.close()

def find_derived_allele_freqs(prefix, ped_lines, ancestral):
    num_sites = (len(ped_lines[0]) - 6) / 2
    num_strains = len(ped_lines)
    f_out = open(prefix + '_daf.txt', 'w')
    for i in range(num_sites):
        a = []
        for line in ped_lines:
            a.append(line[6 + 2 * i])
        a_set = list(set(a))
        if len(a_set) != 2:
            f_out.write('-1 N\n')
        else:
            a_der = a_set[0]
            if a_der == ancestral[i]:
                a_der = a_set[1]
            der_af = float(a.count(a_der)) / num_strains
            f_out.write(str(der_af) + ' ' +  a_der + '\n')
    f_out.close()

def gen_eigen_decomps_af(prefix, ped_lines, map_lines, af, num_bins):
    '''generates eigen decompositions of kinship matrices for each
    allele frequency bin'''

    num_strains = len(ped_lines)
    num_snps = (len(ped_lines[0]) - 6)/2
    assert(len(map_lines) == num_snps)

    lower_bound = 1
    for a in af:
        if a != -1 and a < lower_bound:
            lower_bound = a
    bin_size = float(max(af) - lower_bound) / num_bins
    bins = [[] for b in range(num_bins)]
    for i in range(num_snps):
        if af[i] >= 0:
            bin_ind = int(math.floor((af[i] - lower_bound) / bin_size))
            if bin_ind == num_bins: # happens when we're right at af = .5
                assert(af[i] == .5)
                bin_ind -= 1
            bins[bin_ind].append(i)
    
    for bin_ind in range(num_bins):
        
        print(bin_ind, bin_size, lower_bound)
        prefix_bin = prefix + '_af_bin_' + str(bin_ind * bin_size + lower_bound)
        print(prefix_bin)
        sys.stdout.flush()

        extract_ped_indices(ped_lines, map_lines, prefix_bin, bins[bin_ind])
        
        # assume prefix.fam exists and has all non-missing phenotypes
        os.system('cp ' + prefix + '.fam ' +  prefix_bin + '.fam')

        # generate centered relatedness matrix
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_bin + ' -gk 1 -notsnp -o ' + prefix_bin)

        # eigen decomposition
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_bin + ' -k output/' + prefix_bin + '.cXX.txt -eigen -notsnp -o ' + prefix_bin)

        os.system('mv output/' + prefix_bin + '.eigenD.txt output/D/')
        os.system('mv output/' + prefix_bin + '.eigenU.txt output/U/') 
        
        # remove files we don't need anymore
        os.remove(prefix_bin + '.ped')
        os.remove(prefix_bin + '.map')
        os.remove(prefix_bin + '.bed')
        os.remove(prefix_bin + '.bim')
        os.remove(prefix_bin + '.fam')
        os.remove(prefix_bin + '.log')
        os.remove(prefix_bin + '.nosex')
        os.remove('output/' + prefix_bin + '.cXX.txt')
        os.remove('output/' + prefix_bin + '.log.txt')

        
def gen_eigen_decomps_exclude_chr(prefix, ped_lines, map_lines):
    '''generates eigen decompositions of kinship matrices for each
    chromosome from all snps _not_ on that chromosome'''

    num_strains = len(ped_lines)
    num_snps = (len(ped_lines[0]) - 6)/2
    assert(len(map_lines) == num_snps)

    start_ind = 0
    end_ind = 0
    for chrm in range(1, 17):
        chrm = str(chrm)
        print(chrm)
        sys.stdout.flush()
        prefix_chr = prefix + '_chr' + chrm

        assert(map_lines[start_ind][0] == chrm)
        while end_ind < num_snps and map_lines[end_ind][0] == chrm:
            end_ind += 1
        assert(map_lines[end_ind-1][0] == '16' or map_lines[end_ind][0] == str(int(chrm) + 1))
        inds = range(0, start_ind) + range(end_ind, num_snps)
        
        extract_ped_indices(ped_lines, map_lines, prefix_chr, inds)
        
        # assume prefix.fam exists and has all non-missing phenotypes
        os.system('cp ' + prefix + '.fam ' +  prefix_chr + '.fam')

        # generate centered relatedness matrix
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_chr + ' -gk 1 -notsnp -o ' + prefix_chr)

        # eigen decomposition
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_chr + ' -k output/' + prefix_chr + '.cXX.txt -eigen -notsnp -o ' + prefix_chr)

        os.system('mv output/' + prefix_chr + '.eigenD.txt output/D/')
        os.system('mv output/' + prefix_chr + '.eigenU.txt output/U/') 
        
        # remove files we don't need anymore
        os.remove(prefix_chr + '.ped')
        os.remove(prefix_chr + '.map')
        os.remove(prefix_chr + '.bed')
        os.remove(prefix_chr + '.bim')
        os.remove(prefix_chr + '.fam')
        os.remove(prefix_chr + '.log')
        os.remove(prefix_chr + '.nosex')
        os.remove('output/' + prefix_chr + '.cXX.txt')
        os.remove('output/' + prefix_chr + '.log.txt')

        start_ind = end_ind
    


def gen_eigen_decomps_chr(prefix, ped_lines, map_lines):
    '''generates eigen decompositions of kinship matrices for each
    chromosome'''

    num_strains = len(ped_lines)
    num_snps = (len(ped_lines[0]) - 6)/2
    assert(len(map_lines) == num_snps)

    start_ind = 0
    end_ind = 0
    for chrm in range(1, 17):
        chrm = str(chrm)
        print(chrm)
        sys.stdout.flush()
        prefix_chr = prefix + '_chr' + chrm

        assert(map_lines[start_ind][0] == chrm)
        while end_ind < num_snps and map_lines[end_ind][0] == chrm:
            end_ind += 1
        assert(map_lines[end_ind-1][0] == '16' or map_lines[end_ind][0] == str(int(chrm) + 1))
        inds = range(start_ind, end_ind)
        
        extract_ped_indices(ped_lines, map_lines, prefix_chr, inds)
        
        # assume prefix.fam exists and has all non-missing phenotypes
        os.system('cp ' + prefix + '.fam ' +  prefix_chr + '.fam')

        # generate centered relatedness matrix
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_chr + ' -gk 1 -notsnp -o ' + prefix_chr)

        # eigen decomposition
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_chr + ' -k output/' + prefix_chr + '.cXX.txt -eigen -notsnp -o ' + prefix_chr)

        os.system('mv output/' + prefix_chr + '.eigenD.txt output/D/')
        os.system('mv output/' + prefix_chr + '.eigenU.txt output/U/') 
        
        # remove files we don't need anymore
        os.remove(prefix_chr + '.ped')
        os.remove(prefix_chr + '.map')
        os.remove(prefix_chr + '.bed')
        os.remove(prefix_chr + '.bim')
        os.remove(prefix_chr + '.fam')
        os.remove(prefix_chr + '.log')
        os.remove(prefix_chr + '.nosex')
        os.remove('output/' + prefix_chr + '.cXX.txt')
        os.remove('output/' + prefix_chr + '.log.txt')

        start_ind = end_ind

def gen_eigen_decomps_linkage_groups(prefix, ped_lines, map_lines, group_inds, fn_tag, n = 2000):
    '''generates eigen decompositions of kinship matrices for each
    given snp index based on other snps in the linkage group'''

    num_strains = len(ped_lines)
    num_snps = (len(ped_lines[0]) - 6)/2

    f_tag = open(fn_tag, 'r')
    ld_groups = {}
    line = f_tag.readline()
    while line != '':
        line = line.split()
        tag = line[3]
        snp = line[2]
        if tag in ld_groups:
            ld_groups[tag].append(snp)
        else:
            ld_groups[tag] = [snp]
        line = f_tag.readline()
    f_tag.close()

    rsid_to_ind = {}
    for i in range(len(map_lines)):
        rsid_to_ind[map_lines[i][1]] = i
    
    #for tag in ld_groups:
    if group_inds == None:
        group_inds = range(len(ld_groups.keys()))
    for group_ind in group_inds:
        tag = ld_groups.keys()[group_ind]
        group = ld_groups[tag]
        window = int(round(n / float(len(group)) / 2))
        prefix_rsid = prefix + '_' + tag
        target_chr = map_lines[rsid_to_ind[tag]][0]
        inds = []
        for snp in group:
            i = rsid_to_ind[snp]
            start_ind = max(0, i - window)
            while map_lines[start_ind][0] != target_chr:
                start_ind += 1
            # note that we'd add one here to get an even window on
            # either side but we'll do it this way to keep it more
            # even for different group sizes
            end_ind = min(i + window, num_snps)
            while map_lines[end_ind - 1][0] != target_chr:
                end_ind -= 1
            inds += range(start_ind, end_ind)
        print('*** ' + tag + ' ' + str(len(group)) + ' ' + str(len(inds)))
        # sys.stdout.flush()
        extract_ped_indices(ped_lines, map_lines, prefix_rsid, inds)
        
        # assume prefix.fam exists and has all non-missing phenotypes
        os.system('cp ' + prefix + '.fam ' +  prefix_rsid + '.fam')
        
        # generate centered relatedness matrix
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_rsid + ' -gk 1 -notsnp -o ' + prefix_rsid)

        # eigen decomposition
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_rsid + ' -k output/' + prefix_rsid + '.cXX.txt -eigen -notsnp -o ' + prefix_rsid)

        os.system('mv output/' + prefix_rsid + '.eigenD.txt output/compare_kinships/ld/')
        os.system('mv output/' + prefix_rsid + '.eigenU.txt output/compare_kinships/ld/') 
        
        # remove files we don't need anymore
        os.remove(prefix_rsid + '.ped')
        os.remove(prefix_rsid + '.map')
        os.remove(prefix_rsid + '.bed')
        os.remove(prefix_rsid + '.bim')
        os.remove(prefix_rsid + '.fam')
        os.remove(prefix_rsid + '.log')
        os.remove(prefix_rsid + '.nosex')
        #os.remove('output/' + prefix_rsid + '.cXX.txt')
        os.system('mv output/' + prefix_rsid + '.cXX.txt output/compare_kinships/ld/')
        os.remove('output/' + prefix_rsid + '.log.txt')
        
    #eigenD_zip = prefix + '_ld_' + str(window) + '.eigenD.zip'
    #eigenU_zip = prefix + '_ld_' + str(window) + '.eigenU.zip'
    #os.chdir('output')
    #os.system('zip -r ' + eigenD_zip + ' D/')
    #os.system('rm -r D')
    #os.system('zip -r ' + eigenU_zip + ' U/')
    #os.system('rm -r U')
    #os.chdir('..')


def gen_eigen_decomps(prefix, ped_lines, map_lines, inds, window = 0):
    '''generates eigen decompositions of kinship matrices for each
    given snp index with specified window of snps around it'''

    num_strains = len(ped_lines)
    num_snps = (len(ped_lines[0]) - 6)/2

    os.system('mkdir output/D_window_' + str(window))
    os.system('mkdir output/U_window_' + str(window))
    os.system('mkdir output/kinship_window_' + str(window))

    for i in inds:
        rsid = map_lines[i][1]
        prefix_rsid = prefix + '_' + rsid
        
        # create a ped (and bim/bed/map) file with only the snps in
        # specified window present, not crossing chromosome boundaries
        start_ind = max(0, i - window)
        target_chr = map_lines[i][0]
        while map_lines[start_ind][0] != target_chr:
            start_ind += 1
        end_ind = min(i + window + 1, num_snps)
        while map_lines[end_ind-1][0] != target_chr:
            end_ind -= 1
        extract_ped_indices(ped_lines, map_lines, prefix_rsid, range(start_ind, end_ind))

        # assume prefix.fam exists and has all non-missing phenotypes
        os.system('cp ' + prefix + '.fam ' +  prefix_rsid + '.fam')

        # generate centered relatedness matrix
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_rsid + ' -gk 1 -notsnp -o ' + prefix_rsid)

        # eigen decomposition
        os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix_rsid + ' -k output/' + prefix_rsid + '.cXX.txt -eigen -notsnp -o ' + prefix_rsid)

        os.system('mv output/' + prefix_rsid + '.eigenD.txt output/D_window_' + str(window))
        os.system('mv output/' + prefix_rsid + '.eigenU.txt output/U_window_' + str(window)) 

        # remove files we don't need anymore
        os.remove(prefix_rsid + '.ped')
        os.remove(prefix_rsid + '.map')
        os.remove(prefix_rsid + '.bed')
        os.remove(prefix_rsid + '.bim')
        os.remove(prefix_rsid + '.fam')
        os.remove(prefix_rsid + '.log')
        os.remove(prefix_rsid + '.nosex')
        #os.remove('output/' + prefix_rsid + '.cXX.txt')
        os.system('mv output/' + prefix_rsid + '.cXX.txt output/kinship_window_' + str(window))
        os.remove('output/' + prefix_rsid + '.log.txt')

    #eigenD_zip = prefix + '_' + str(window) + '.eigenD.zip'
    #eigenU_zip = prefix + '_' + str(window) + '.eigenU.zip'
    #os.chdir('output')
    #os.system('zip -r ' + eigenD_zip + ' D_window_' + str(window))
    #os.system('rm -r D')
    #os.system('zip -r ' + eigenU_zip + ' U_window' )
    #os.system('rm -r U')
    #os.chdir('..')

def gen_marker_beds(prefix, ped_lines, map_lines, inds):
    '''create a bed file for each snp in inds with only that snp (for
    testing associations at single snp)'''

    num_strains = len(ped_lines)
    num_snps = (len(ped_lines[0]) - 6)/2

    for i in inds:
        rsid = map_lines[i][1]
        prefix_rsid = prefix + '_' + rsid

        # create ped/bed/bim files with only this snp
        extract_ped_indices(ped_lines, map_lines, prefix_rsid, [i])

        os.system('mv ' + prefix_rsid + '.bed bed/')
        # note that we don't need to worry about keeping bim files
        # around because we can just index into the one we create for
        # the whole genome

        # remove files we don't need anymore
        os.remove(prefix_rsid + '.bim')
        os.remove(prefix_rsid + '.ped')
        os.remove(prefix_rsid + '.fam')
        os.remove(prefix_rsid + '.map')
        os.remove(prefix_rsid + '.log')
        os.remove(prefix_rsid + '.nosex')
