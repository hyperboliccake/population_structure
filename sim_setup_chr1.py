import os
from helpers import *

# things that we only need to do once across all simulations

prefix = 'subset_302_homozygotes'

# indices into original ped/map
tag_inds = [int(x) for x in open(prefix + '_tag_inds.txt', 'r').readline().split()]
ped_lines = [line.split() for line in open(prefix + '.ped', 'r').readlines()]
map_lines = [line.split() for line in open(prefix + '.map', 'r').readlines()]

chr1_inds = []
for i in range(len(map_lines)):
    if map_lines[i][0] == '1':
        chr1_inds.append(i)

chr1_tag_inds = []
for i in tag_inds:
    if map_lines[i][0] == '1':
        chr1_tag_inds.append(i)

f_tag_inds = open(prefix + '_tag_inds_chr1.txt', 'w')
for i in chr1_tag_inds:
    f_tag_inds.write(str(i) + ' ')
f_tag_inds.write('\n')
f_tag_inds.close()

# make ped/bed file with only tag snps to go with tag map file
print('making ped/bed files with only tag snps on chr1')
#extract_ped_indices(ped_lines, map_lines, prefix + '_tag_chr1', chr1_tag_inds)

# filter ped file by MAF and biallecism to determine which snps we can
# pick as causal (this is unrelated to tag snps)
print('filtering all snps by allele frequency')
#filter_ped(prefix, ped_lines, .05, chr1_inds, '_chr1')
"""
print('generating eigen decompositions of kinship matrices for tag snps')
# window method
print('===10')
gen_eigen_decomps(prefix, ped_lines, map_lines, chr1_tag_inds, 10)
print('===100')
gen_eigen_decomps(prefix, ped_lines, map_lines, chr1_tag_inds, 100)
print('===1000')
gen_eigen_decomps(prefix, ped_lines, map_lines, chr1_tag_inds, 1000)
print('===5000')
gen_eigen_decomps(prefix, ped_lines, map_lines, chr1_tag_inds, 5000)
print('===10000')
gen_eigen_decomps(prefix, ped_lines, map_lines, chr1_tag_inds, 10000)

"""
print('===20000')
gen_eigen_decomps(prefix, ped_lines, map_lines, chr1_tag_inds, 20000)
"""
print('generating bed files for tag snps')
#gen_marker_beds(prefix, ped_lines, map_lines, tag_inds)
#gen_bed_chr(prefix, ped_lines, map_lines, tag_inds)

# eigen decomposition of kinship matrix for whole genome (along with
# bed/bim files)
print('generating eigen decomposition of whole-genome kinship matrix')
#os.system('plink --noweb --file ' + prefix + ' --make-bed --out ' + prefix)
#copy_fam_nonmissing(prefix, prefix)
#os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix + ' -gk 1 -notsnp -o ' + prefix)
#os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix + ' -k output/' + prefix + '.cXX.txt -eigen -notsnp -o ' + prefix)


#os.system('cp ' + prefix + '.map ' + prefix + '_shuffled.map')
#os.system('plink --noweb --file ' + prefix + '_shuffled --make-bed --out ' + prefix + '_shuffled')
#copy_fam_nonmissing(prefix, prefix + '_shuffled')
#os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix + '_shuffled' + ' -gk 1 -notsnp -o ' + prefix + '_shuffled')
#os.system('/net/gs/vol1/home/aclark4/software/gemma -bfile ' + prefix + '_shuffled -k output/' + prefix + '_shuffled.cXX.txt -eigen -notsnp -o ' + prefix + '_shuffled')
#ped_lines = [line.split() for line in open(prefix + '_shuffled.ped', 'r').readlines()]
#extract_ped_indices(ped_lines, map_lines, prefix + '_tag_shuffled', tag_inds)
"""
