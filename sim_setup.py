import os
from helpers import *

# things that we only need to do once across all simulations

prefix = 'subset_302_homozygotes'

#tag_inds = [int(x) for x in open(prefix + '_tag_inds.txt', 'r').readline().split()]
ped_lines = [line.split() for line in open(prefix + '.ped', 'r').readlines()]
map_lines = [line.split() for line in open(prefix + '.map', 'r').readlines()]

# make ped/bed file with only tag snps to go with tag map file
print('making ped/bed files with only tag snps')
#extract_ped_indices(ped_lines, map_lines, prefix + '_tag', tag_inds)

# filter ped file by MAF and biallecism to determine which snps we can
# pick as causal
print('filtering all snps by allele frequency')
filter_ped(prefix, ped_lines, 3.0/302, 1, None, '_all')
filter_ped(prefix, ped_lines, 3.0/302, .05, None, '_rare')
filter_ped(prefix, ped_lines, .05, .1, None, '_intermediate')
filter_ped(prefix, ped_lines, .1, 1, None, '_common')

print('generating eigen decompositions of kinship matrices for tag snps')
#gen_eigen_decomps(prefix, ped_lines, map_lines, tag_inds, 0)
#gen_eigen_decomps_exclude_chr(prefix, ped_lines, map_lines)

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
