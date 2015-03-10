import os
import random
from helpers import *

prefix = 'subset_302_homozygotes'

tag_inds = [int(x) for x in open(prefix + '_tag_inds.txt', 'r').readline().split()]
ped_lines = [line.split() for line in open(prefix + '.ped', 'r').readlines()]
map_lines = [line.split() for line in open(prefix + '.map', 'r').readlines()]

random_inds = random.sample(range(73053), 100)

gen_eigen_decomps_linkage_groups(prefix, ped_lines, map_lines, random_inds, 'subset_302_homozygotes_snps_tagged.txt')

random_inds_into_orig = [tag_inds[i] for i in random_inds]
gen_eigen_decomps(prefix, ped_lines, map_lines, random_inds_into_orig, 1000)

f_inds = open('output/compare_kinships/random_inds.txt', 'w')
for i in random_inds:
    f_inds.write(str(i) + ' ')
f_inds.write('\n')
f_inds.close()
