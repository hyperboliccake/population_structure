map_lines = [line.split() for line in open('subset_302_homozygotes.map', 'r').readlines()]
loc_to_rsid = {}
rsid_to_ind = {}
for i in range(len(map_lines)):
    line = map_lines[i]
    loc_to_rsid[(line[0], line[3])] = line[1]
    rsid_to_ind[line[1]] = i
chrs = [str(x) for x in range(1, 17)]

f_snps_tagged = open('subset_302_homozygotes_snps_tagged.txt', 'w')
f_map_tag = open('subset_302_homozygotes_tag.map', 'w')
f_tag_inds = open('subset_302_homozygotes_tag_inds.txt', 'w')

for c in chrs:
    print(c)
    f = open('/net/gs/vol1/home/aclark4/population_structure/code/output/ldSelect/ldSelect_chr_' + c + '.txt', 'r')
    line = f.readline()
    while line != '':
        if line[:3] == 'Bin':
            #print('Bin ' + line.split()[1])
            tag_snps = f.readline().split()[3:]
            chosen_tag = tag_snps[int(len(tag_snps)/2)]
            chosen_id = loc_to_rsid[(c, chosen_tag)]
            f_map_tag.write(c + ' ' + chosen_id + ' 0 ' + chosen_tag + '\n')
            tag_inds.append(rsid_to_ind[chosen_id])
            all_snps = tag_snps + f.readline().split()[3:]
            for s in all_snps:
                f_snps_tagged.write(c + ' ' + s + ' ')
                f_snps_tagged.write(loc_to_rsid[(c, s)] + ' ' + chosen_id + '\n')
        line = f.readline()
    f.close()

f_snps_tagged.close()
f_map_tag.close()

tag_inds.sort()
for i in tag_inds:
    f_tag_inds.write(str(i) + ' ')
f_tag_inds.write('\n')
f_tag_inds.close()

