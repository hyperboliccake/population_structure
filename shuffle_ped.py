import random

prefix = 'subset_302_homozygotes'
f_ped = open(prefix + '.ped', 'r')
line = f_ped.readline()
num_sites = int((len(line.split()) - 6) / 2)
alleles_by_site = [[] for site in range(num_sites)]
infos = []
while line != '':
    line = line.split()
    infos.append(line[:6])
    alleles = line[6::2]
    for i in range(num_sites):
        alleles_by_site[i].append(alleles[i])
    line = f_ped.readline()
f_ped.close()

print('shuffling')

alleles_shuffled = []
for a in alleles_by_site:
    random.shuffle(a)
    alleles_shuffled.append(a)

f_ped_shuffled = open(prefix + '_shuffled.ped', 'w')
num_strains = len(alleles_by_site[0])
for strain in range(num_strains):
    print(strain)
    for info in infos[strain]:
        f_ped_shuffled.write(info + '  ')
    for site in range(num_sites):
        a = alleles_shuffled[site][strain]
        f_ped_shuffled.write(a + ' ' + a + '  ')
    f_ped_shuffled.write('\n')
f_ped_shuffled.close()
