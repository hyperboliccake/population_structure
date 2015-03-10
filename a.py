import os
for i in range(1,32):
    try:
        os.remove('subset_302_homozygotes_' + str(i) + '.fam')
        os.remove('subset_302_homozygotes_' + str(i) + '.bed')
        os.remove('subset_302_homozygotes_' + str(i) + '.bim')
    except:
        pass
