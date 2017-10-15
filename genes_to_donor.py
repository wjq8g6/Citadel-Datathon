import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


inter = [[0 for _ in range(8)] for a in range(8)]
tcga = pd.read_csv('tcga.csv')
tcga = tcga.loc[tcga['fpkm_expression'] > 0]
print(tcga.loc[tcga['gene_id'] == 'ENSG00000238191'])
print(len(tcga['gene_id'].unique()))
list_organs = tcga['organ'].unique()
gene_sets = []
for org in list_organs:
    temp = tcga.loc[tcga['organ'] == org]
    ids = temp['gene_id'].unique()
    print("Organ: " + org + "  ", len(ids))
    gene_sets.append(set(ids))

for i in range(8):
    for j in range(8):
        inter[i][j] = len(set.intersection(gene_sets[i],gene_sets[j]))

for row in inter:
    print(row)

diff = []
print(len(set.intersection(*gene_sets)))

sample_exp = pd.read_csv('gtex_sample_expression.csv')
sample_exp = sample_exp.loc[sample_exp['rpkm_expression'] > 0]
temp = sample_exp.loc[sample_exp['gene_id'] == "ENSG00000237467"]
donor_data = pd.read_csv('gtex_donor.csv')
print(temp['sample_id'])

for k in range(8):
    temp = gene_sets[:k] + gene_sets[(k+1):]
    print("Organ " + list_organs[k]+":", (set.difference(gene_sets[k],set.union(*temp))))
    diff.append(set.difference(gene_sets[k],set.union(*temp)))

donors_gene = [dict() for p in range(8)]
for j in range(8):
    print(list_organs[j])
    dict = donors_gene[j]
    for gene in diff[j]:
        temp = sample_exp.loc[sample_exp['gene_id'] == gene]
        tissues = set()
        for samp in temp['sample_id'].values:
            lst = samp.split('-')
            donor = '-'.join(lst[:2])
            tissues.add(donor)
        te = {}
        male = 0
        female = 0
        age = 0.0
        ind = donor_data.index[donor_data['donor'].isin(tissues)].tolist()
        for don in ind:
            gend = donor_data.get_value(don, 'gender')
            ag = donor_data.get_value(don, 'age')
            age += ag
            if gend == 'M':
                male += 1
            else:
                female += 1

        if (male + female) > 0:
            age = age / (male + female)
        te['male'] = male
        te['female'] = female
        te['age'] = age
        dict[gene] = te

for i in range(8):
    print("-----" + list_organs[i] + '-----')
    count = 1
    fig = plt.figure()
    fig.suptitle('Organ: ' + list_organs[i])
    for key in donors_gene[i].keys():
        ax = plt.subplot(len(donors_gene[i].keys()),1, count)
        dict = donors_gene[i].get(key)
        print("Gene : " + key)
        print("Male : ",dict['male'])
        print("Female : ", dict['female'])
        print("Age : ", dict['age'])
        ax.bar([0,1], [dict['male'],dict['female']], align='center')
        ax.set_xticklabels(['', '', 'Male','','','','Female'])
        ax.set_ylabel('Frequency')
        ax.set_title("Gene: " + key)
        count+= 1
    plt.show()

