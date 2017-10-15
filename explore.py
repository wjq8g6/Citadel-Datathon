import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



donor_data = pd.read_csv('gtex_donor.csv')
donor_list = donor_data['donor']
sample_data = pd.read_csv('gtex_sample.csv')
sample_exp = pd.read_csv('gtex_sample_expression.csv')
gene_dict = {}
sample_dict = {}
for donor in donor_list:
    temp = sample_data.loc[sample_data['donor'] == donor]
    sample_list = temp['sample_id'].values
    print(donor)
    for sample in sample_list:
        genes = sample_exp.loc[sample_exp['sample_id'] == sample]
        gene_dict[donor] = genes['gene_id'].values






gene_model = pd.read_csv('gtex_gene_model.csv')
gene_model['gene_id'].unique()
genes = pd.read_csv('genes.csv')
print(len(gene_model['gene_id'].unique()))





death_class = donor_data['death_class'].unique()
death_class.sort()
donors = {}
for i in death_class:
    temp = donor_data.loc[donor_data['death_class'] == i]
    print("Death class ", i)
    print("Count: ", len(temp))
    print("\n")
    donors[i] = temp





breast_data = pd.read_csv('breast_tcga.csv')
print(breast_data['gene_id'].value_counts())
ens_ids = breast_data['gene_id'].unique()
print(len(ens_ids))
breast_genes = genes.loc[genes['ensembl_id'].isin(ens_ids)]
breast_genes.to_csv('breast_genes.csv', index=False)
