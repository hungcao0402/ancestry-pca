import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
from optparse import OptionParser
import numpy as np
from sklearn.neighbors import NearestNeighbors

# def ancestry_pca(p=None, sp=None,
#                 mapping_data='database/mapping_population.csv',
#                 data_vcf='database/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz',
#                 output_temp='output_temp',
#                 output_figure='output_figure'):
#     # check p and sp
#     if p == None and sp == None:
#         raise(Exception('No population was selected'))
#     elif p != None and sp != None:
#         raise(Exception('Choose only population (-p) or super population (-s)'))
    
#     # extract sub-samples in mapping_population to temporary_sample.txt
#     mapping_population = pd.read_csv(mapping_data)
#     if p != None:
#         samples = mapping_population[mapping_population['Population'].isin(p.split())]['Individual ID']
#     else:
#         samples = mapping_population[mapping_population['Super population'].isin(sp.split())]['Individual ID']
    
#     if samples.shape[0] == 0:
#         raise(Exception('Population or super population not found'))
    
#     samples.to_csv(f'{output_temp}/temporary_sample.txt', index=None, header=False)
    
#     # run bcftools to extract vcf
#     print("Extracting subsamples ...")
#     os.system(command=f'bcftools view -S {output_temp}/temporary_sample.txt {data_vcf} -Oz -o {output_temp}/temporary_sample.vcf.gz')

#     # run plink pca with .vcf file to obtain eigenvec and eigenval files
#     os.system(command=f'plink --vcf {output_temp}/temporary_sample.vcf.gz \
#                             --maf 0.05 --geno 0.05 --mind 0.05 \
#                             --pca 10 --out {output_temp}/temporary_pca')
    
#     # read eigenvec and eigenval files
#     eigenvec = pd.read_table(f'{output_temp}/temporary_pca.eigenvec', header=None, sep=' ')
#     eigenval = pd.read_table(f'{output_temp}/temporary_pca.eigenval', header=None)

#     # calculate percentile of eigen values (% variance explained)
#     eigen_per = round(eigenval.iloc[:, 0] / sum(eigenval.iloc[:, 0]) *100, 2)

#     # mapping population
#     eigenvec = eigenvec.loc[:, [1, 2, 3]].rename(columns={1:'Individual ID', 2:'PC1', 3:'PC2'})
#     eigenvec = pd.merge(eigenvec, mapping_population)

#     # plot
#     ax = plt.figure(figsize=(8, 6))
#     if p != None:
#         output_file = '_'.join(p.split())
#         sns.scatterplot(x='PC1', y='PC2', data=eigenvec, hue='Population', s=12).set(
#                         title='Population with PC1 and PC2', 
#                         xlabel=f'PC1 ({eigen_per[0]}% variance explained)',
#                         ylabel=f'PC2 ({eigen_per[1]}% variance explained)')
#     else:
#         output_file = '_'.join(sp.split())
#         sns.scatterplot(x='PC1', y='PC2', data=eigenvec, hue='Super population', s=12).set(
#                         title='Super population with PC1 and PC2', 
#                         xlabel=f'PC1 ({eigen_per[0]}% variance explained)',
#                         ylabel=f'PC2 ({eigen_per[1]}% variance explained)')
    
#     plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
#     # plt.show()

#     ax.savefig(f'{output_figure}/{output_file}' + '.png', dpi=150, bbox_inches='tight')

# fix by Sang 22/08/23
def percentile_ancestry(eigenvecs, user, p=None, sp=None):
    if user == None:
        return eigenvecs, None, None
    eigenvec = eigenvecs.copy()
    idx = []
    for id in range(eigenvec.shape[0]):
        if eigenvec.loc[id, 'Individual ID'] in user.split():
            eigenvec.loc[id, 'Population'] = eigenvec.loc[id, 'Individual ID']
            eigenvec.loc[id, 'Super population'] = eigenvec.loc[id, 'Individual ID']
            idx.append(id)

    users = eigenvec.iloc[idx, :]
    eigenvec = eigenvec.drop(idx).reset_index(drop=True)
    eigenvec_array = np.array(eigenvec.loc[:, ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']])

    # KNN with K = 100
    neigh = NearestNeighbors(n_neighbors=100)
    neigh.fit(eigenvec_array)
    id_neigh = neigh.kneighbors(np.array(users.loc[idx, ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']]), return_distance=False)

    output = []

    # check ancestry
    for i, id in enumerate(idx):
        ancestry_array = eigenvec.loc[id_neigh[i], ['Population Details', 'Super population Details']]
        string_pop = []; string_sup = []
        a = ancestry_array['Population Details'].value_counts()
        b = ancestry_array['Super population Details'].value_counts()
        for i, x in enumerate(a):
            string_pop.append(str(a.index[i]) + ' ' + str(x))
        for i, x in enumerate(b):
            string_sup.append(str(b.index[i]) + ' ' + str(x))
        string_pop = '%, '.join(string_pop) + '%'
        string_sup = '%, '.join(string_sup) + '%'
        users.loc[id, 'Population Details'] = string_pop
        users.loc[id, 'Super population Details'] = string_sup
        
        if p != None:
            print('Ancestry of sample ' + users.loc[id, 'Individual ID']  + ': ' + string_pop)
            output.append('Ancestry of sample ' + users.loc[id, 'Individual ID']  + ': ' + string_pop)
        else:
            print('Ancestry of sample', users.loc[id, 'Individual ID'], ': ', string_sup)
            output.append('Ancestry of sample ' + users.loc[id, 'Individual ID']  + ': ' + string_sup)
    
    # eigenvec = pd.concat([eigenvec, users], axis=0, ignore_index=True)
    return eigenvec, users, output

def ancestry_pca(p=None, sp=None, user=None,
                mapping_data='database/mapping_population.csv',
                data_vcf='database/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz',
                output_temp='output_temp',
                output_figure='output_figure'):
    # check p and sp
    if p == None and sp == None:
        raise(Exception('No population was selected'))
    elif p != None and sp != None:
        raise(Exception('Choose only population (-p) or super population (-sp)'))
    
    # extract sub-samples in mapping_population to temporary_sample.txt
    mapping_population = pd.read_csv(mapping_data)
    if p != None:
        samples = mapping_population[mapping_population['Population'].isin(p.split())]['Individual ID']
    else:
        samples = mapping_population[mapping_population['Super population'].isin(sp.split())]['Individual ID']
    
    if samples.shape[0] == 0:
        raise(Exception('Population or super population not found'))

    # extract user(s) (one or multiple individual IDs)
    if user != None:
        samples = pd.concat([samples, pd.Series(user.split())], axis=0, ignore_index=True).drop_duplicates(keep='first', ignore_index=True)

    # write samples to txt file
    samples.to_csv(f'{output_temp}/temporary_sample.txt', index=None, header=False)
    
    # run bcftools to extract vcf
    print('Extracting samples...\n')
    os.system(command=f'bcftools view -S {output_temp}/temporary_sample.txt {data_vcf} -Oz -o {output_temp}/temporary_sample.vcf.gz')

    # run plink pca with .vcf file to obtain eigenvec and eigenval files
    print(f'Samples extracted in {output_temp}/temporary_sample.vcf.gz.\nRun PCA with plink tool...\n')
    os.system(command=f'plink --vcf {output_temp}/temporary_sample.vcf.gz \
                            --maf 0.05 --geno 0.05 --mind 0.05 \
                            --pca 10 --out {output_temp}/temporary_pca')
    
    # read eigenvec and eigenval files
    eigenvec = pd.read_table(f'{output_temp}/temporary_pca.eigenvec', header=None, sep=' ')
    eigenval = pd.read_table(f'{output_temp}/temporary_pca.eigenval', header=None)

    # calculate percentile of eigen values (% variance explained)
    eigen_per = round(eigenval.iloc[:, 0] / sum(eigenval.iloc[:, 0]) *100, 2)

    # mapping population
    eigenvec = eigenvec.loc[:, [1, 2, 3, 4, 5, 6]].rename(columns={1:'Individual ID', 2:'PC1', 3:'PC2', 4:'PC3', 5:'PC4', 6:'PC5'})
    eigenvec = pd.merge(eigenvec, mapping_population)

    # check percentile ancestry
    eigenvec, users, output = percentile_ancestry(eigenvecs=eigenvec, user=user, p=p, sp=sp)

    # plot
    ax = plt.figure(figsize=(8, 6))
    if p != None:
        output_file = '_'.join(p.split())
        len_color1 = eigenvec['Population'].nunique()
        sns.scatterplot(x='PC1', y='PC2', data=eigenvec, hue='Population', s=12, palette=sns.color_palette('pastel')[:len_color1])
        
        if user != None:
            len_color2 = users['Population'].nunique()
            sns.scatterplot(x='PC1', y='PC2', data=users, hue='Population', s=150, 
                                palette=sns.color_palette('bright')[len_color1:len_color1+len_color2], marker='*').set(
                                title='Population with PC1 and PC2', 
                                xlabel=f'PC1 ({eigen_per[0]}% variance explained)',
                                ylabel=f'PC2 ({eigen_per[1]}% variance explained)')

    else:
        output_file = '_'.join(sp.split())
        len_color1 = eigenvec['Super population'].nunique()
        sns.scatterplot(x='PC1', y='PC2', data=eigenvec, hue='Super population', s=12, palette=sns.color_palette('pastel')[:len_color1])
        
        if user != None:
            len_color2 = users['Super population'].nunique()
            sns.scatterplot(x='PC1', y='PC2', data=users, hue='Super population', s=150, 
                                palette=sns.color_palette('bright')[len_color1:len_color1+len_color2], marker='*').set(
                                title='Super population with PC1 and PC2', 
                                xlabel=f'PC1 ({eigen_per[0]}% variance explained)',
                                ylabel=f'PC2 ({eigen_per[1]}% variance explained)')

    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    # plt.show()
    ax.savefig(f'{output_figure}/{output_file}' + '.png', dpi=150, bbox_inches='tight')

    return plt.gcf(), output
 
    
# main
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-p", "--population", type=str, default=None,
                    help='Population to plot')
    parser.add_option('-s', '--superpopulation', type=str, default=None,
                    help='Super population to plot')
    parser.add_option('-u', '--user', type=str, default=None,
                    help='Users to check ancestry')

    (opt, args) = parser.parse_args()
    a, b = ancestry_pca(p=opt.population, sp=opt.superpopulation, user=opt.user)
    print(b)