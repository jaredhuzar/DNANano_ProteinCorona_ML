#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns


outputDir = 'Results\\'
mappingPath=sys.argv[1]
interactionsPath = sys.argv[2]
dataPresentPath = sys.argv[3]

dataPresent = pd.read_csv(dataPresentPath)
abundance_present = dataPresent['Present']
abundance_present = np.asarray(abundance_present)
overallRes3 = dataPresent[['Present', 'Sample', 'Protein']]
overallRes3=overallRes3.rename(columns={'Present':'Abundance','Protein':'ID'})
dataPresent.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0'],inplace = True, axis =1)

translate=pd.read_csv(mappingPath, sep = '\t')
endoPros=pd.read_csv(interactionsPath, sep = '\t')

endoProteins=translate[translate['preferredName'].isin(list(set(list(endoPros['#node1']) + list(endoPros['node2']))))]
endoProteins.reset_index(drop=True, inplace = True)
endoProteins=endoProteins['queryItem']

samples = list()
endoCountList = list()
for x in np.unique(overallRes3['Sample']):
    samples.append(x)
    endoCount = 0
    for z in endoProteins:
        subRes=overallRes3[overallRes3['ID'] == z]
        subRes=subRes[subRes['Sample'] == x]
        if(subRes['Abundance'].values[0] == 1):
            endoCount = endoCount+1
    endoCountList.append(endoCount)
    
endo=pd.DataFrame({'Sample':samples, 'Count':endoCountList})
endocytosis=endo
endo.index = endo['Sample']

endocytosis=endocytosis.reindex(['BoxPL','I-10PL','I-12 PL','I-14 PL','I-16 PL','tube PL','tetrahedron PL','tetrahedron+cho PL','tetrahedron+cho','tetrahedron','tube','I-16','I-14','I-12','I-10','tile','tensegrity'])
endocytosis.transpose()

hueVec = list()
for z in endocytosis.index:
    if('PL' in z):
        hueVec.append('Modified')
        continue
    elif('cho' in z):
        hueVec.append('Modified')
        continue
    else:
        hueVec.append('Bare')
        
endo=endocytosis.melt(value_vars=['Count'])
endo.iloc[:,0] = endocytosis.index
endo=endo.rename({'variable':'Sample'},axis = 1)
endo=endo.rename({'value':'Protein Count'},axis = 1)
endo['Hue'] = hueVec
endo=endo.rename({'Hue':'Nanostructure'},axis = 1)

fig = plt.figure(figsize=(8, 6))
plt.rc('font', family='serif')
plt.rc('font', family='serif')

ax = sns.barplot(data=endo,y ='Sample',x = 'Protein Count', orient='h',hue = 'Nanostructure')
#plt.axvline(36,linewidth=4, color='r')
plt.title('Endocytic Proteins in Corona')
#plt.yaxis('Enrichment Value')
ax.set_xlabel('Protein Count')
plt.savefig(outputDir + 'EndocyticProteinEnrichment.png',bbox_inches='tight')

endo.to_csv(outputDir + 'EndocyticProteinCount.csv')

