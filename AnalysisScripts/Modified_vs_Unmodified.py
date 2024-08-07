#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted

datapresentPath=sys.argv[1]
outputDir = 'Results\\'

dataPresent = pd.read_csv(datapresentPath)
abundance_present = dataPresent['Present']
abundance_present = np.asarray(abundance_present)
overallRes3 = dataPresent[['Present', 'Sample', 'Protein']]
dataPresent.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0'],inplace = True, axis =1)

dataNew=dataPresent
abundance_binary=abundance_present

overallRes3=overallRes3.rename(columns={'Present':'Abundance','Protein':'ID'})

polCoated = dataNew[dataNew['PolyL'] == 1]
bare = dataNew[dataNew['PolyL'] == 0]
coatedAbundance = abundance_binary[dataNew['PolyL'] == 1]
bareAbundance =  abundance_binary[dataNew['PolyL'] == 0]

nonPlList = list()
plList = list()
bothList = list()
structures = list()
for c in np.unique(overallRes3['Sample']):
    if('PL' in c):
        if(c !='BoxPL'):
            both=0
            tetOnly = 0
            tetCholOnly=0

            tet = overallRes3[overallRes3['Sample'] == c]
            if(c == 'I-10PL'):
                tetChol = overallRes3[overallRes3['Sample'] == c[:-2]]
            else:
                tetChol = overallRes3[overallRes3['Sample'] == c[:-3]]
            for z in np.unique(overallRes3['ID']):
                if((tet['ID'].str.contains(z).any()) and (tetChol['ID'].str.contains(z).any())):
                    tetVal = tet[tet['ID'] == z]['Abundance'].values[0]
                    tetCholVal = tetChol[tetChol['ID'] == z]['Abundance'].values[0]
                    if(tetVal == 1 and tetCholVal == 1):
                        both = both +1
                    if(tetVal == 0 and tetCholVal == 1):
                        tetCholOnly=tetCholOnly+1
                    if(tetVal ==1 and tetCholVal == 0):
                        tetOnly = tetOnly+1
                    else:
                        continue
            nonPlList.append(tetCholOnly)
            if(c == 'I-10PL'):
                structures.append(c[:-2])
            else:
                structures.append(c[:-3])
            plList.append(tetOnly)
            bothList.append(both)
            print(both/(both+tetOnly+tetCholOnly))
            #venn2(subsets = (tetOnly, both, tetCholOnly), set_labels = (c, c[:-2]))
            plt.show()
            #plt.savefig(c+'PLvsNonPL.png'bbox_inches='tight')
            #plt.clf()
dfStack=pd.DataFrame({'Coated':plList,'Both':bothList,'Bare':nonPlList}, index = structures)

fig = plt.figure(figsize=(8, 6))
plt.rc('font', family='serif')
plt.rc('font', family='serif')

dfStack.plot(kind='bar', stacked=True, color = sns.color_palette("colorblind", 3),legend=False)
heights = list()
    
plt.xlabel('Structures')
plt.ylabel('Number of Proteins in Corona')
plt.xticks(rotation=45)
#plt.legend(bbox_to_anchor=(8, 1.0), loc='upper left')
#plt.tight_layout()
# title of plot
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('Coated vs Uncoated Structures')
plt.savefig(outputDir + 'coatedvsUncoated.png',bbox_inches='tight')

both=0
tetOnly = 0
tetCholOnly=0

tet = overallRes3[overallRes3['Sample'] == 'tetrahedron']
tetChol = overallRes3[overallRes3['Sample'] == 'tetrahedron+cho']
for z in np.unique(overallRes3['ID']):
    if(z not in tet['ID'].values):
        tetVal = 0
    else:
        tetVal = tet[tet['ID'] == z]['Abundance'].values[0]
    if(z not in tetChol['ID'].values):
        tetCholVal = 0
    else:
        tetCholVal = tetChol[tetChol['ID'] == z]['Abundance'].values[0]
    if(tetVal == 1 and tetCholVal == 1):
        both = both +1
    if(tetVal == 0 and tetCholVal == 1):
        tetCholOnly=tetCholOnly+1
    if(tetVal ==1 and tetCholVal == 0):
        tetOnly = tetOnly+1
    else:
        continue
venn2(subsets = (tetCholOnly,both,tetOnly), set_labels = ('Tetrahedron w/Chol','Tetrahedron Only'))
plt.savefig(outputDir + 'tetvstetChol.png',bbox_inches='tight')

