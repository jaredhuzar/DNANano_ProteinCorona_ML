import pandas as pd
import numpy as np
import sys

outputDir = 'AnalyzedDataTest\\'


origami=pd.read_excel(sys.argv[1])
origami.index = origami['Unnamed: 0']
origami.drop(['Unnamed: 0'], axis=1,inplace = True)
origami=origami.transpose()

newPresent = pd.read_csv(sys.argv[2]) #'TranslatedPresent8_6.csv'
newPresent.drop('Unnamed: 0',axis =1 ,inplace = True)

mergedPros = pd.read_csv(sys.argv[3],index_col=0)

newPresent=newPresent.loc[newPresent['ID'].isin(set(newPresent['ID']).intersection(set(mergedPros.index)))]
newPresent.reset_index(drop=True, inplace=True)
newPresent = newPresent[newPresent.Sample != 'MB']
newPresent.reset_index(drop=True, inplace=True)

origami.index = np.unique(newPresent['Sample']) 

rows = list()
for z in range(len(newPresent)):
    prot=newPresent.loc[z,'ID']
    nano=newPresent.loc[z,'Sample']
    if(nano !='MB'):
        a = list(origami.loc[nano])
        b = list(mergedPros.loc[prot])
        b.extend(a)
    rows.append(b)
data=pd.DataFrame(rows)
#data.columns = mergedPros.columns
cols1=np.asarray([np.asarray(mergedPros.columns),np.asarray(origami.columns)], dtype='object')
data.columns = np.concatenate(cols1).ravel()
data['Present'] = newPresent['Abundance']
data['Sample'] = newPresent['Sample']
data['Protein'] = newPresent['ID']

data.to_csv(outputDir + 'overallResPresent8_6.csv')

rows = list()
for z in range(len(newPresent)):
    prot=newPresent.loc[z,'ID']
    nano=newPresent.loc[z,'Sample']
    protOld = newPresent.loc[z,'Protein']
    if(nano !='MB'):
        a = list(origami.loc[nano])
        b = list(mergedPros.loc[prot])
        cVal = allEnrich.loc[(allEnrich['Sample'] == nano) & (allEnrich['Protein'] == protOld), 'Abundance'].values[0]
        b.extend(a)
        b.append(cVal)
    rows.append(b)
dataTotEnrich=pd.DataFrame(rows)
enrichColumns = list(np.concatenate(cols1).ravel())
enrichColumns.append('Enrichment')
dataTotEnrich.columns = enrichColumns

dataTotEnrich['Sample'] = data['Sample']
dataTotEnrich['Protein'] = data['Protein']

dataTotEnrich.to_csv(outputDir + 'overallResEnriched8_6.csv')

