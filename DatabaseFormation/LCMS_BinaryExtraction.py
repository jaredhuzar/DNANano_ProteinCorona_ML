#!/usr/bin/env python
# coding: utf-8

# # Getting LCMS Data

# In[ ]:


import sys
import pandas as pd
import numpy as np
import math

lcmsDataPath = sys.argv[1]
outputDir = 'AnalyzedDataTest\\'

df = pd.read_excel(lcmsDataPath, None)

df1 = pd.read_excel(lcmsDataPath,sheet_name = list(df.keys())[0])
df1.columns = df1.iloc[0,:]
df1.drop([0],axis = 0, inplace = True)
df1.index = df1['Protein']
df1.drop(['Protein'], axis = 1, inplace = True)

samples = list()
allPres = list()
for z in df1.columns:
    proteinPres = list()
    for x in df1[z]:
        if(math.isnan(x)):
            proteinPres.append(0)
        else:
            proteinPres.append(1)
    samples.append(z)
    allPres.append(proteinPres)
    
allPresdf=pd.DataFrame(allPres).transpose()
allPresdf.columns = samples
allPresdf.index = df1.index #Any amount of protein in the corona

df1 = pd.read_excel(lcmsDataPath,sheet_name = list(df.keys())[0])
df1.columns = df1.iloc[0,:]
df1.drop([0],axis = 0, inplace = True)
df1.index = df1['Protein']
df1.drop(['Protein'], axis = 1, inplace = True)

samples = list()
allPres = list()
for z in df1.columns:
    proteinPres = list()
    for x in df1[z]:
        if(math.isnan(x)):
            proteinPres.append(0)
        else:
            if(x > 0):
                proteinPres.append(1)
            else:
                proteinPres.append(0)
    samples.append(z)
    allPres.append(proteinPres)
    
allPresdfEnrich=pd.DataFrame(allPres).transpose()
allPresdfEnrich.columns = samples
allPresdfEnrich.index = df1.index #Only more in the corona than in serum

df2 = pd.read_excel(lcmsDataPath,sheet_name = list(df.keys())[1])
samples = list()
allPros=pd.DataFrame({'Protein':allPresdf.index})
for z in df2.columns:
    allPros[z] = 0
    listPros=df2[z].dropna()
    print(set(listPros).difference(set(allPresdf.index)))
    allPros.loc[allPros['Protein'].isin(listPros),z] = 1
    samples.append(z)
allPros.index  = allPros['Protein']
allPros.drop(['Protein'], axis = 1, inplace = True) #Not found in serum, only on MB and in NPs

df3 = pd.read_excel(lcmsDataPath,sheet_name = list(df.keys())[2])
samples = list()
allProsCorona=pd.DataFrame({'Protein':allPresdf.index})
for z in df3.columns:
    allProsCorona[z] = 0
    listPros=df3[z].dropna()
    print(set(listPros).difference(set(allPresdf.index)))
    allProsCorona.loc[allProsCorona['Protein'].isin(listPros),z] = 1
    samples.append(z)
allProsCorona.index  =allProsCorona[ 'Protein']
allProsCorona.drop(['Protein'], axis = 1, inplace = True) #Not on MB or in serum, only nanostructure corona

mbandCorona=allProsCorona + allPros #Only found in nanostructures (no MB or serum) and not found in serum (just MB and NPs)
allPresdf.columns=allProsCorona.columns 
coronaandfc= allPresdf + allProsCorona #Only found in nanostructure (not MB or serum), and any amount of protein in corona (enriched or depleted)

coronaPresent = allPros + allPresdf + allProsCorona
coronaPresent.index = allPresdf.index
coronaPresent=coronaPresent.replace(2,1)
coronaPresent=coronaPresent.replace(3,1)

samples = list()
proteins = list()
abundanceVals = list()
for z in coronaPresent.columns:
    for x in range(len(coronaPresent[z])):
        proteins.append(coronaPresent.index[x])
        samples.append(z)
        abundanceVals.append(coronaPresent.loc[coronaPresent.index[x],z])
allPresent = pd.DataFrame({'Protein':proteins, 'Abundance':abundanceVals, 'Sample':samples})
allPresent.to_csv(outputDir + 'PresentProteins8_4.csv')

allPresdfEnrich.columns=allProsCorona.columns
coronaEnrich = allPresdfEnrich+allProsCorona+allPros
coronaEnrich.index = allPresdfEnrich.index
coronaEnrich=coronaEnrich.replace(2,1)
coronaEnrich=coronaEnrich.replace(3,1)
samples = list()
proteins = list()
abundanceVals = list()
for z in coronaEnrich.columns:
    for x in range(len(coronaEnrich[z])):
        proteins.append(coronaEnrich.index[x])
        samples.append(z)
        abundanceVals.append(coronaEnrich.loc[coronaEnrich.index[x],z])
allEnrich = pd.DataFrame({'Protein':proteins, 'Abundance':abundanceVals, 'Sample':samples})
allEnrich.to_csv(outputDir + 'EnrichedProteins8_4.csv')

