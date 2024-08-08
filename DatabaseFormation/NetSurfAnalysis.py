#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
import os
import numpy as np
import warnings
import sys
warnings.filterwarnings('ignore')

inputDir = sys.argv[1] #"Data\\netSurf\\"
outputDir = 'AnalyzedDataTest\\'

netSurf1=pd.read_csv(inputDir + '10thNetSurf59.csv')
for fil in os.listdir(inputDir):
    if(fil != '10thNetSurf59.csv'):
        netSurf2=pd.read_csv(inputDir + fil)
        netSurf1 = netSurf1.append(netSurf2)
netSurf2=pd.read_csv(inputDir + '11thNetSurf.csv')
netSurf1 = netSurf1.append(netSurf2)
netSurf1.reset_index(inplace=True, drop = True)
netSurf = netSurf1
netSurf=netSurf.drop_duplicates(subset=['id', 'n'])


netSurf.to_csv(outputDir + 'AllNetSurf.csv')

disMean = list()
psiMean = list()
phiMean = list()
totList = list()
totListaa = list()
pros = list()
for x in range(len(np.unique(netSurf['id']))):
    protein1=np.unique(netSurf['id'])[x]
    pros.append(protein1)
    subSurf =netSurf[netSurf['id'] == protein1]
    disMean.append(np.mean(subSurf['disorder']))
    psiMean.append(np.mean(subSurf['psi']))
    phiMean.append(np.mean(subSurf['phi']))
    perProteinList = list()
    for y in np.unique(netSurf['q8']):
        if(y in subSurf['q8'].value_counts()):
            perProteinList.append(subSurf['q8'].value_counts()[y]/sum(subSurf['q8'].value_counts()))
        else:
            perProteinList.append(0)
    totList.append(perProteinList)

    perProteinListaa = list()
    for y in np.unique(netSurf['seq']):
        if(y in subSurf['seq'].value_counts()):
            perProteinListaa.append(subSurf['seq'].value_counts()[y]/sum(subSurf['seq'].value_counts()))
        else:
            perProteinListaa.append(0)
    totListaa.append(perProteinListaa)
    
totListaaExpo = list()
for x in range(len(np.unique(netSurf['id']))):
    protein1=np.unique(netSurf['id'])[x]
    subSurf =netSurf[netSurf['id'] == protein1]
    perProteinListaaExpo = list()
    for y in np.unique(netSurf['seq']):
        if(y in subSurf['seq'].value_counts()):
            perProteinListaaExpo.append(sum(subSurf[subSurf['seq'] == y]['asa'])/sum(subSurf['asa']))
        else:
            perProteinListaaExpo.append(0)
    totListaaExpo.append(perProteinListaaExpo)
    
aminoAcids=np.unique(netSurf['seq'])
i=0
for x in aminoAcids:
    aa1=x + 'Expo'
    aminoAcids[i] = aa1
    i=i+1

totaa=pd.DataFrame(totListaa)
totaa.index = pros
totaa.columns = np.unique(netSurf['seq'])

tot=pd.DataFrame(totList)
q8list = list()
for a in np.unique(netSurf['q8']):
    q8list.append(a + 'q8')
tot.index = pros
tot.columns = q8list


totListaaExpo=pd.DataFrame(totListaaExpo)
totListaaExpo.index =  pros
totListaaExpo.columns = aminoAcids

#allProteinData=pd.concat([allProteinData, totaa], axis = 1)
#allProteinData=pd.concat([allProteinData, tot], axis = 1)

#allProteinData['DisorderMean'] = disMean
#allProteinData['psiMean'] = psiMean
#allProteinData['phiMean'] = phiMean

moreNetSurf=pd.DataFrame({'disMean':disMean, 'phiMean':phiMean,'psiMean':psiMean}, index = pros)

combos = list()
cho =pd.DataFrame(netSurf[['q3','seq']].value_counts()).index
for s in range(len(cho)):
    combos.append(cho[s][0]+cho[s][1])

combos=pd.DataFrame(combos, columns = ['Combos'])
combos.index = combos['Combos']

for z in range(len(np.unique(netSurf['id']))):
    protein1=np.unique(netSurf['id'])[z]
    subSurf =netSurf[netSurf['id'] == protein1]
    valCount=pd.DataFrame(subSurf[['q3','seq']].value_counts())
    totSum=sum(valCount[0])
    combos[protein1] = np.zeros(len(combos))
    for t in range(len(valCount)):
        valKey=valCount.index[t][0]+valCount.index[t][1]
        val = valCount[0][t]
        combos.loc[valKey,protein1] = val/totSum
        
combos=combos.drop(['Combos'], axis = 1)
combos=combos.transpose()
combos = combos.reset_index(drop = True)
combos.index = pros

combos.to_csv(outputDir + 'combos.csv')
moreNetSurf.to_csv(outputDir + 'moreNetSurf.csv')
totListaaExpo.to_csv(outputDir + 'exposedAminoAcid.csv')
totaa.to_csv(outputDir+'aminoAcid.csv')
tot.to_csv(outputDir + 'secondaryAminoAcid.csv')

