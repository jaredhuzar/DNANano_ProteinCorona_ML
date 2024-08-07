#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import upsetplot
from upsetplot import plot

dataPresentPath=sys.argv[1]
dataPresent=pd.read_csv(dataPresentPath)
abundance_present = dataPresent['Present']
abundance_present = np.asarray(abundance_present)
overallResPres = data[['Present', 'Sample', 'Protein']]
dataPresent.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0'],inplace = True, axis =1)


def commonProteins(listPros,overallRes3, label):
    allPros = np.unique(overallRes3['ID'])
    setNames = list()
    listSets = list()
    #for c in np.unique(overallRes3['Sample']):
    #for c in ['tube','tile','tensegrity','I-12']:
    #for c in ['tube','tile','tensegrity','I-12']:
    for c in listPros:
        subRes = overallRes3[overallRes3['Sample'] == c]
        subRes = subRes[subRes['Abundance'] == 1]
        listSets.append(set(subRes['ID']))
        setNames.append(c)

    setNums=pd.DataFrame([[e in setA for setA in listSets] for e in allPros], columns = setNames)

    yVals = list()
    xVals = list()
    resVals = list()
    for y in range(len(listSets)):
        for x in range(len(listSets)):
            if(x == y):
                resVals.append(1)
            else:
                common=listSets[x].intersection(listSets[y])
                totVals = listSets[x].union(listSets[y])
                resVals.append(len(common)/len(totVals))
            xVals.append(setNames[x])
        yVals.append(setNames[y])
        
    percentSim=pd.DataFrame(np.asarray(resVals).reshape((len(listPros),len(listPros))), index =yVals, columns = np.unique(xVals))
    percentSim = percentSim*100
    figsize=(6, 6)
    sns.set_context("paper", font_scale = 1.8)
    plt.rc('font', family='serif')
    cg = sns.clustermap(percentSim.loc[listPros,listPros], annot=True,cmap="Blues",fmt ='.0f',cbar_pos=(0, .4, .08, .4),annot_kws={"size": 8})
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xmajorticklabels(), fontsize = 16)
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize = 16)
    #plt.savefig('StructureCoronaSimilarityHeatmap.png',bbox_inches='tight')

    allPros = np.unique(overallRes3['ID'])
    df = pd.DataFrame([[e in setA for setA in listSets] for e in allPros], columns = setNames)
    df=df.loc[:,listPros]
    df_up = df.value_counts()
    fig1=upsetplot.UpSet(df_up,sort_by ='cardinality', min_degree = 1)
    #plot(fig1)
    fig1.plot()
    plt.savfig(outputDir + label +'CommonProteins.png')
    plt.clf()
    return pd.DataFrame(df_up)


overallResPres=overallResPres.rename(columns={"Protein":'ID'})
overallResPres=overallResPres.rename(columns={"Present":'Abundance'})
dfAllCommon=commonProteins(np.unique(overallResPres['Sample']),overallResPres, 'AllSamples')

overallPL=overallResPres[overallResPres.Sample.str.contains("PL")]
dfAllPL=commonProteins(np.unique(overallPL['Sample']),overallPL, 'PLSamples')

overallBare=overallResPres[~ overallResPres.Sample.str.contains("PL")]
overallBare=overallBare[~ overallBare.Sample.str.contains("cho")]
dfAllBare=commonProteins(np.unique(overallBare['Sample']),overallBare, 'BareSamples')

overall2D=overallResPres.loc[overallResPres['Sample'].isin(['tile','I-10','I-12','I-14','I-16'])]
df2D=commonProteins(np.unique(overall2D['Sample']),overall2D, '2DSamples')

overallDimensions=overallResPres.loc[overallResPres['Sample'].isin(['tube','tensegrity','tile','I-12'])]
dfDimensions=commonProteins(np.unique(overallDimensions['Sample']),overallDimensions, '2d_3dSamples')

