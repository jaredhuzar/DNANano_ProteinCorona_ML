#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from RunModelAnalysisPresent import train_run_model
import pandas as pd
import os 
import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns

sns.set_style('white')
sns.set_context("paper", font_scale = 2)

datapresentPath=sys.argv[1]
outputDir = 'Results\\'

dataPresent = pd.read_csv(datapresentPath)
abundance_present = dataPresent['Present']
abundance_present = np.asarray(abundance_present)
overallRes3 = dataPresent[['Present', 'Sample', 'Protein']]
#dataPresent.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0'],inplace = True, axis =1)

dataNew=dataPresent
abundance_binary=abundance_present

overallRes3=overallRes3.rename(columns={'Present':'Abundance','Protein':'ID'})

polCoated = dataNew[dataNew['PolyL'] == 1]
bare = dataNew[dataNew['PolyL'] == 0]
coatedAbundance = abundance_binary[dataNew['PolyL'] == 1]
bareAbundance =  abundance_binary[dataNew['PolyL'] == 0]



def plotDataSubsets(dataPresent):
    abundance_binary = dataPresent['Present']
    abundance_binary = np.asarray(abundance_present)
    polCoated = dataNew[dataPresent['PolyL'] == 1]
    bare = dataNew[dataPresent['PolyL'] == 0]
    coatedAbundance = abundance_binary[dataPresent['PolyL'] == 1]
    bareAbundance =  abundance_binary[dataPresent['PolyL'] == 0]
    polCoated.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0','PolyL'],inplace = True, axis =1)
    bare.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0','PolyL'],inplace = True, axis =1)
    polCoated=polCoated.reset_index(drop=True)
    bare=bare.reset_index(drop=True)
    structMet=train_run_model(polCoated,coatedAbundance,'Coated')
    functMet=train_run_model(bare,bareAbundance,'Bare')
    
    metricComparison=pd.concat([structMet[0],functMet[0]])
    metricComparison.reset_index(drop=True, inplace = True)

    fig = plt.figure(figsize=(4, 3))
    plt.rc('font', family='serif')
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='large')
    plt.rc('ytick', labelsize='large')
    g = sns.catplot(
        data=metricComparison, kind="bar",
        x="Model", y="Scores", hue="Metric", palette='colorblind', alpha=.6, height=6,saturation = .7
    )
    #g.despine(left=True)
    g.set(ylim=(0, 1), title='Performance Coated and Uncoated')
    g.set_axis_labels("", "Score")
    plt.xticks(rotation=45)
    g.legend.set_title("")
    plt.savefig(outputDir + 'PerformanceCoatedUncoated.png',bbox_inches='tight')
    plt.clf()
    return metricComparison
    
    
metrics=plotDataSubsets(dataPresent)
metrics.to_csv(outputDir + 'CoatedUncoatedComparison.csv')

