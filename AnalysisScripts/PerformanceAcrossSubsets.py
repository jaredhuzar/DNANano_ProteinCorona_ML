from RunModelAnalysisPresent import train_run_model
import pandas as pd
import os 
import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns

datafilePath=sys.argv[1]
origamidataPath=sys.argv[2]
funcListPath=sys.argv[3]

inputDir = 'AnalyzedData\\'
outputDir = 'Results\\'

dataPresent = pd.read_csv(datafilePath)
abundance_present = dataPresent['Present']
abundance_present = np.asarray(abundance_present)
dataPresent.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0'],inplace = True, axis =1)

origami=pd.read_excel(origamidataPath)
origami.index = origami['Unnamed: 0']
origami.drop(['Unnamed: 0'], axis=1,inplace = True)
origami=origami.transpose()

functions = pd.read_csv(funcListPath)
functions=functions.transpose()
functions.drop('Unnamed: 0', inplace = True, axis = 0)
functions.columns = functions.iloc[0,:]


def subset_data(inputArr,param):
    if(param == 'Origami'):
        outArr = inputArr[origami.columns]
    elif(param == 'Structure'):
        a2 = list(origami.columns)
        a2.append('Interactivity')
        a3 = list(functions.columns)
        a4=a2+a3
        outArr = inputArr.drop(np.array(a4),axis=1)
        return(outArr)
    elif(param == 'Functions'):
        a2 = list(functions.columns)
        a2.append('Interactivity')
        outArr=inputArr[a2]
    elif(param == 'All'):
        outArr = inputArr
    return(outArr)

def plotDataSubsets(dataPresent):
    Xstruct=subset_data(dataPresent,'Structure')
    Xfunct=subset_data(dataPresent,'Functions')
    Xorigami=subset_data(dataPresent,'Origami')

    structMet=train_run_model(Xstruct,abundance_present,'Structure')
    functMet=train_run_model(Xfunct,abundance_present,'Functional')
    origMet=train_run_model(Xorigami,abundance_present,'Origami')
    statsPres=train_run_model(dataPresent,abundance_present,'All')
    
    metricComparison=pd.concat([statsPres[0],origMet[0],functMet[0], structMet[0]])
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
    g.set(ylim=(0, 1), title='Performance Across Data Subsets')
    g.set_axis_labels("", "Score")
    plt.xticks(rotation=45)
    g.legend.set_title("")
    plt.savefig(outputDir + 'PerformanceSubsets.png',bbox_inches='tight')
    plt.clf()
    
    
plotDataSubsets(dataPresent)

