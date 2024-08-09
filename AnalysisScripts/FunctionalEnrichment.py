import pandas as pd
import numpy as np
import sys
import scipy
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style('white')
sns.set_context("paper", font_scale = 2)

dataPresentPath= sys.argv[1]
funcListPath = sys.argv[2]
outputDir = 'Results\\'

dataPresent = pd.read_csv(dataPresentPath)
abundance_present = dataPresent['Present']
abundance_present = np.asarray(abundance_present)
dataPresent.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0'],inplace = True, axis =1)

functions = pd.read_csv(funcListPath)
functions=functions.transpose()
functions.drop('Unnamed: 0', inplace = True, axis = 0)
functions.columns = functions.iloc[0,:]

def plotFunctionalEnrichment(dataDF,abundanceArray, functions):
    listVals = list()
    enrichVals = list()
    pVals = list()
    functs = list()
    pValDe = list()
    for z in functions:
        if(sum(dataDF[z]) > 390): # 23 or more proteins with the function
            funcArray=dataDF[z]
            functs.append(z)
            totSamples = len(abundanceArray)
            totfunctionalSamples = len(abundanceArray[funcArray == 1])
            totNonfunctionalSamples = len(abundanceArray[funcArray == 0])
            funcandAbundant = pd.DataFrame(abundanceArray)[funcArray == 1].value_counts()[1]
            totAbundant = pd.DataFrame(abundanceArray).value_counts()[1]
            defuncandAbundant = pd.DataFrame(abundanceArray)[funcArray == 0].value_counts()[1]
            pValDe.append(1-stats.hypergeom.cdf(defuncandAbundant-1, totSamples, totAbundant, totNonfunctionalSamples))
            enrichVals.append(funcandAbundant*totSamples/(totfunctionalSamples*totAbundant))
            pVals.append(1 - stats.hypergeom.cdf(funcandAbundant-1, totSamples, totAbundant, totfunctionalSamples))
            resArray = abundanceArray
            listVals.append(pd.DataFrame(abundanceArray)[funcArray == 1].value_counts()[1]/(pd.DataFrame(abundanceArray)[funcArray == 1].value_counts()[1] + pd.DataFrame(abundanceArray)[funcArray == 1].value_counts()[0]))


    pValSig = list()
    structSig = list()
    enrichSig = list()
    for z in range(len(pVals)):
        if(pVals[z]<0.05):
            pValSig.append(pVals[z])
            structSig.append(functs[z])
            enrichSig.append(enrichVals[z])

    for z in range(len(pValDe)):
        if(pValDe[z]<0.05):
            pValSig.append(pValDe[z])
            structSig.append(functs[z])
            enrichSig.append(enrichVals[z])

    functionalVals=pd.DataFrame({'Function':structSig,'Enrichment Value':enrichSig})
    funcVals2=functionalVals.transpose()
    funcVals2.columns = funcVals2.iloc[0]
    funcVals2.drop(['Function'], axis =0, inplace = True)
    funcVals2.sort_values(by = ['Enrichment Value'], axis = 1)

    ax = sns.barplot(data = funcVals2.sort_values(by = ['Enrichment Value'], axis = 1,ascending = False),orient='h')
    plt.axvline(1,linewidth=4, color='r')
    plt.title('Functional Protein Enrichment')
    #plt.yaxis('Enrichment Value')
    ax.set_xlabel('Enrichment Value')
    #plt.plot()
    plt.savefig(outputDir + 'FunctionalProteinEnrichment.png',bbox_inches='tight')
    return(funcVals2.sort_values(by = ['Enrichment Value'], axis = 1,ascending = False))

functionalTable=plotFunctionalEnrichment(dataPresent,abundance_present, functions)
functionalTable.to_csv(outputDir + 'FunctionalEnrichment.csv')

