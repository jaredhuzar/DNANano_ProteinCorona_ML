import scipy
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sns.set_style('white')
sns.set_context("paper", font_scale = 2)

proStatsPath=sys.argv[1]


outputDir = 'Results\\'

dataPresentPath=sys.argv[2]
dataPresent=pd.read_csv(dataPresentPath)
overallResPres = dataPresent[['Present', 'Sample', 'Protein']]
overallResPres=overallResPres.rename(columns={"Protein":'ID'})
overallResPres=overallResPres.rename(columns={"Present":'Abundance'})

proDatabase=pd.read_csv(proStatsPath)
proDatabase.index = proDatabase['Unnamed: 0']
proDatabase.drop('Unnamed: 0', axis =1, inplace = True)

pros = np.unique(overallResPres['ID'])
presPros = list()
absentPros = list()
for pro in pros:
    subRes=overallResPres[overallResPres['ID']==pro]
    if(subRes['Abundance'].sum() == 17):
        presPros.append(pro)
    elif(subRes['Abundance'].sum() == 0):
        absentPros.append(pro)
pd.DataFrame(presPros).to_csv(outputDir + 'UniversallyPresentProteins.csv')
pd.DataFrame(absentPros).to_csv(outputDir + 'UniversallyAbsentProteins.csv')
#allposPros=pd.read_csv(posProteinsPath)
#allnegPros=pd.read_csv(negProteinsPath)

posproData=proDatabase.loc[presPros]
negproData=proDatabase.loc[absentPros]

posproData['Dataset'] = 'Present'
negproData['Dataset'] = 'Absent'

comparison=pd.concat([posproData,negproData], ignore_index = True)
comparison['Dummy'] = 0
pVals = list()
zVals = list()
count = 0
for z in comparison.columns:
    if(z != 'Dataset'):
        if(z !='Dummy'):
            zVals.append(z)
            pVals.append(stats.ttest_ind(comparison[comparison['Dataset'] == 'Present'][z], comparison[comparison['Dataset'] == 'Absent'][z], equal_var = False)[1])
            if(stats.ttest_ind(comparison[comparison['Dataset'] == 'Present'][z], comparison[comparison['Dataset'] == 'Absent'][z], equal_var = False)[1] < .01):
                count = count+1
                
pvalDF=pd.DataFrame(pVals, index = zVals)
pVal=pvalDF.sort_values(by=0)
compSmall=comparison[comparison.columns.intersection(list(pVal[pVal.iloc[:,0] < .01].index))]
#compSmall.drop(['antigen binding'], axis = 1, inplace = True)
compSmall['Dataset'] = comparison['Dataset']
compSmall['Dummy'] = 0
valCols = list()
variables = list()
datasets = list()
for z in range(len(compSmall.columns)):
    if(compSmall.columns[z] != 'Dummy' and compSmall.columns[z] !='Dataset'):
        for a in range(len(compSmall)):
            valCol=compSmall.iloc[a,z]
            variable = compSmall.columns[z]
            datasetName = compSmall.loc[a,'Dataset']
            valCols.append(valCol)
            variables.append(variable)
            datasets.append(datasetName)
compPlot=pd.DataFrame({'Values':valCols,'Variables':variables,'Protein':datasets})

pVal=pvalDF.sort_values(by=0)

compSmall=comparison[comparison.columns.intersection(list(pVal[pVal.iloc[:,0] < .01].index))]
compSmall['Dummy'] = 0
compSmall['Dataset'] = comparison['Dataset']

valCols = list()
variables = list()
datasets = list()
for z in range(len(compSmall.columns)):
    if(compSmall.columns[z] != 'Dummy' and compSmall.columns[z] !='Dataset'):
        #if(compSmall.columns[z] in pVal[0:13].index):
        for a in range(len(compSmall)):
            valCol=compSmall.iloc[a,z]
            variable = compSmall.columns[z]
            datasetName = compSmall.loc[a,'Dataset']
            valCols.append(valCol)
            variables.append(variable)
            datasets.append(datasetName)
            
compPlot=pd.DataFrame({'Values':valCols,'Variables':variables,'Protein':datasets})

group = list()
remove = list()
for z in range(len(compPlot)):
    var = compPlot.loc[z,'Variables']
    if(len(var) == 5):
        group.append('ExposedAminoAcid')
    elif(len(var) == 1):
        group.append('SingleAminoAcid')
    elif(len(var) == 2):
        group.append('AminoSecondaryStructure')
    elif((var == 'flexMin') or (var == 'Net Charge') or (var == 'Aromaticity')or (var=='flexMed') or (var =='IsoelectricPoint')):
        group.append('Bulk')
    elif((var == 'Molecular Weight') or (var == 'LengthAminoAcids')):
        group.append('Size')
    else:
        remove.append(var)

compPlot = compPlot.loc[~compPlot['Variables'].isin(remove)]
compPlot.reset_index(drop=True, inplace = True)        
compPlot['Group'] = group

sns.violinplot(data=compPlot[compPlot['Group']=='SingleAminoAcid'], x='Variables', hue='Protein', y='Values', split=True,inner = 'quart', cut = 0, pallete = 'colorblind')
plt.xticks(rotation=45)
plt.savefig(outputDir + 'SingleAminoAcid_PresvsAbsent.png',bbox_inches='tight')
plt.clf()
#plt.show()

sns.violinplot(data=compPlot[compPlot['Group']=='AminoSecondaryStructure'], x='Variables', hue='Protein', y='Values', split=True,inner = 'quart', cut = 0, pallete = 'colorblind')
plt.xticks(rotation=45)
plt.savefig(outputDir + 'SecondaryAcid_PresvsAbsent.png',bbox_inches='tight')
plt.clf()
#plt.show()

sns.violinplot(data=compPlot[compPlot['Group']=='ExposedAminoAcid'], x='Variables', hue='Protein', y='Values', split=True,inner = 'quart', cut = 0, pallete = 'colorblind')
plt.xticks(rotation=45)
plt.savefig(outputDir + 'ExposedAminoAcid_PresvsAbsent.png',bbox_inches='tight')
#plt.clf()
#plt.show()

sns.violinplot(data=compPlot[compPlot['Group']=='Size'], x='Variables', hue='Protein', y='Values', split=True,inner = 'quart', cut = 0, pallete = 'colorblind')
plt.xticks(rotation=45)
plt.savefig(outputDir + 'size_PresvsAbsent.png',bbox_inches='tight')
plt.clf()
#plt.show()

sns.violinplot(data=compPlot[compPlot['Group']=='Bulk'], x='Variables', hue='Protein', y='Values', split=True,inner = 'quart', cut = 0, pallete = 'colorblind')
plt.xticks(rotation=45)
plt.savefig(outputDir + 'Bulk_PresvsAbsent.png',bbox_inches='tight')
plt.clf()

outputDir1 = 'Results\\allVariables\\'
vars1=np.unique(compPlot['Variables'])
for var in vars1:
    sns.violinplot(data=compPlot[compPlot['Variables']==var], x='Variables', hue='Protein', y='Values', split=True,inner = 'quart', cut = 0, pallete = 'colorblind')
    plt.xticks(rotation=45)
    plt.savefig(outputDir1 + var + '_PresvsAbsent.png',bbox_inches='tight')
    plt.clf()
    
    
sns.violinplot(data=compPlot[compPlot['Variables'].isin(['L','F','W'])], x='Variables', hue='Protein', y='Values', split=True,inner = 'quart', cut = 0, pallete = 'colorblind')
plt.xticks(rotation=45)
plt.title('Hydrophobic Amino Acids')
plt.savefig(outputDir + 'HydrophobicAA_PresvsAbsent.png',bbox_inches='tight')
plt.clf()
#plt.show()

sns.violinplot(data=compPlot[compPlot['Variables'].isin(['K'])], x='Variables', hue='Protein', y='Values', split=True,inner = 'quart', cut = 0, pallete = 'colorblind')
plt.xticks(rotation=45)
plt.title('Positively Charged Amino Acids')
plt.savefig(outputDir + 'PositiveAA_PresvsAbsent.png',bbox_inches='tight')
plt.clf()
#plt.show()

subPlot = compPlot[compPlot['Group']=='SingleAminoAcid']
subPlot = subPlot[subPlot['Variables'].isin(['C','P'])]
f, ax0 = plt.subplots()
ax1 = ax0.twinx()
var_order = ["C", "P"]
hue_order = ['Present', 'Absent']
i=0
for ax, var_name in zip([ax0, ax1], var_order):
    ax0.set_ylim(0,0.125)
    ax1.set_ylim(0,0.3)
    sns.violinplot(x=subPlot["Variables"], y=subPlot[subPlot['Variables'] == var_order[i]]['Values'], hue=subPlot["Protein"],ax=ax,split=True,inner = 'quart', cut = 0)
    i=i+1
plt.title('Special Amino Acids')
plt.savefig(outputDir + 'SingleAminoAcid_PresvsAbsent.png',bbox_inches='tight')

subPlot = compPlot[compPlot['Group']=='SingleAminoAcid']
subPlot = subPlot[subPlot['Variables'].isin(['C','K'])]
f, ax0 = plt.subplots()
ax1 = ax0.twinx()
var_order = ["C", "K"]
hue_order = ['Present', 'Absent']
i=0
for ax, var_name in zip([ax0, ax1], var_order):
    ax0.set_ylim(0,0.125)
    ax1.set_ylim(0,0.3)
    sns.violinplot(x=subPlot["Variables"], y=subPlot[subPlot['Variables'] == var_order[i]]['Values'], hue=subPlot["Protein"],ax=ax,split=True,inner = 'quart', cut = 0)
    i=i+1
plt.title('Plot in Paper')
plt.savefig(outputDir + 'PaperPlot_PresvsAbsent.png',bbox_inches='tight')

subPlot = compPlot[compPlot['Group']=='Size']
f, ax0 = plt.subplots()
ax1 = ax0.twinx()
var_order = ["LengthAminoAcids", "Molecular Weight"]
hue_order = ['Present', 'Absent']
i=0
for ax, var_name in zip([ax0, ax1], var_order):
    ax0.set_ylim(0,6000)
    ax1.set_ylim(0,600000)
    sns.violinplot(x=subPlot["Variables"], y=subPlot[subPlot['Variables'] == var_order[i]]['Values'], hue=subPlot["Protein"],ax=ax,split=True,inner = 'quart', cut = 0)
    i=i+1
plt.savefig(outputDir + 'Size_PresvsAbsent.png',bbox_inches='tight')

subPlot = compPlot[compPlot['Group']=='Bulk']
subPlot = subPlot[subPlot['Variables'].isin(['Net Charge','flexMin'])]
f, ax0 = plt.subplots()
ax1 = ax0.twinx()
var_order = ["Net Charge", "flexMin"]
hue_order = ['Present', 'Absent']
i=0
for ax, var_name in zip([ax0, ax1], var_order):
    ax0.set_ylim(0,0.5)
    ax1.set_ylim(.9,1)
    sns.violinplot(x=subPlot["Variables"], y=subPlot[subPlot['Variables'] == var_order[i]]['Values'], hue=subPlot["Protein"],ax=ax,split=True,inner = 'quart', cut = 0)
    i=i+1
plt.title('Plot in Paper')
plt.savefig(outputDir + 'Bulk_PresvsAbsent.png',bbox_inches='tight')

subPlot = compPlot[compPlot['Variables'].isin(['HC','HR'])]
f, ax0 = plt.subplots()
ax1 = ax0.twinx()
var_order = ["HC", "HR"]
hue_order = ['Present', 'Absent']
i=0
for ax, var_name in zip([ax0, ax1], var_order):
    ax0.set_ylim(0,0.05)
    ax1.set_ylim(0,0.12)
    sns.violinplot(x=subPlot["Variables"], y=subPlot[subPlot['Variables'] == var_order[i]]['Values'], hue=subPlot["Protein"],ax=ax,split=True,inner = 'quart', cut = 0)
    i=i+1
plt.title('Different Seconday Structure')
plt.savefig(outputDir + 'secondStructure_PresvsAbsent.png',bbox_inches='tight')

