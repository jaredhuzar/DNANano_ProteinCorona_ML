import pandas as pd
import numpy as np
import re
import requests as r
import sys

presentProteinsPath=sys.argv[1] #"PresentProteins8_6.csv"
outputDir='AnalyzedDataTest\\'

idTransfers = {'Ubiquitin-60S ribosomal protein L40 ':'P62987','Ubiquitin-40S ribosomal protein S27a ':'P62979',
               'RPL22 ':'P35268','Protein KIAA0100 ':'Q14667','ITIH4 protein ':'Q14624','Heat shock 70kDa protein 9B (Mortalin-2) ':'H0Y8K0',
               'HCG40295 ':'F6WE04','HCG2029799 ':'Q5TEH5','GTP-binding protein SAR1b ':'Q9Y6B6','GTP-binding protein SAR1a ':'Q9NR31',
               'Complement subcomponent C1r ':'P00736', 'Alpha-2-macroglobulin ':'P01023', 'Profilin ':'P07737','Sushi domain-containing protein':'Q6UWL2'
              }

tots=pd.read_csv(presentProteinsPath)
tots.drop(['Unnamed: 0'],inplace = True, axis = 1)

namesUnique = np.unique(tots['Protein'])
proteinIDs = namesUnique
Accessions = list()
ids = list()
noIDs = list()
proteinIDsNew = list()
for pro9 in proteinIDs:
    url = 'https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+AND(protein_name:"' + pro9 + '")AND(reviewed:true)&fields=id&format=json&compressed=false&size=1'
    response = r.post(url)
    Data=''.join(response.text)
    
    if(pro9 in list(idTransfers.keys())):
        Accessions.append(idTransfers[pro9])
        proteinIDsNew.append(pro9)
        ids.append('Manual')
    
    else:
    
        test_str = Data
        test_str=test_str.replace('"primaryAccession":"',"*")
        test_str=test_str.replace('","',"*")
        re=test_str.split("*")
        if(len(re)>1):
            Accession=re[1]
            Accessions.append(Accession)
            proteinIDsNew.append(pro9)
        else:
            noIDs.append(pro9)

        test_str = Data
        test_str=test_str.replace('uniProtkbId":"',"*")
        test_str=test_str.replace('"}',"*")
        re=test_str.split("*")
        if(len(re) > 1):
            uniprotID=re[1]
            ids.append(uniprotID)

Translate=pd.DataFrame({'Old':proteinIDsNew,'New':Accessions})

tots['ID'] = np.nan
for z in range(len(Translate)):
    tots.loc[tots['Protein'] == Translate['Old'][z], 'ID'] = Translate['New'][z]
tots=tots.dropna(subset = ['ID'], axis =0)

allPresent = tots
for z in np.unique(allPresent['Sample']):
    subPresent = allPresent[allPresent['Sample'] == z]
    subsubPres=subPresent.loc[subPresent.duplicated('ID', keep = False)]
    #print(subsubPres.index)
    dropIndices = list()
    keepIndices = list()
    for a in np.unique(subsubPres['ID']):
        perIDpres=subsubPres[subsubPres['ID'] ==a]
        if(len(perIDpres) == 3):
            if(sum(perIDpres['Abundance']) == 0):
                dropIndices.append(perIDpres.index[1])
                dropIndices.append(perIDpres.index[0])
                keepIndices.append(perIDpres.index[2])
            else:
                if(perIDpres.iloc[0,1] == 1):
                    dropIndices.append(perIDpres.index[1])
                    dropIndices.append(perIDpres.index[2])
                    keepIndices.append(perIDpres.index[0])
                elif(perIDpres.iloc[1,1] == 1):
                    dropIndices.append(perIDpres.index[0])
                    dropIndices.append(perIDpres.index[2])
                    keepIndices.append(perIDpres.index[1])
                else:
                    dropIndices.append(perIDpres.index[0])
                    dropIndices.append(perIDpres.index[1])
                    keepIndices.append(perIDpres.index[2])
        if(perIDpres.iloc[0,1] == 1):
            dropIndices.append(perIDpres.index[1])
            keepIndices.append(perIDpres.index[0])
        elif(perIDpres.iloc[1,1] == 1):
            dropIndices.append(perIDpres.index[0])
            keepIndices.append(perIDpres.index[1])
        else:
            dropIndices.append(perIDpres.index[0])
            keepIndices.append(perIDpres.index[1])
        #else:
        #    dropIndices.append(perIDpres.index[0])
        #    dropIndices.append(perIDpres.index[1])
    if(z=='BoxPL'):
        allPresNew=subPresent.drop(dropIndices, axis =0)
    else:
        addArray=subPresent.drop(dropIndices, axis =0)
        allPresNew=allPresNew.append(addArray)
        
allPresNew=allPresNew.reset_index(drop = True)

allPresNew.to_csv(outputDir + 'TranslatedPresent8_6.csv')

