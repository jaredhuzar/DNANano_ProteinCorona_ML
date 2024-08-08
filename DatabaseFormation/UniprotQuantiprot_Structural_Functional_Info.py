import re
import pandas as pd
import numpy as np
import requests as r
import warnings
warnings.filterwarnings("ignore")

import quantiprot
from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.metrics.aaindex import get_aaindex_file
from quantiprot.metrics.basic import average
from quantiprot.metrics.aaindex import get_aa2charge, get_aa2hydropathy
from quantiprot.metrics.basic import average, average_absolute
import Bio
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import sys
from io import StringIO
from Bio import SeqIO

outputDir = 'AnalyzedDataTest\\'
allPresentPath = sys.argv[1] #'TranslatedPresent8_4.csv'

allPresNew = pd.read_csv(allPresentPath)
proteinIDs = np.unique(allPresNew['ID'])

functions = list()
num_interactions = list()

proteinSeqs = list()
gravyVals = list()
monotropicBools = list()
secondarySturctures = list()
molarExtinctions = list()
mws = list()
instabilities = list()
lens = list()
flexes = list()
isoelectrics = list()
neutralCharges = list()
aromaticities = list()


for ID in proteinIDs:
    url = "http://www.uniprot.org/uniprot/"
    proteinURL = url + ID + '.txt'
    response = r.post(proteinURL,verify=False)
    Data=''.join(response.text)
    string1 = ID + '; (.*); NbExp'
    len(re.findall(string1, Data))
    num_interactions.append(len(re.findall(string1, Data)))
    functions.append(re.findall('; F:(.*); ', Data))

    proteinURL = url + ID + '.fasta'
    response = r.post(proteinURL,verify=False)
    Data=''.join(response.text)
    proteinSeq=StringIO(Data)
    proteinSeq=SeqIO.read(proteinSeq,'fasta')
    proteinSeqs.append(proteinSeq.seq)
    if(ID == 'P22352'): #Removed U
        analysed_seq  = ProteinAnalysis('MARLLQASCLLSLLLAGFVSQSRGQEKSKMDCHGGISGTIYEYGALTIDGEEYIPFKQYAGKYVLFVNVASYUGLTGQYIELNALQEELAPFGLVILGFPCNQFGKQEPGENSEILPTLKYVRPGGGFVPNFQLFEKGDVNGEKEQKFYTFLKNSCPPTSELLGTSDRLFWEPMKVHDIRWNFEKFLVGPDGIPIMRWHHRTTVSNVKMDILSYMRRQAALGVKRK')
        analysed_seq2  = ProteinAnalysis('MARLLQASCLLSLLLAGFVSQSRGQEKSKMDCHGGISGTIYEYGALTIDGEEYIPFKQYAGKYVLFVNVASYGLTGQYIELNALQEELAPFGLVILGFPCNQFGKQEPGENSEILPTLKYVRPGGGFVPNFQLFEKGDVNGEKEQKFYTFLKNSCPPTSELLGTSDRLFWEPMKVHDIRWNFEKFLVGPDGIPIMRWHHRTTVSNVKMDILSYMRRQAALGVKRK')
    elif(ID == 'P49908'):
        analysed_seq  = ProteinAnalysis('MWRSLGLALALCLLPSGGTESQDQSSLCKQPPAWSIRDQDPMLNSNGSVTVVALLQASUYLCILQASKLEDLRVKLKKEGYSNISYIVVNHQGISSRLKYTHLKNKVSEHIPVYQQEENQTDVWTLLNGSKDDFLIYDRCGRLVYHLGLPFSFLTFPYVEEAIKIAYCEKKCGNCSLTTLKDEDFCKRVSLATVDKTVETPSPHYHHEHHHNHGHQHLGSSELSENQQPGAPNAPTHPAPPGLHHHHKHKGQHRQGHPENRDMPASEDLQDLQKKLCRKRCINQLLCKLPTDSELAPRSUCCHCRHLIFEKTGSAITUQCKENLPSLCSUQGLRAEENITESCQURLPPAAUQISQQLIPTEASASURUKNQAKKUEUPSN')
        analysed_seq2  = ProteinAnalysis('MWRSLGLALALCLLPSGGTESQDQSSLCKQPPAWSIRDQDPMLNSNGSVTVVALLQASYLCILQASKLEDLRVKLKKEGYSNISYIVVNHQGISSRLKYTHLKNKVSEHIPVYQQEENQTDVWTLLNGSKDDFLIYDRCGRLVYHLGLPFSFLTFPYVEEAIKIAYCEKKCGNCSLTTLKDEDFCKRVSLATVDKTVETPSPHYHHEHHHNHGHQHLGSSELSENQQPGAPNAPTHPAPPGLHHHHKHKGQHRQGHPENRDMPASEDLQDLQKKLCRKRCINQLLCKLPTDSELAPRSCCHCRHLIFEKTGSAITQCKENLPSLCSQGLRAEENITESCQRLPPAAQISQQLIPTEASASRKNQAKKEPSN')
    else:
        analysed_seq  = ProteinAnalysis(proteinSeq.seq)
        analysed_seq2 = ProteinAnalysis(proteinSeq.seq)

    monotropicBools.append(analysed_seq.monoisotopic)
    secondarySturctures.append(analysed_seq.secondary_structure_fraction())
    gravyVals.append(analysed_seq2.gravy())
    molarExtinctions.append(analysed_seq.molar_extinction_coefficient())
    mws.append(analysed_seq.molecular_weight())
    instabilities.append(analysed_seq2.instability_index())
    lens.append(analysed_seq.length)
    flexes.append(analysed_seq2.flexibility())
    isoelectrics.append(analysed_seq.isoelectric_point())
    neutralCharges.append(analysed_seq.charge_at_pH(7))
    aromaticities.append(analysed_seq.aromaticity())
    
a1=np.asarray(functions, dtype=object)
funcList=sum(functions, [])
functions1=pd.DataFrame(index = proteinIDs,columns = np.unique(funcList))
#functions1 = functions1.drop('Q9Y6R7')

for x in range(len(a1)):
    pro1=a1[x]
    if(len(pro1) > 0):
        for y in range(len(pro1)):
            func=pro1[y]
            functions1.loc[functions1.index[x],func] = 1
functions1=functions1.fillna(0)

functions1.head()

flexMed = list()
flexMax = list()
flexMin = list()
for z in range(len(flexes)):
    flex=flexes[z]
    flexMed.append(np.median(flex))
    flexMax.append(max(flex))
    flexMin.append(min(flex))
goFlex=pd.DataFrame({'flexMed':flexMed, 'flexMin':flexMin,'flexMax':flexMax},index = proteinIDs)

Interactions=pd.DataFrame({'Interactivity':num_interactions}, index = proteinIDs)
functions=list(functions1.columns)

pd.DataFrame(functions).to_csv(outputDir + 'FunctionList8_6.csv')
functions1.to_csv(outputDir + 'functionalInfo.csv')
Interactions.to_csv(outputDir + 'interactivityData.csv')
goFlex.to_csv(outputDir + 'flexibility.csv')


ofile = open(outputDir + "allProteins8_6.txt", "w")

for i in range(len(proteinSeqs)):
    ofile.write(">" + proteinIDs[i] + "\n" +str(ProteinAnalysis(proteinSeqs[i]).sequence)+ "\n")

ofile.close()


uversky_fs = FeatureSet("uversky")
net_abs_charge = Feature(get_aa2charge(default=0)).then(average_absolute)
mean_hydropathy = Feature(get_aa2hydropathy(default=0)).then(average)
uversky_fs.add(mean_hydropathy, name="mean_hydropathy")
uversky_fs.add(net_abs_charge, name="net_abs_charge")

allProteins = load_fasta_file(outputDir + 'allProteins8_6.txt')
allProteins_seq = uversky_fs(allProteins)
hydropathy=allProteins_seq.columns(feature="mean_hydropathy")[0]
netCharge=allProteins_seq.columns(feature="net_abs_charge")[0]
proteinDatabase=pd.DataFrame({'Protein':proteinIDs,'Hydropathy':hydropathy,'Net Charge':netCharge,'Aromaticity':aromaticities,'Neutral Charge':neutralCharges,'IsoelectricPoint':isoelectrics,'LengthAminoAcids':lens, 'Instability':instabilities, 'Gravy':gravyVals, 'Molecular Weight':mws})
proteinDatabase[['Helix','Turn','Sheets']] = secondarySturctures
proteinDatabase.index = proteinDatabase['Protein']
proteinDatabase.drop(['Protein'], axis = 1, inplace = True)

proteinDatabase.to_csv(outputDir + 'ProData8_6.csv')

