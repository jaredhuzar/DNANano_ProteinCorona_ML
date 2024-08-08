import pandas as pd
import sys

outputDir = 'AnalyzedDataTest\\'

proDataPath=sys.argv[1]
proteinDatabase=pd.read_csv(proDataPath,index_col=0)

functions1Path=sys.argv[2]
functions1=pd.read_csv(functions1Path,index_col=0)

goFlexPath=sys.argv[3]
goFlex=pd.read_csv(goFlexPath,index_col=0)

InteractionsPath=sys.argv[4]
Interactions=pd.read_csv(InteractionsPath,index_col=0)

moreNetSurfPath=sys.argv[5]
moreNetSurf=pd.read_csv(moreNetSurfPath,index_col=0)

totPath=sys.argv[6]
tot=pd.read_csv(totPath,index_col=0)

totaaPath=sys.argv[7]
totaa=pd.read_csv(totaaPath,index_col=0)

totListaaExpoPath=sys.argv[8]
totListaaExpo=pd.read_csv(totListaaExpoPath,index_col=0)

combosPath=sys.argv[9]
combos=pd.read_csv(combosPath,index_col=0)


proteinDatabase=proteinDatabase.loc[proteinDatabase.index.isin(set(combos.index).intersection(set(proteinDatabase.index)))]
functions1=functions1.loc[functions1.index.isin(set(combos.index).intersection(set(functions1.index)))]
goFlex=goFlex.loc[goFlex.index.isin(set(combos.index).intersection(set(goFlex.index)))]
Interactions=Interactions.loc[Interactions.index.isin(set(combos.index).intersection(set(Interactions.index)))]
moreNetSurf=moreNetSurf.loc[moreNetSurf.index.isin(set(combos.index).intersection(set(moreNetSurf.index)))]
tot=tot.loc[tot.index.isin(set(proteinDatabase.index).intersection(set(tot.index)))]
totaa=totaa.loc[totaa.index.isin(set(proteinDatabase.index).intersection(set(proteinDatabase.index)))]
totListaaExpo=totListaaExpo.loc[totListaaExpo.index.isin(set(proteinDatabase.index).intersection(set(proteinDatabase.index)))]
combos=combos.loc[combos.index.isin(set(proteinDatabase.index).intersection(set(combos.index)))]
moreNetSurf=moreNetSurf.loc[moreNetSurf.index.isin(set(proteinDatabase.index).intersection(set(moreNetSurf.index)))]

mergedPros=pd.concat([proteinDatabase, goFlex], axis=1)
mergedPros=pd.concat([mergedPros, functions1], axis=1)
mergedPros=pd.concat([mergedPros,totaa], axis = 1)
mergedPros=pd.concat([mergedPros, totListaaExpo], axis =1)
mergedPros=pd.concat([mergedPros,combos], axis =1)
mergedPros=pd.concat([mergedPros,moreNetSurf], axis =1)
mergedPros=pd.concat([mergedPros,tot], axis =1)
mergedPros=pd.concat([mergedPros,Interactions], axis =1)

mergedPros.to_csv(outputDir + 'ProStatistics8_7.csv')

