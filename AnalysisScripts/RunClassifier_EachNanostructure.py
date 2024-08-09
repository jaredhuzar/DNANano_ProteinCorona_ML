import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sklearn
import pandas as pd
from sklearn.tree import DecisionTreeClassifier # Import Decision Tree Classifier
from sklearn.model_selection import train_test_split # Import train_test_split function
from sklearn import metrics
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import KFold
import math
import xgboost as xgb
import seaborn as sns
import shap
from scipy import stats
import imblearn
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay 

sns.set_style('white')
sns.set_context("paper", font_scale = 2)

outputDir = 'Results\\'

dataPresentPath=sys.argv[1]
dataPresent=pd.read_csv(dataPresentPath)

overallRes3 = dataPresent[['Present', 'Sample', 'Protein']]
overallRes3=overallRes3.rename(columns={"Protein":'ID'})
overallRes3=overallRes3.rename(columns={"Present":'Abundance'})
dataPresent.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0'],inplace = True, axis =1)

classifier = xgb.XGBClassifier()
structs = list()
f1_all = list()
acc_structs = list()
prec_all = list()
rec_all = list()
allScores = list()
featureimportancesTot = list()
allAUCs = list()
for z in np.unique(overallRes3['Sample']):
    subRes = overallRes3[overallRes3['Sample'] == z]
    abundanceTest1 = subRes['Abundance']
    subX=dataPresent.loc[subRes.index]
    subX = subX.reset_index(drop = True)
    subY= np.asarray(abundanceTest1)
    #pros = list()
    #for x in range(len(abundanceBin)):
    #    bool1=abundanceBin[x]
    #    if(bool1 == 1):
    #        pros.append(subRes.loc[x,'Proteins'])
    #pd.DataFrame(pros).to_csv(z + '_PresentPros.csv', index = None, header = None)
    kFold=KFold(n_splits=10,shuffle=True, random_state=42)
    clf = classifier

    features = list()
    scores = list()
    f1_structs = list()
    prec_structs = list()
    rec_structs = list()
    aucScores = list()
    for train_index,test_index in kFold.split(subX):
        X_train, X_test, y_train, y_test = subX.loc[train_index], subX.loc[test_index], subY[train_index], subY[test_index]
        clf = clf.fit(X_train,y_train)
        y_pred = clf.predict(X_test)
        scores.append(metrics.accuracy_score(y_test, y_pred))
        f1_structs.append(metrics.f1_score(y_test, y_pred))
        prec_structs.append(metrics.precision_score(y_test, y_pred))
        rec_structs.append(metrics.recall_score(y_test, y_pred))
        features.append(clf.feature_importances_)
        y_pred = clf.predict_proba(X_test)
        aucScores.append(roc_auc_score(y_test,y_pred[:,1]))
        structs.append(z)
    #explainer = shap.TreeExplainer(clf, X_train)
    #shap_values = explainer.shap_values(X_test)
    #plt.figure()
    #plt.title(z)
    #shap.summary_plot(shap_values, X_test.astype("float"), show=False)
    #plt.savefig(z + ' SHAP.png')
    allScores.append(scores)
    allAUCs.append(aucScores)
    f1_all.append(f1_structs)
    featureimportancesTot.append(features)
    prec_all.append(prec_structs)
    rec_all.append(rec_structs)

allStructsMet=pd.DataFrame({'Structures':structs,'AUC':np.array(allAUCs).flatten(),'F1':np.array(f1_all).flatten(),'Accuracy':np.array(allScores).flatten(),'Precision':np.array(prec_all).flatten(),'Recall':np.array(rec_all).flatten()})

allStructsMet=pd.DataFrame({'Structures':structs,'AUC':np.array(allAUCs).flatten(),'F1':np.array(f1_all).flatten(),'Accuracy':np.array(allScores).flatten(),'Precision':np.array(prec_all).flatten(),'Recall':np.array(rec_all).flatten()})

vals = list()
metricname  = list()
DataSet = list()
for i in allStructsMet.columns[1:6]:
    vals1=allStructsMet[i]
    for z in range(len(vals1)):
        DataSet.append(allStructsMet.loc[z,'Structures'])
        vals.append(vals1[z])
        metricname.append(i)

aucVal = list()
f1Val = list()
accVal = list()
precVal = list()
recVal = list()
for z in np.unique(allStructsMet['Structures']):
    structDF=allStructsMet[allStructsMet['Structures'] == z]
    f1Val.append(np.mean(structDF['F1']))
    accVal.append(np.mean(structDF['Accuracy']))
    precVal.append(np.mean(structDF['Precision']))
    aucVal.append(np.mean(structDF['AUC']))
    recVal.append(np.mean(structDF['Recall']))

allStructsMet=pd.DataFrame({'Structures':np.unique(allStructsMet['Structures']),'AUC':aucVal,'F1':f1Val,'Accuracy':accVal,'Precision':precVal,'Recall':recVal})

fig = plt.figure(figsize=(4, 3))
plt.rc('font', family='serif')
plt.rc('font', family='serif')
#plt.rc('xtick', labelsize='large')
#plt.rc('ytick', labelsize='large')

#vals = list()
#metricname  = list()
#DataSet = list()
#for i in allStructsMet.columns[1:6]:
#    vals1=allmet[i]
#    for z in range(len(vals1)):
#        DataSet.append(allStructsMet.loc[z,'Structures'])
#        vals.append(vals1[z])
#        metricname.append(i)
everymet=pd.DataFrame({'Scores':vals, 'Metric':metricname, 'Structure':DataSet})

#met = allmet
labs=["{:.2f}".format(np.mean(allStructsMet['AUC'])),"{:.2f}".format(np.mean(allStructsMet['F1'])),"{:.2f}".format(np.mean(allStructsMet['Accuracy'])),"{:.2f}".format(np.mean(allStructsMet['Recall'])),"{:.2f}".format(np.mean(allStructsMet['Precision']))]
fig, ax = plt.subplots()
g=sns.barplot(data = everymet,x='Metric', y = 'Scores',palette='colorblind', saturation = .7)
ax.set_ylim(0, 1)
for i in ax.containers:
    ax.bar_label(i,labels = labs,padding=-30)
#plt.title('All Information Enriched vs Not Enriched')
g.set(ylim=(0, 1), xlabel ="Metrics", ylabel = "Scores", title = 'Per Nanostructure Performance')
plt.xticks(rotation = 45)
plt.savefig(outputDir + 'PerNanostructurePerformance.png', bbox_inches='tight')
plt.show()

