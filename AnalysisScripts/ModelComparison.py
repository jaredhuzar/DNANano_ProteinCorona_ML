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
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_regression
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
from sklearn.metrics import RocCurveDisplay, auc

sns.set_style('white')
sns.set_context("paper", font_scale = 2)

outputDir = 'Results\\'

inputFile=sys.argv[1]

dataPresent = pd.read_csv(inputFile)
abundance_present = dataPresent['Present']
abundance_present = np.asarray(abundance_present)
dataPresent.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0'],inplace = True, axis =1)
dataNew = dataPresent
classifier = xgb.XGBClassifier()

Y = np.asarray(abundance_present)

clf = classifier
scores = list()
features = list()
acc=list()
prec=list()
rec = list()
aucScores=list()
kFold=KFold(n_splits=10,shuffle=True, random_state=1)
for train_index,test_index in kFold.split(dataNew):
    X_train, X_test, y_train, y_test = dataNew.loc[train_index], dataNew.loc[test_index], Y[train_index], Y[test_index]
    clf = clf.fit(X_train,y_train)
    y_pred = clf.predict(X_test)
    #scores.append(metrics.accuracy_score(y_test, y_pred))
    scores.append(metrics.f1_score(y_test, y_pred))
    acc.append(metrics.accuracy_score(y_test, y_pred))
    features.append(clf.feature_importances_)
    prec.append(metrics.precision_score(y_test, y_pred))
    rec.append(metrics.recall_score(y_test, y_pred))
    y_pred = clf.predict_proba(X_test)
    aucScores.append(roc_auc_score(y_test,y_pred[:,1]))

allmet=pd.DataFrame({'AUC':aucScores,'F1':scores, 'Accuracy':acc, 'Recall':rec, 'Precision':prec})
allmet['Model'] = 'XGBoost'

classifier = RandomForestClassifier()
clf = classifier
scores = list()
features = list()
acc=list()
prec=list()
rec = list()
aucScores=list()
kFold=KFold(n_splits=10,shuffle=True, random_state=1)
for train_index,test_index in kFold.split(dataNew):
    X_train, X_test, y_train, y_test = dataNew.loc[train_index], dataNew.loc[test_index], Y[train_index], Y[test_index]
    clf = clf.fit(X_train,y_train)
    y_pred = clf.predict(X_test)
    #scores.append(metrics.accuracy_score(y_test, y_pred))
    scores.append(metrics.f1_score(y_test, y_pred))
    acc.append(metrics.accuracy_score(y_test, y_pred))
    features.append(clf.feature_importances_)
    prec.append(metrics.precision_score(y_test, y_pred))
    rec.append(metrics.recall_score(y_test, y_pred))
    y_pred = clf.predict_proba(X_test)
    aucScores.append(roc_auc_score(y_test,y_pred[:,1]))

rfMet=pd.DataFrame({'AUC':aucScores,'F1':scores, 'Accuracy':acc, 'Recall':rec, 'Precision':prec})
rfMet['Model'] = 'Random Forest'

classifier = GradientBoostingClassifier()
clf = classifier
scores = list()
features = list()
acc=list()
prec=list()
rec = list()
aucScores=list()
kFold=KFold(n_splits=10,shuffle=True, random_state=1)
for train_index,test_index in kFold.split(dataNew):
    X_train, X_test, y_train, y_test = dataNew.loc[train_index], dataNew.loc[test_index], Y[train_index], Y[test_index]
    clf = clf.fit(X_train,y_train)
    y_pred = clf.predict(X_test)
    #scores.append(metrics.accuracy_score(y_test, y_pred))
    scores.append(metrics.f1_score(y_test, y_pred))
    acc.append(metrics.accuracy_score(y_test, y_pred))
    features.append(clf.feature_importances_)
    prec.append(metrics.precision_score(y_test, y_pred))
    rec.append(metrics.recall_score(y_test, y_pred))
    y_pred = clf.predict_proba(X_test)
    aucScores.append(roc_auc_score(y_test,y_pred[:,1]))

gradientMet=pd.DataFrame({'AUC':aucScores,'F1':scores, 'Accuracy':acc, 'Recall':rec, 'Precision':prec})
gradientMet['Model'] = 'Gradient Boosting'

metricComparison = pd.concat([allmet, rfMet, gradientMet])
metricComparison.reset_index(drop=True, inplace = True)

vals = list()
metricname  = list()
DataSet = list()
for i in metricComparison.columns[0:5]:
    vals1=metricComparison[i]
    for z in range(len(vals1)):
        DataSet.append(metricComparison.loc[z,'Model'])
        vals.append(vals1[z])
        metricname.append(i)
everymet=pd.DataFrame({'Scores':vals, 'Metric':metricname, 'Model':DataSet})

fig = plt.figure(figsize=(4, 3))
plt.rc('font', family='serif')
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')
g = sns.catplot(
    data=everymet, kind="bar",
    x="Model", y="Scores", hue="Metric",
    palette='colorblind', alpha=.6, height=6,saturation = .7
)
#g.despine(left=True)
g.set(ylim=(0, 1), title='Performance Across Models')
g.set_axis_labels("", "Score")
plt.xticks(rotation=45)
g.legend.set_title("")
plt.savefig(outputDir + 'comparingModels.png',bbox_inches='tight')
everymet.to_csv(outputDir + 'ModelComparison.csv')

