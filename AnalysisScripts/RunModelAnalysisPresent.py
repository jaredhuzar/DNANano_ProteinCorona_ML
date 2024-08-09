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
from sklearn.metrics import RocCurveDisplay, auc

sns.set_style('white')
sns.set_context("paper", font_scale = 2)

inputDir = 'AnalyzedData\\'
outputDir = 'Results\\'

inputFile=sys.argv[1]

dataPresent = pd.read_csv(inputFile)
abundance_present = dataPresent['Present']
abundance_present = np.asarray(abundance_present)
dataPresent.drop(['Present', 'Protein', 'Sample', 'Unnamed: 0'],inplace = True, axis =1)



def train_run_model(xData,yData, label1):
    classifier = xgb.XGBClassifier()
    data = xData
    Y = yData

    clf = classifier
    scores = list()
    features = list()
    acc=list()
    prec=list()
    rec = list()
    aucScores = list()
    b = 0
    kFold=KFold(n_splits=10,shuffle=True, random_state=1)
    for train_index,test_index in kFold.split(data):
        X_train, X_test, y_train, y_test = data.loc[train_index], data.loc[test_index], Y[train_index], Y[test_index]
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
        if(b == 0):
            fets = pd.DataFrame(clf.feature_importances_)
        else:
            fets=pd.concat([fets,pd.DataFrame(clf.feature_importances_)], ignore_index = True, axis = 1)
        b=b+1
        
    allmet=pd.DataFrame({'AUC':aucScores,'F1':scores, 'Accuracy':acc, 'Recall':rec, 'Precision':prec})
    allmet['DataSet'] = label1
    fig = plt.figure(figsize=(4, 3))
    plt.rc('font', family='serif')
    plt.rc('font', family='serif')
    #plt.rc('xtick', labelsize='large')
    #plt.rc('ytick', labelsize='large')

    vals = list()
    metricname  = list()
    DataSet = list()
    for i in allmet.columns[0:5]:
        vals1=allmet[i]
        for z in range(len(vals1)):
            DataSet.append(allmet.loc[z,'DataSet'])
            vals.append(vals1[z])
            metricname.append(i)
    everymet=pd.DataFrame({'Scores':vals, 'Metric':metricname, 'Model':DataSet})

    met = allmet
    labs=["{:.2f}".format(np.mean(met['AUC'])),"{:.2f}".format(np.mean(met['F1'])),"{:.2f}".format(np.mean(met['Accuracy'])),"{:.2f}".format(np.mean(met['Recall'])),"{:.2f}".format(np.mean(met['Precision']))]
    fig, ax = plt.subplots()
    g=sns.barplot(data = everymet,x='Metric', y = 'Scores',palette='colorblind', saturation = .7)
    ax.set_ylim(0, 1)
    for i in ax.containers:
        ax.bar_label(i,labels = labs,padding=-30)
    #plt.title('All Information Enriched vs Not Enriched')
    g.set(ylim=(0, 1), xlabel ="Metrics", ylabel = "Scores", title ='Model Performance')
    plt.xticks(rotation = 45)
    plt.savefig(outputDir + label1 + 'PresentPerformance.png', bbox_inches = 'tight')
    plt.clf()
    plt.show()
    
    explainer = shap.Explainer(clf, data)
    shap_values = explainer.shap_values(data)
    shap.summary_plot(shap_values, data.astype("float"), show = False)
    plt.savefig(outputDir + label1+'ShapPresent.png', bbox_inches = 'tight')
    plt.clf()
    plt.show()
    
    clf = classifier
    scores = list()
    features = list()
    acc=list()
    cv=KFold(n_splits=10,shuffle=True, random_state=1)
    y=np.asarray(yData)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    n_splits=10
    fig, ax = plt.subplots(figsize=(6, 6))
    for fold, (train, test) in enumerate(cv.split(data, y)):
        classifier.fit(data.loc[train], y[train])
        viz = RocCurveDisplay.from_estimator(
            classifier,
            data.loc[test],
            y[test],
            name=f"ROC fold {fold+1}",
            alpha=0.3,
            lw=1,
            ax=ax,
            #plot_chance_level=(fold == n_splits - 1),
        )
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(
        mean_fpr,
        tprs_lower,
        tprs_upper,
        color="grey",
        alpha=0.2,
        label=r"$\pm$ 1 std. dev.",
    )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
        title=f"Mean ROC curve with variability Present Protein Prediction",
    )

    plt.plot([0, 1], [0, 1], color="navy", lw=1, linestyle="--")
    ax.axis("square")
    ax.legend(loc="lower right")
    plt.savefig(outputDir + label1 + 'PresentAUC.png',bbox_inches='tight')
    plt.clf()
    return [everymet,fets,shap_values]

def getFeatureImportances(featureImportance, data):
    medianFeatures=featureImportance.median(axis =1)
    retFeatures=pd.DataFrame(medianFeatures)
    retFeatures.index = data.columns
    return retFeatures

statsPres,featureImportancePres,shapsPres =train_run_model(dataPresent,abundance_present,'Present')

presentImportance=getFeatureImportances(featureImportancePres,dataPresent)

presentImportance.to_csv(outputDir + 'PresentFeatureImportance.csv')
statsPres.to_csv(outputDir + 'PresentStats.csv')

