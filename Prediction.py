# using the files obtained from DEG analysis in R

import pandas as pd
import numpy as np
import seaborn as sns

de = pd.read_csv('/content/Differential_expression_result.csv', index_col='Unnamed: 0')
data = pd.read_csv('/content/normalized_counts.csv',index_col='Unnamed: 0')

data.shape

data = data[data.index.isin(de.index)]

gene = pd.read_csv('/content/Gene_MetaData.csv', index_col='Unnamed: 0')
cov = pd.read_csv('/content/Covariates_Data.csv', index_col='Unnamed: 0')

genes = gene.loc[list(data.index)]['gene_name']
data.set_index(genes, inplace=True)
for i in data.columns:
  if i not in list(cov.index):
    data.drop(i, axis = 1, inplace=True)

cov.drop('definition', axis = 1, inplace = True)
cov = cov.loc[list(data.columns)]

cov['gender'] = [1 if i =='male' else 0 for i in cov['gender']]
cov['race'] = [1 if i == 'white' else 2 if i == 'asian' else 3 if i == 'black or african american' else 4 if i == 'american indian or alaska native' else 0 for i in cov['race']]

merged_data = data.transpose().merge(cov, left_index=True, right_index=True)
merged_data['ajcc_pathologic_stage'].unique()
merged_data['ajcc_pathologic_stage'] = [2 if i == 'Stage IIIA' or i == 'Stage IIIB' else 1 if i == 'Stage IIB' or i == 'Stage IIA' or i == 'Stage II'  else 3 if i == 'Stage IV' else 0 for i in merged_data['ajcc_pathologic_stage']]

x = merged_data.drop(['ajcc_pathologic_stage'], axis=1)
y = merged_data['ajcc_pathologic_stage']

from sklearn.feature_selection import RFECV
from sklearn.linear_model import LogisticRegression

logreg = LogisticRegression()
rfecv = RFECV(estimator=logreg, step=1, cv=5, min_features_to_select=30)
rfecv.fit(x, y)
ranking = rfecv.ranking_

feature_names_rfe = [name for i, name in enumerate(list(x)) if rfecv.support_[i]]

len(feature_names_rfe)

from sklearn.linear_model import Lasso

# create a LassoCV object to perform feature selection based on L1 regularization
lasso = Lasso()

# fit the object to the data and obtain the coefficients of the features
lasso.fit(x, y)
coefficients = lasso.coef_

coefficients

feature_names = [name for i, name in enumerate(list(x)) if lasso.coef_[i]]
len(feature_names)

feature_selected_data = merged_data[list(feature_names_rfe)]

X = x[feature_names_rfe]
y = y

from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix

xtrain, xtest, ytrain, ytest = train_test_split(X, y, test_size=0.33, random_state=42)

stage = pd.DataFrame(y)
stage1 = stage.loc[stage['ajcc_pathologic_stage']==0]
stage2 = stage.loc[stage['ajcc_pathologic_stage']==1]
stage3 = stage.loc[stage['ajcc_pathologic_stage']==2]
stage4 = stage.loc[stage['ajcc_pathologic_stage']==3]

stage1_train = merged_data.loc[stage1.index[:round(278*60/100)]]
stage1_test = merged_data.loc[stage1.index[round(278*60/100):]]

stage2_train = merged_data.loc[stage2.index[:round(len(stage2)*60/100)]]
stage2_test = merged_data.loc[stage2.index[round(len(stage2)*60/100):]]

stage3_train = merged_data.loc[stage3.index[:round(len(stage3)*60/100)]]
stage3_test = merged_data.loc[stage3.index[round(len(stage3)*60/100):]]

stage4_train = merged_data.loc[stage4.index[:round(len(stage4)*60/100)]]
stage4_test = merged_data.loc[stage4.index[round(len(stage4)*60/100):]]

merged_df = pd.merge(stage1_train.transpose(), stage2_train.transpose(), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, stage3_train.transpose(), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, stage4_train.transpose(), left_index=True, right_index=True)

merged_df_test = pd.merge(stage1_test.transpose(), stage2_test.transpose(), left_index=True, right_index=True)
merged_df_test = pd.merge(merged_df_test, stage3_test.transpose(), left_index=True, right_index=True)
merged_df_test = pd.merge(merged_df_test, stage4_test.transpose(), left_index=True, right_index=True)

x_train = merged_df.transpose().drop('ajcc_pathologic_stage', axis = 1)
x_test = merged_df_test.transpose().drop('ajcc_pathologic_stage', axis = 1)

y_train = list(stage.loc[(x_train.index)]['ajcc_pathologic_stage'])
y_test = list(stage.loc[(x_test.index)]['ajcc_pathologic_stage'])

x_train = x_train[feature_names_rfe]
x_test = x_test[feature_names_rfe]

svc = SVC()
svc.fit(xtrain, ytrain)
pred = svc.predict(xtest)

svc = SVC()
svc.fit(x_train, y_train)
pred = svc.predict(x_test)

print(classification_report(y_test, pred))

print(classification_report(y_test, pred))

print(confusion_matrix(y_test,pred))

from sklearn.tree import DecisionTreeClassifier
dt = DecisionTreeClassifier()
dt.fit(xtrain,ytrain)
pred1 = dt.predict(xtest)
print(classification_report(ytest, pred1))

from sklearn.tree import DecisionTreeClassifier
dt = DecisionTreeClassifier()
dt.fit(x_train,y_train)
pred1 = dt.predict(x_test)
print(classification_report(y_test, pred1))

from sklearn.ensemble import RandomForestClassifier
rf = RandomForestClassifier()
rf.fit(xtrain,ytrain)
pred3 = rf.predict(xtest)
print(confusion_matrix(ytest,pred3))
print(classification_report(ytest, pred3))

rf.fit(x_train,y_train)
pred3 = rf.predict(x_test)
print(confusion_matrix(y_test,pred3))
print(classification_report(y_test, pred3))

log = LogisticRegression()
log.fit(x_train,y_train)
pred5 = log.predict(x_test)
print(confusion_matrix(y_test,pred5))
print(classification_report(y_test, pred5))

log = LogisticRegression()
log.fit(xtrain,ytrain)
pred5 = log.predict(xtest)
print(confusion_matrix(ytest,pred5))
print(classification_report(ytest, pred5))

from sklearn.ensemble import GradientBoostingClassifier
boost = GradientBoostingClassifier()
boost.fit(x_train, y_train)
pred6 = boost.predict(x_test)
print(confusion_matrix(y_test,pred6))
print(classification_report(y_test, pred6))

boost.fit(xtrain, ytrain)
pred6 = boost.predict(xtest)
print(confusion_matrix(ytest,pred6))
print(classification_report(ytest, pred6))

from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier

# Define the hyperparameter search space
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [10, 20, 30, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'max_features': ['auto', 'sqrt']
}

# Create a random forest classifier object
rfc = RandomForestClassifier()

# Create a GridSearchCV object and fit it to the data
grid_search = GridSearchCV(estimator=rfc, param_grid=param_grid, cv=5)
grid_search.fit(x_train, y_train)

# Print the best parameters and score
print("Best Parameters:", grid_search.best_params_)
print("Best Score:", grid_search.best_score_)

from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression

# Define the hyperparameter search space
param_grid = {
    'penalty': ['l1', 'l2', 'elasticnet', 'none'],
    'C': [0.1, 1, 10, 100],
    'solver': ['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga']
}

# Create a logistic regression object
lr = LogisticRegression(max_iter=1000)

# Create a GridSearchCV object and fit it to the data
grid_search = GridSearchCV(estimator=lr, param_grid=param_grid, cv=5)
grid_search.fit(x_train, y_train)

# Print the best parameters and score
print("Best Parameters:", grid_search.best_params_)
print("Best Score:", grid_search.best_score_)

fpr1, tpr1, thresholds1 = roc_curve(y_test, boost.predict_proba(x_test)[:,1], pos_label=1)

x_train.shape

from sklearn.metrics import roc_curve, roc_auc_score, auc
from scipy import interp
from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt

cv = StratifiedKFold(n_splits = 5, shuffle= False)
np.random.seed(100)
tprs = []
aucs = []

mean_fpr = np.linspace(0,1,100)
boost.fit(x_train, y_train)
fpr,tpr,t = roc_curve(y_test, boost.predict_proba(x_test)[:,1], pos_label=1)

roc_auc = auc(fpr,tpr)

plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'Manual split (AUC = %0.2f)' % roc_auc)

i= 1
for train,test in cv.split(X,y):
  prediction = boost.fit(X.iloc[train],y.iloc[train]).predict_proba(X.iloc[test]) 
  fpr,tpr,t = roc_curve(list(y.iloc[test]), prediction[:,1], pos_label = 1)
  roc_auc = auc(fpr,tpr)
  tprs.append(interp(mean_fpr,fpr,tpr))
  aucs.append(roc_auc)
  plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
  i = i+1

plt.plot([0,1],[0,1], linestyle = '--', lw = 2, color = 'black')
mean_tpr = np.mean(tprs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)

plt.plot(mean_fpr, mean_tpr, color = 'blue', 
         label = r'Mean ROC (AUC = %0.2f)' %(mean_auc), lw = 2, alpha=1)


plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Gradient Boosting Classifier ROC')
plt.legend(loc = 'lower right')

plt.savefig("Gradient Boosting Classifier ROC.png")
plt.show()

from sklearn.metrics import roc_curve, roc_auc_score, auc
from scipy import interp
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits = 5, shuffle= False)
np.random.seed(100)
tprs = []
aucs = []

mean_fpr = np.linspace(0,1,100)
dt.fit(x_train, y_train)
fpr,tpr,t = roc_curve(y_test, dt.predict_proba(x_test)[:,1], pos_label=1)
roc_auc = auc(fpr,tpr)
plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'Manual split (AUC = %0.2f)' % roc_auc)

i= 1
for train,test in cv.split(X,y):
  prediction = dt.fit(x.iloc[train],y.iloc[train]).predict_proba(x.iloc[test]) 
  fpr,tpr,t = roc_curve(list(y.iloc[test]), prediction[:,1], pos_label = 1)
  roc_auc = auc(fpr,tpr)
  tprs.append(interp(mean_fpr,fpr,tpr))
  aucs.append(roc_auc)
  plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
  i = i+1

plt.plot([0,1],[0,1], linestyle = '--', lw = 2, color = 'black')
mean_tpr = np.mean(tprs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)

plt.plot(mean_fpr, mean_tpr, color = 'blue', 
         label = r'Mean ROC (AUC = %0.2f)' %(mean_auc), lw = 2, alpha=1)


plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Decision Tree Classifier ROC')
plt.legend(loc = 'lower right')

plt.savefig("Decision Tree Classifier ROC.png")
plt.show()

from sklearn.metrics import roc_curve, roc_auc_score, auc
from scipy import interp
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits = 5, shuffle= False)
np.random.seed(100)
tprs = []
aucs = []


log = LogisticRegression(C= 0.1, penalty= 'l1', solver= 'saga')

mean_fpr = np.linspace(0,1,100)
log.fit(x_train, y_train)
fpr,tpr,t = roc_curve(y_test, log.predict_proba(x_test)[:,1], pos_label=1)
#tprs.append(interp(mean_fpr,fpr,tpr))
roc_auc = auc(fpr,tpr)
#aucs.append(roc_auc)
plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'Manual split (AUC = %0.2f)' % roc_auc)

i= 1
for train,test in cv.split(X,y):
  prediction = log.fit(x.iloc[train],y.iloc[train]).predict_proba(x.iloc[test]) 
  fpr,tpr,t = roc_curve(list(y.iloc[test]), prediction[:,1], pos_label = 1)
  roc_auc = auc(fpr,tpr)
  tprs.append(interp(mean_fpr,fpr,tpr))
  aucs.append(roc_auc)
  plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
  i = i+1

plt.plot([0,1],[0,1], linestyle = '--', lw = 2, color = 'black')
mean_tpr = np.mean(tprs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)

plt.plot(mean_fpr, mean_tpr, color = 'blue', 
         label = r'Mean ROC (AUC = %0.2f)' %(mean_auc), lw = 2, alpha=1)


plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Logistic Regression Classifier ROC')
plt.legend(loc = 'lower right')

plt.savefig("Logistic RegressionClassifier ROC.png")
plt.show()

from sklearn.metrics import roc_curve, roc_auc_score, auc
from scipy import interp
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits = 5, shuffle= False)
np.random.seed(100)
tprs = []
aucs = []
rf = RandomForestClassifier(max_depth=30, max_features='auto', min_samples_leaf=2, min_samples_split=5, n_estimators=300)


mean_fpr = np.linspace(0,1,100)
rf.fit(x_train, y_train)
fpr,tpr,t = roc_curve(y_test, rf.predict_proba(x_test)[:,1], pos_label=1)
#tprs.append(interp(mean_fpr,fpr,tpr))
roc_auc = auc(fpr,tpr)
#aucs.append(roc_auc)
plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'Manual split (AUC = %0.2f)' % roc_auc)

i= 1
for train,test in cv.split(X,y):
  prediction = rf.fit(x.iloc[train],y.iloc[train]).predict_proba(x.iloc[test]) 
  fpr,tpr,t = roc_curve(list(y.iloc[test]), prediction[:,1], pos_label = 1)
  roc_auc = auc(fpr,tpr)
  tprs.append(interp(mean_fpr,fpr,tpr))
  aucs.append(roc_auc)
  plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
  i = i+1

plt.plot([0,1],[0,1], linestyle = '--', lw = 2, color = 'black')
mean_tpr = np.mean(tprs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)

plt.plot(mean_fpr, mean_tpr, color = 'blue', 
         label = r'Mean ROC (AUC = %0.2f)' %(mean_auc), lw = 2, alpha=1)


plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Random Forest Classifier ROC')
plt.legend(loc = 'lower right')

plt.savefig("Random Forest Classifier ROC.png")
plt.show()

from sklearn.ensemble import AdaBoostClassifier
ada = AdaBoostClassifier()

from sklearn.metrics import roc_curve, roc_auc_score, auc
from scipy import interp
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits = 5, shuffle= False)
np.random.seed(100)
tprs = []
aucs = []

mean_fpr = np.linspace(0,1,100)
ada.fit(x_train, y_train)
fpr,tpr,t = roc_curve(y_test, ada.predict_proba(x_test)[:,1], pos_label=1)
#tprs.append(interp(mean_fpr,fpr,tpr))
roc_auc = auc(fpr,tpr)
#aucs.append(roc_auc)
plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'Manual split (AUC = %0.2f)' % roc_auc)

i= 1
for train,test in cv.split(X,y):
  prediction = ada.fit(x.iloc[train],y.iloc[train]).predict_proba(x.iloc[test]) 
  fpr,tpr,t = roc_curve(list(y.iloc[test]), prediction[:,1], pos_label = 1)
  roc_auc = auc(fpr,tpr)
  tprs.append(interp(mean_fpr,fpr,tpr))
  aucs.append(roc_auc)
  plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
  i = i+1

plt.plot([0,1],[0,1], linestyle = '--', lw = 2, color = 'black')
mean_tpr = np.mean(tprs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)

plt.plot(mean_fpr, mean_tpr, color = 'blue', 
         label = r'Mean ROC (AUC = %0.2f)' %(mean_auc), lw = 2, alpha=1)


plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('AdaBoostClassifier ROC')
plt.legend(loc = 'lower right')

plt.savefig("AdaBoostClassifier  ROC.png")
plt.show()

from sklearn.naive_bayes import GaussianNB
NB = GaussianNB()

from sklearn.metrics import roc_curve, roc_auc_score, auc
from scipy import interp
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits = 5, shuffle= False)
np.random.seed(100)
tprs = []
aucs = []

mean_fpr = np.linspace(0,1,100)
NB.fit(x_train, y_train)
fpr,tpr,t = roc_curve(y_test, NB.predict_proba(x_test)[:,1], pos_label=1)
#tprs.append(interp(mean_fpr,fpr,tpr))
roc_auc = auc(fpr,tpr)
#aucs.append(roc_auc)
plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'Manual split (AUC = %0.2f)' % roc_auc)

i= 1
for train,test in cv.split(X,y):
  prediction = NB.fit(x.iloc[train],y.iloc[train]).predict_proba(x.iloc[test]) 
  fpr,tpr,t = roc_curve(list(y.iloc[test]), prediction[:,1], pos_label = 1)
  roc_auc = auc(fpr,tpr)
  tprs.append(interp(mean_fpr,fpr,tpr))
  aucs.append(roc_auc)
  plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
  i = i+1

plt.plot([0,1],[0,1], linestyle = '--', lw = 2, color = 'black')
mean_tpr = np.mean(tprs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)

plt.plot(mean_fpr, mean_tpr, color = 'blue', 
         label = r'Mean ROC (AUC = %0.2f)' %(mean_auc), lw = 2, alpha=1)


plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Gaussian NB ROC')
plt.legend(loc = 'lower right')

plt.savefig("Gaussian NB ROC.png")
plt.show()

from sklearn.metrics import roc_curve, roc_auc_score, auc
from scipy import interp
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits = 5, shuffle= False)
np.random.seed(100)
tprs = []
aucs = []

mean_fpr = np.linspace(0,1,100)
svc.fit(x_train, y_train)
fpr,tpr,t = roc_curve(y_test, svc.predict(x_test), pos_label=1)
#tprs.append(interp(mean_fpr,fpr,tpr))
roc_auc = auc(fpr,tpr)
#aucs.append(roc_auc)
plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'Manual split (AUC = %0.2f)' % roc_auc)

i= 1
for train,test in cv.split(X,y):
  prediction = svc.fit(x.iloc[train],y.iloc[train]).predict(x.iloc[test]) 
  fpr,tpr,t = roc_curve(list(y.iloc[test]), prediction, pos_label = 1)
  roc_auc = auc(fpr,tpr)
  tprs.append(interp(mean_fpr,fpr,tpr))
  aucs.append(roc_auc)
  plt.plot(fpr, tpr, lw=2, alpha = 0.3, label = 'ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
  i = i+1

plt.plot([0,1],[0,1], linestyle = '--', lw = 2, color = 'black')
mean_tpr = np.mean(tprs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)

plt.plot(mean_fpr, mean_tpr, color = 'blue', 
         label = r'Mean ROC (AUC = %0.2f)' %(mean_auc), lw = 2, alpha=1)


plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Support Vector Machine ROC')
plt.legend(loc = 'lower right')

plt.savefig("Support Vector Machine ROC.png")
plt.show()
