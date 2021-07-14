#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""    MicrobeID classification model 
       
   
"""

# Futures
from __future__ import print_function


# Built-in/Generic Imports
import pdb

#Libs
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle

from sklearn.cluster import AffinityPropagation
from sklearn.cross_validation import LeaveOneOut
from sklearn.metrics import classification_report
from sklearn.cross_validation import cross_val_predict
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score


__author__ = 'Alvaro Perdones-Montero'
__copyright__ = 'Copyright 2018, MicrobeID'
__credits__ = ['Alvaro Perdones-Montero','Simon Cameron']
__license__ = 'GPL'
__version__ = '1.0.0'
__maintainer__ = 'Alvaro Perdones-Montero'
__email__ = 'a.perdones-montero@imperial.ac.uk'
__status__ = 'Development'

# Save file with the probability of the classification with the prediction and the actual class
def save_probabilities(probabilities_max,index,prediction,output_file):
    if len(probabilities_max) == len(index) == len(prediction):
        target = open(output_file, 'a')
	target.write("================ PROBABILITIES ================\n")
        for i in range(0,len(probabilities_max)):
            target.write(str(index[i])+ "\t" + str(probabilities_max[i]) + "\t" + str(prediction[i]) + "\n")
        target.close()

# Save file with the model in pickle format that can be loaded for classification of unknown samples
def save_model(model,modelname):

    with open(modelname, 'wb') as output:
        pickle.dump(model, output, pickle.HIGHEST_PROTOCOL)

# Save file with the model and the affinity clusters
def save_model_and_affinityclusters(model_and_affinityclusters,modelname):

    with open(modelname, 'wb') as output:
        pickle.dump(model_and_affinityclusters, output, pickle.HIGHEST_PROTOCOL)
		
# Save file with the precision, recall and F1 for each class
def save_results(results,corrected,output_file,report_string):
   target = open(output_file, 'a')
   target.write("================ PRECISION AND RECALL ================\n")
   target.write(results + "\n")
   target.write("================ " + report_string + " ================\n")
   target.write(str(corrected) + "\n")
   target.close()
   
# Calculates TPR (True positive rate)
def correct_classification_rate(y_actual, y_predicted):
   number_passes = 0
   for i in range(len(y_actual)):
      if y_actual[i] == y_predicted[i]:
         number_passes = number_passes + 1
   total_samples = len(y_actual)
   print(total_samples)
   print(number_passes)
   correct_perc = 100*(float(number_passes)/float(total_samples))

   return correct_perc

# Plot confusion matrix   
def plot_confusion_matrix(cm, classes, title='Confusion matrix', cmap=plt.cm.Blues):
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=90, fontsize ='small')
    plt.yticks(tick_marks, classes, fontsize = 'small')
    #plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.tight_layout()

# Calculates all paired ratios in the input matrix 	
def parallel_ratio(X_filtered_veryintense, init, end):

    X_filtered_ratio = X_filtered_veryintense.mean(axis=1,numeric_only=True)
    for i in range(init, end -1):
	partdf = X_filtered_veryintense.ix[:,i+1:].div(X_filtered_veryintense.ix[:,i],axis="rows")
        mz_intense = X_filtered_veryintense.columns[i]
        partdf.columns = X_filtered_veryintense.ix[:,i+1:].columns + "00000" + str(int(float(mz_intense) * 10))
        X_filtered_ratio = pd.concat([X_filtered_ratio, partdf], axis=1)

    return X_filtered_ratio

# Calculates all paired ratios in the input matrix.
# The column name format is A__vs__B
def parallel_ratio_2(X_filtered_veryintense, init, end):

    X_filtered_ratio = X_filtered_veryintense.mean(axis=1,numeric_only=True)
    for i in range(init, end):
	partdf = X_filtered_veryintense.div(X_filtered_veryintense.ix[:,i],axis="rows")
        mz_intense = X_filtered_veryintense.columns[i]
	partdf.columns = X_filtered_veryintense.columns + "__vs__" + mz_intense
	X_filtered_ratio = pd.concat([X_filtered_ratio, partdf], axis=1)

    return X_filtered_ratio

##################
#                #
#    MAIN CODE   #
#                # 
##################	


# Read csv file OMB 997 format
X = pd.read_csv("matrix_for_postprocessing.val.csv", delimiter=",",index_col=0)

print("Input validation file readed")

# Select scan for burn
# It depends on the OMB burn selection
burns_start_scan = [1]


header = X.columns.values

# OMB 997 format includes 5 columns as metadata
# Class, Filename, Start scan, End scan, Sum.
# Metadata is selected
X_metadata = X.iloc[:,0:4]
classes = X_metadata.index
samples = X_metadata.ix[:,'Sample']

# Creates cross-reference between Filename and Class
# Shouldn't be the same file name several times in the input file
sample_class_cross = dict()
for s,c in zip(samples,classes):
	if s not in sample_class_cross:
		sample_class_cross[s] = c
	else:
		if c != sample_class_cross[s]:
			print("ERROR PLEASE CHECK DATA MATRIX")
			break

samples_grouped = []
classes_grouped = []			
for s in sample_class_cross:
	samples_grouped.append(s)
	classes_grouped.append(sample_class_cross[s])

sample_class_cross_grouped = pd.DataFrame(classes_grouped,index=samples_grouped)	


# Select data only from input file
X = X.iloc[:,4:]
X_rows = list(X.index)
X_filtered_row = X

X_filtered = X_filtered_row

model = pickle.load(open("RFModel.model","rb"))

# Detection of bins representating the data set using Affinity Propagation
# Used for reducing the dataset for ratio analysis
#af = AffinityPropagation()
#clusters = af.fit(X_filtered.transpose()).cluster_centers_indices_
#clusters_features = X_filtered.columns.values[clusters]


X_filtered = X_filtered.ix[:,model[1]]


# Select the entries in the input file that match the scans selected
# Can be used for having different regions in the input files and select
# only the ones that we want.
# As example, we could have an input file including a background region
# in that case, we should remove those entries for the analysis
X_filtered_2 = X_filtered.ix[X_metadata.loc[:,"Start scan"].isin(burns_start_scan),:]
X_metadata_2 = X_metadata.ix[X_metadata.loc[:,"Start scan"].isin(burns_start_scan),:]

# Includes again the metadata
X_filtered_3 = pd.concat([X_metadata_2, X_filtered_2], axis=1)

# If we have more than one regions in the input file, per example, if 
# we have the 3 burns, we do the sum up of all of them
X_filtered_4 = X_filtered_3.groupby(['Sample'],as_index=False).sum()

# Uses the cross reference for assign a class to each entry
# after doing the sum up
classes_aftergrouping = []
for l in X_filtered_4['Sample']:
	classes_aftergrouping.append(sample_class_cross[l])
	
X_filtered_4.index = classes_aftergrouping

# Removes the metadata
X_filtered = X_filtered_4.iloc[:,4:]

samples_names = X_filtered_4['Sample']
classes_actual = X_filtered_4.index
# Ratio analysis
#X_filtered_ratio = X_filtered.mean(axis=1,numeric_only=True)
#X_filtered_veryintense = X_filtered.copy()
#X_filtered_ratio = parallel_ratio_2(X_filtered_veryintense,0,X_filtered_veryintense.shape[1])
#X_filtered = X_filtered_ratio.replace(np.nan,0).replace([np.inf,-np.inf],10000000000000000).copy()
X_scikit_filtered = X_filtered


# Starts classification model building
y = X_filtered.index
y_scikit = np.ravel(y)
number_samples = len(y_scikit)
groups = y_scikit
classes = list(set(groups))


predicted = model[0].predict(X_scikit_filtered)
predicted_proba = model[0].predict_proba(X_scikit_filtered)
predicted_proba_df = pd.DataFrame(data = predicted_proba, index = samples_names, columns = model[0].classes_)
predicted_proba_df.ix[:,"Actual"] = classes_actual
predicted_proba_df.to_csv("probabilities.csv")
# Calculate accuracy
results = accuracy_score(groups, predicted)
# Calculate precision, recall, F1
prec_recall = classification_report(groups, predicted)
# Calculate confusion matrix
confusion = confusion_matrix(groups, predicted,labels=classes)
# Calcalate TPR
corrected = correct_classification_rate(groups, predicted)

# Save confusion matrix as png
output_file_confusion = "confusion_matrix.png"
confusion_normalized = confusion.astype('float') / confusion.sum(axis=1)[:, np.newaxis]
plt.figure()
plot_confusion_matrix(confusion_normalized, classes, title='Normalized confusion matrix')
plt.savefig(output_file_confusion)

# Save classification results for the CV
output_report_filename = "classification_results.txt"
save_results(prec_recall,corrected,output_report_filename,"TPR")


