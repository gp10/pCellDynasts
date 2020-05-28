#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 14:51:05 2017

@author: gp10
"""

import pandas as pd
import seaborn as sns
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import scipy.io as sio
import os


######################################################################################
""" SCRIPT USED FOR PLOTTING EXPERIMENTAL H2BGFP INTENSITY DISTRIBUTIONS OVER TIME """
######################################################################################
## An experimental data set is selected and histograms of individual-cell H2BGFP intensities plotted


""" SELECTION OF EXPERIMENTAL DATA SET AND OUTPUT SETTINGS: """
selectDataSet = 1 # ( 1=OE_R26rtTA | 2=Paw_R26rtTA | 3=Ear_R26rtTA | 4=Back_R26rtTA | 5=TailScale_R26rtTA | 6=TailInterscale_R26rtTA )
typePlot = "Boxplot" # ( "Boxplot" | "Violinplot" )


""" LOAD EXPERIMENTAL H2BGFP DATA """
switcher = {
    1: sio.loadmat(os.path.join(os.getcwd(),'Datasets/H2BGFPdil_OE_R26rtTA_dataset.mat')),
    2: sio.loadmat(os.path.join(os.getcwd(),'Datasets/H2BGFPdil_Paw_R26rtTA_dataset.mat')),
    3: sio.loadmat(os.path.join(os.getcwd(),'Datasets/H2BGFPdil_Ear_R26rtTA_dataset.mat')),
    4: sio.loadmat(os.path.join(os.getcwd(),'Datasets/H2BGFPdil_Back_R26rtTA_dataset.mat')),
    5: sio.loadmat(os.path.join(os.getcwd(),'Datasets/H2BGFPdil_TailScale_R26rtTA_dataset.mat')),
    6: sio.loadmat(os.path.join(os.getcwd(),'Datasets/H2BGFPdil_TailInterscale_R26rtTA_dataset.mat')),
}
matVars = switcher.get(selectDataSet, 'invalid data set')
nmice = matVars['nmice'][0]
rtime = matVars['rtime'][0]
H2BGFP_all = matVars['H2BGFP_all'][0]
H2BGFP = matVars['H2BGFP']
FieldView = matVars['FieldView']
if selectDataSet != 4: # All data sets except Back_R26rtTA come from the same experiment with same collection time points:
    ntime = ['0d','7d','12d','18d']
    experim = ['A','A','A','A']
else: # Back_R26rtTA data set comes from 3 experiments structured in different time points:
    ntime = ['0d_a','5d','14d','0d_b','11d','0d_c','21d']
    experim = ['A','A','A','B','B','C','C']


""" log2 TRANSFORM AND REARRANGE DATA FORMAT """
# log2 transform and reorganize dataset as a list (removing possible NaNs):
list_log2_H2BGFP = []
list_fieldNo = []
miceName = []
timeName = []
timeCount = []
experimName = []
for aj in range(len(rtime)):
    for ej in range(nmice[aj]):
        #Data:
        mycol = H2BGFP[ej,aj]
        mycol = mycol.astype(float)
        mycol = np.array(np.log2(mycol[np.isnan(mycol)!=1]))
        list_log2_H2BGFP.append(mycol)
        
        #FieldView:
        mycol2 = FieldView[ej,aj]
        mycol2 = np.array(mycol2[np.isnan(mycol2.astype(float))!=1])
        list_fieldNo.append(mycol2)
        
        #Mouse:
        miceName.append( ntime[aj] + '.' + str(ej+1) )
        
        #Time:
        timeName.append(ntime[aj])
        timeCount.append(rtime[aj]*7) # time is originally given in weeks and here converted into days
        
        #Experiment:
        experimName.append(experim[aj])
    
# Normalize data by the mean H2BGFP intensity at time = 0
list_nlog2_H2BGFP = []
if selectDataSet != 4: # All data sets except Back_R26rtTA come from the same experiment with same collection time points:
    mean_log2_H2BGFP_t0 = np.mean([np.mean(x) for x in list_log2_H2BGFP[0:nmice[0]]])
    list_nlog2_H2BGFP[:] = [x - mean_log2_H2BGFP_t0 for x in list_log2_H2BGFP]
else: # Back_R26rtTA data set comes from 3 experiments structured in different time points:
    cumsum_nmice = np.cumsum(nmice)
    mean_log2_H2BGFP_t0_setA = np.mean([np.mean(x) for x in list_log2_H2BGFP[0:nmice[0]]]) # experiment KM41 (t=0, t=5, t=14)
    list_nlog2_H2BGFP[0:cumsum_nmice[2]] = [x - mean_log2_H2BGFP_t0_setA for x in list_log2_H2BGFP[0:cumsum_nmice[2]]]
    mean_log2_H2BGFP_t0_setB = np.mean([np.mean(x) for x in list_log2_H2BGFP[cumsum_nmice[2]:cumsum_nmice[2]+nmice[3]]]) # experiment KM29 (t=0, t=11)
    list_nlog2_H2BGFP[cumsum_nmice[2]:cumsum_nmice[4]] = [x - mean_log2_H2BGFP_t0_setB for x in list_log2_H2BGFP[cumsum_nmice[2]:cumsum_nmice[4]]]
    mean_log2_H2BGFP_t0_setC = np.mean([np.mean(x) for x in list_log2_H2BGFP[cumsum_nmice[4]:cumsum_nmice[4]+nmice[5]]]) # experiment KM24 (t=0, t=21)
    list_nlog2_H2BGFP[cumsum_nmice[4]:cumsum_nmice[6]] = [x - mean_log2_H2BGFP_t0_setC for x in list_log2_H2BGFP[cumsum_nmice[4]:cumsum_nmice[6]]]
    
# Compact/rearrange data as a single-vector list:
All_list_nlog2_H2BGFP = []
All_list_fieldNo = []
All_miceName = []
All_timeName = []
All_timeCount = []
All_experimName = []
for aj in range(len(list_nlog2_H2BGFP)):
    All_list_nlog2_H2BGFP.extend(list_nlog2_H2BGFP[aj])
    All_list_fieldNo.extend(list_fieldNo[aj])
    All_miceName.extend(np.repeat(miceName[aj],len(list_nlog2_H2BGFP[aj])))
    All_timeName.extend(np.repeat(timeName[aj],len(list_nlog2_H2BGFP[aj])))
    All_timeCount.extend(np.repeat(timeCount[aj],len(list_nlog2_H2BGFP[aj])))
    All_experimName.extend(np.repeat(experimName[aj],len(list_nlog2_H2BGFP[aj])))
All_timeCount[:] = [int(x) for x in All_timeCount]

# Reallocate data set into a dictionary:
all_Data = {"H2BGFP intensity": All_list_nlog2_H2BGFP, "field": All_list_fieldNo, "animal": All_miceName, "timeName": All_timeName, "time": All_timeCount, "experiment": All_experimName}
my_df = pd.DataFrame(all_Data)


""" PLOTTING EXPERIMENTAL H2BGFP INTENSITY DISTRIBUTIONS """
# Plot the whole data set (pooled per animal):
sns.set()
fig = plt.figure()
ax = plt.gca()
if selectDataSet != 4: # All data sets except Back_R26rtTA
    if typePlot == "Boxplot":
        # boxplot:
        flierprops = dict(marker='o', markerfacecolor='grey', markersize=1, linestyle='none')
        ax = sns.boxplot(x="time", y="H2BGFP intensity", hue="animal", data=my_df, dodge=True, color="green", showfliers=True, flierprops=flierprops) #fliersize=1
    else:
        # violinplot:
        ax = sns.factorplot(x="time", y="H2BGFP intensity", data=my_df, hue = "animal", kind="violin", color="green", split=False, inner="quartile", scale="area", size=4, aspect=2, cut=0)
    # swarmplot:
    #ax = sns.swarmplot(x="time", y="H2BGFP intensity", data=my_df, hue="animal", dodge=True, size=1, color="lightgray", linewidth=0)
else: # Back_R26rtTA data set
    if typePlot == "Boxplot":
        # boxplot:
        flierprops = dict(marker='o', markerfacecolor='grey', markersize=1, linestyle='none')
        ax = sns.factorplot(x="time", y="H2BGFP intensity", data=my_df, hue="animal", row="experiment", kind="box", color="green", showfliers=True, flierprops=flierprops, size=4, aspect=2)
    else:
        # violinplot:
        ax = sns.factorplot(x="time", y="H2BGFP intensity", data=my_df, hue="animal", row="experiment", kind="violin", color="green", split=False, inner="quartile", scale="area", size=4, aspect=2, cut=0)
    # alternative violinplot:
    #fig, axes = plt.subplots(nrows=3,ncols=1, figsize=(8,8))
    #sns.factorplot(x="time", y="H2BGFP intensity", data=my_df[my_df.experiment == 'A'], hue = "animal", kind="violin", color="green", split=False, inner="quartile", scale="area", cut=0, size=4, aspect=2, ax=axes[0])
    #sns.factorplot(x="time", y="H2BGFP intensity", data=my_df[my_df.experiment == 'B'], hue = "animal", kind="violin", color="green", split=False, inner="quartile", scale="area", cut=0, size=4, aspect=2, ax=axes[1])
    #sns.factorplot(x="time", y="H2BGFP intensity", data=my_df[my_df.experiment == 'C'], hue = "animal", kind="violin", color="green", split=False, inner="quartile", scale="area", cut=0, size=4, aspect=2, ax=axes[2])
    #ax.set_xticklabels(rotation=45)
plt.ylim((-12, 2))

# Plot the whole data set (separated per field of view):
fig = plt.figure()
ax = plt.gca()
ax = sns.factorplot(x="field", y="H2BGFP intensity", col="animal", data=my_df, kind="violin", color="grey", split=False, inner="quartile", scale="area", size=4, aspect=.7, cut =0)
plt.ylim((-12, 2))

# Detailed H2BGFP distributions (separated per field of view) in the only animal where the monomodality test yielded sig. (case where "Scale" data set is selected):
if selectDataSet == 5:
    fig = plt.figure()
    ax = plt.gca()
    ax = sns.factorplot(x="field", y="H2BGFP intensity", col="animal", data=my_df[my_df.animal == "18d.1"], kind="violin", color="grey", split=False, inner="quartile", scale="area", size=4, aspect=1.5, cut =0)
    plt.ylim((-5, 1))

