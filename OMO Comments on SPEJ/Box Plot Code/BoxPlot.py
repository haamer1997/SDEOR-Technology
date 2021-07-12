# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 23:20:05 2021

@author: ahass16
"""
#%% Import libraries
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
#%% Plot Conductive
#dataset = pd.read_excel('BoxPlot.xlsx', sheet_name='Conductive')
#fig_dims = (15, 8)
#fig, ax = plt.subplots(figsize=fig_dims)
#
#g = sns.boxplot(x="NFs", y="RF", data=dataset, ax = ax)
#
#ax.tick_params(labelsize=22)
#ax.figure.tight_layout()
#
#ax.set_ylabel('Oil EOR RF (%)',fontsize=38,fontweight='bold',fontname="Times New Roman Bold")
#ax.set_xlabel('Number of Natural Fractures',fontsize=34,fontweight='bold',fontname="Times New Roman Bold")
#plt.setp(ax.get_xticklabels(), fontsize=34,fontname="Times New Roman")
#plt.setp(ax.get_yticklabels(), fontsize=34,fontname="Times New Roman")
#ax.set_axisbelow(True)
#ax.grid(linestyle='--', linewidth=0.5)
#plt.savefig('ConductiveNFs.png',dpi=300,bbox_inches='tight')
#%% Plot Non-Conductive
#dataset = pd.read_excel('BoxPlot.xlsx', sheet_name='Non-Conductive')
#fig_dims = (15, 8)
#fig, ax = plt.subplots(figsize=fig_dims)
#
#g = sns.boxplot(x="NFs", y="RF", data=dataset, ax = ax)
#
#ax.tick_params(labelsize=22)
#ax.figure.tight_layout()
#
#ax.set_ylabel('Oil EOR RF (%)',fontsize=38,fontweight='bold',fontname="Times New Roman Bold")
#ax.set_xlabel('Number of Natural Fractures',fontsize=34,fontweight='bold',fontname="Times New Roman Bold")
#plt.setp(ax.get_xticklabels(), fontsize=34,fontname="Times New Roman")
#plt.setp(ax.get_yticklabels(), fontsize=34,fontname="Times New Roman")
#ax.set_axisbelow(True) #plot grids behind graph elements
#ax.grid(linestyle='--', linewidth=0.5)
#plt.savefig('Non_ConductiveNFs.png',dpi=300,bbox_inches='tight')
#%% Plot Non-Conductive
dataset = pd.read_excel('BoxPlot.xlsx', sheet_name='Mixed_Conductive')
fig_dims = (13, 6)
fig, ax = plt.subplots(figsize=fig_dims)

g = sns.boxplot(x="Conductivity", y="RF", data=dataset, ax = ax)

ax.tick_params(labelsize=22)
ax.figure.tight_layout()

ax.set_ylabel('Oil EOR RF (%)',fontsize=38,fontweight='bold',fontname="Times New Roman Bold")
ax.set_xlabel('Natural Fracture Conductivity',fontsize=34,fontweight='bold',fontname="Times New Roman Bold")
plt.setp(ax.get_xticklabels(), fontsize=34,fontname="Times New Roman")
plt.setp(ax.get_yticklabels(), fontsize=34,fontname="Times New Roman")
ax.set_axisbelow(True) #plot grids behind graph elements
ax.grid(linestyle='--', linewidth=0.5)
#plt.savefig('mixed_conductivity.png',dpi=300,bbox_inches='tight')
plt.savefig('mixed_conductivity.png',dpi=300,bbox_inches='tight')

