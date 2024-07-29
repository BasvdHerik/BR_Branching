from openalea.lpy import *
import matplotlib; matplotlib.use("TkAgg");import matplotlib.pyplot as plt;
import numpy as np;import random as ran;
import pandas as pd;import seaborn as sns;


exp_data = pd.read_excel('~/Projects/BR_Branching/Model_3/Tomato_graft_data_all.xlsx')

def Lsystem_calc(runs):
    sim = 'Model3.lpy';
    #Division Rate
    # C_part = 0.75;
    # bas1= (0.88*ran.uniform(0.2, 0.25 )*C_part) + 0.88*(ran.uniform(0.2, 0.25) *(1-C_part))      #WT/WT
    # bas2= (0.88*ran.uniform(0.2, 0.25 )*C_part) + 0.88*(ran.uniform(0.1, 0.125)*(1-C_part))      #WT/bri
    # bas3= (0.88*ran.uniform(0.1, 0.125)*C_part) + 0.88*(ran.uniform(0.2, 0.25) *(1-C_part))      #bri/WT
    # bas4= (0.88*ran.uniform(0.1, 0.125)*C_part) + 0.88*(ran.uniform(0.1, 0.125)*(1-C_part))      #bri/bri
    #
    # #meristem size
    # ml1=ran.randint(30,35)*1.75 #WT - root
    # ml2=ran.randint(20,25)*1.75 #bri - root
    #
    # em_age1 = 48        #WT - shoot
    # em_age2 = 120;      #bri - shoot


    C_part = 0.75;
    bas1= (ran.uniform(0.2, 0.25 )*0.3*C_part) + (ran.uniform(0.2, 0.25) *(1-C_part)*0.3) * 0.5     #WT/WT
    bas2= (ran.uniform(0.2, 0.25 )*0.3*C_part) + (ran.uniform(0.1, 0.125)*(1-C_part)*0.3) * 0.5      #WT/bri
    bas3= (ran.uniform(0.1, 0.125)*0.3*C_part) + (ran.uniform(0.2, 0.25) *(1-C_part)*0.3) * 0.5      #bri/WT
    bas4= (ran.uniform(0.1, 0.125)*0.3*C_part) + (ran.uniform(0.1, 0.125)*(1-C_part)*0.3) * 0.5      #bri/bri
    #meristem size
    ml1=ran.randint(150,175)
    ml2=ran.randint(100,125)


    em_age1 = 48
    em_age2 = 120
    # l = Lsystem(sim,{'priming': 0.25, 'maxOrder' :2,  'em' :1/2, 'em_age1': em_age1, 'em_age2': em_age2, 'g_l' :0.66, 'basal1': bas1, 'basal2':bas2, 'ml1':ml1, 'ml2':ml2}); #, 'maxMeristemSize': ml})
    l = Lsystem(sim,{ 'em' :1/2, 'em_age1': em_age1, 'em_age2': em_age2, 'g_l' :0.8, 'basal1': bas1, 'basal2':bas2, 'ml1':ml1, 'ml2':ml2, 'basal3':bas3, 'basal4':bas4}); #, 'maxMeristemSize': ml})


    #l = Lsystem(sim,{'priming': 0.2, 'maxOrder' :2,  'em' :1/2, 'em_age1': em_age1, 'em_age2': em_age2, 'g_l' :0.75, 'basal1': bas1, 'basal2':bas2, 'basal3':bas3, 'basal4':bas4, 'ml1':ml1, 'ml2':ml2}); #, 'maxMeristemSize': ml})
    lstring = l.derive()
    try :
        lat1_l = l.Roots1.secondary_length.item(); lat2_l = l.Roots2.secondary_length.item();lat3_l = l.Roots3.secondary_length.item();lat4_l = l.Roots4.secondary_length.item()
    except AttributeError:
        lat1_l = 0; lat2_l=0; lat3_l=0;lat4_l=0; print('No lats formed')

    new_row_wt_wt    = [l.FW1, bas1, runs, 'WT/WT',   'Simulation' , l.LRn1 , l.TL1 , l.Roots1.primary_length.item() ,    lat1_l, l.BI1, l.LRn1_d7 , l.TL1_d7 , l.PR1_d7 ,  l.LR1_d7, l.BI1_d7, l.FW1_d7,l.LRn1_d0 , l.TL1_d0 , l.PR1_d0 ,  l.LR1_d0, l.BI1_d0, l.FW1_d0]
    new_row_wt_cpd   = [l.FW2, bas1, runs, 'WT/bri',  'Simulation' , l.LRn2 , l.TL2 , l.Roots2.primary_length.item() ,    lat2_l, l.BI2, l.LRn2_d7 , l.TL2_d7 , l.PR2_d7 ,  l.LR2_d7, l.BI2_d7, l.FW2_d7,l.LRn2_d0 , l.TL2_d0 , l.PR2_d0 ,  l.LR2_d0, l.BI2_d0, l.FW2_d0 ]
    new_row_cpd_wt   = [l.FW3, bas2, runs, 'bri/WT', 'Simulation' , l.LRn3 , l.TL3 , l.Roots3.primary_length.item() ,    lat3_l, l.BI3 , l.LRn3_d7 , l.TL3_d7 , l.PR3_d7 ,  l.LR3_d7, l.BI3_d7, l.FW3_d7,l.LRn3_d0 , l.TL3_d0 , l.PR3_d0 ,  l.LR3_d0, l.BI3_d0, l.FW3_d0]
    new_row_cpd_cpd  = [l.FW4, bas2, runs, 'bri/bri', 'Simulation' , l.LRn4 , l.TL4 , l.Roots4.primary_length.item() ,    lat4_l, l.BI4, l.LRn4_d7 , l.TL4_d7 , l.PR4_d7 ,  l.LR4_d7, l.BI4_d7, l.FW4_d7,l.LRn4_d0 , l.TL4_d0 , l.PR4_d0 ,  l.LR4_d0, l.BI4_d0, l.FW4_d0]


    Data_multi = pd.DataFrame(columns = ['FW','TAdiv','ID', 'Genotype', 'Type','LRn','TL','PR','LR','BI','LRn_d7','TL_d7','PR_d7','LR_d7','BI_d7','FW_d7','LRn_d0','TL_d0','PR_d0','LR_d0','BI_d0','FW_d0']);

    Data_multi.loc[len(Data_multi)] =new_row_wt_wt;
    Data_multi.loc[len(Data_multi)] =new_row_wt_cpd;
    Data_multi.loc[len(Data_multi)] =new_row_cpd_wt;
    Data_multi.loc[len(Data_multi)] =new_row_cpd_cpd;

    return Data_multi

import multiprocessing
runs = np.arange(250)
with multiprocessing.Pool() as pool:
    Data_multi = pool.map(Lsystem_calc, runs)
    Data_runs = pd.concat(Data_multi)

Data = Data_runs;
Data.reindex();

path = 'Output/Graft_Fit_new__';


Data_plot = pd.concat([exp_data,Data],ignore_index=True)
Data_plot.to_csv(path + '.csv')

#Plots
colors = ['Blue', 'Green', 'Purple','Orange']; label  = ['Data','Turing','Clock','Reflux']
ylabel = ['Lateral Root Number','Total Root Length (cm)','Primary Root Length (cm)','Lateral Root Length (cm)', 'Branching Index', 'Priming Site Production (h-1)', 'Priming Site Spacing (cm)']
ylim   = [100, 50, 7];
ylabel = ['Lateral Root Number','Total Root Length (cm)','Branching Index']
ydata  = ['LRn','TL','BI']

fig, axis = plt.subplots(nrows=3, ncols =6, figsize=(16,8))
sns.boxplot(data=Data_plot, x = 'Genotype', y='FW_d0', hue='Type', ax=axis[0,5])
sns.boxplot(data=Data_plot, x = 'Genotype', y='TL_d0', hue='Type', ax=axis[0,2])
sns.boxplot(data=Data_plot, x = 'Genotype', y='BI_d0', hue='Type', ax=axis[0,4])
sns.boxplot(data=Data_plot, x = 'Genotype', y='LRn_d0',hue='Type', ax=axis[0,3])
sns.boxplot(data=Data_plot, x = 'Genotype', y='PR_d0', hue='Type', ax=axis[0,0])
sns.boxplot(data=Data_plot, x = 'Genotype', y='LR_d0', hue='Type', ax=axis[0,1])

sns.boxplot(data=Data_plot, x = 'Genotype', y='FW_d7', hue='Type', ax=axis[1,5])
sns.boxplot(data=Data_plot, x = 'Genotype', y='TL_d7', hue='Type', ax=axis[1,2])
sns.boxplot(data=Data_plot, x = 'Genotype', y='BI_d7', hue='Type', ax=axis[1,4])
sns.boxplot(data=Data_plot, x = 'Genotype', y='LRn_d7',hue='Type', ax=axis[1,3])
sns.boxplot(data=Data_plot, x = 'Genotype', y='PR_d7', hue='Type', ax=axis[1,0])
sns.boxplot(data=Data_plot, x = 'Genotype', y='LR_d7', hue='Type', ax=axis[1,1])

sns.boxplot(data=Data_plot, x = 'Genotype', y='FW', hue='Type', ax=axis[2,5])
sns.boxplot(data=Data_plot, x = 'Genotype', y='TL', hue='Type', ax=axis[2,2])
sns.boxplot(data=Data_plot, x = 'Genotype', y='BI', hue='Type', ax=axis[2,4])
sns.boxplot(data=Data_plot, x = 'Genotype', y='LRn',hue='Type', ax=axis[2,3])
sns.boxplot(data=Data_plot, x = 'Genotype', y='PR', hue='Type', ax=axis[2,0])
sns.boxplot(data=Data_plot, x = 'Genotype', y='LR', hue='Type', ax=axis[2,1])

# axis[0,0].set_ylim(0,20);axis[1,0].set_ylim(0,20);axis[2,0].set_ylim(0,20)
# axis[0,1].set_ylim(0,100);axis[1,1].set_ylim(0,100);axis[2,1].set_ylim(0,100)
# axis[0,2].set_ylim(0,100);axis[1,2].set_ylim(0,100);axis[2,2].set_ylim(0,100)
# axis[0,3].set_ylim(0,50);axis[1,3].set_ylim(0,50);axis[2,3].set_ylim(0,50)
# axis[0,4].set_ylim(0,2.5);axis[1,4].set_ylim(0,2.5);axis[2,4].set_ylim(0,2.5)
# axis[0,5].set_ylim(0,100);axis[1,5].set_ylim(0,100);axis[2,5].set_ylim(0,100)
plt.tight_layout()

# fig, axis = plt.subplots(nrows=1, ncols =6, figsize=(12,3))
# sns.boxplot(data=Data_plot, x = 'Genotype', y='FW', hue='Type', ax=axis[5],legend=False)
# sns.boxplot(data=Data_plot, x = 'Genotype', y='TL', hue='Type', ax=axis[2],legend=False)
# sns.boxplot(data=Data_plot, x = 'Genotype', y='BI', hue='Type', ax=axis[4],legend=False)
# sns.boxplot(data=Data_plot, x = 'Genotype', y='LRn',hue='Type', ax=axis[3],legend=False)
# sns.boxplot(data=Data_plot, x = 'Genotype', y='PR', hue='Type', ax=axis[0],legend=False)
# sns.boxplot(data=Data_plot, x = 'Genotype', y='LR', hue='Type', ax=axis[1],legend=False)
#
# axis[0].set_ylim(0,20);
# axis[1].set_ylim(0,100);
# axis[2].set_ylim(0,100);
# axis[3].set_ylim(0,50);
# axis[4].set_ylim(0,2.5);
# axis[5].set_ylim(0,100);
#
# plt.tight_layout()
#
plt.savefig(path + '_overview.png')
plt.savefig(path + '_overview.svg')

#
# fig, axis = plt.subplots(figsize=(4,3))
# sns.boxplot(data=Data_plot, x = 'Genotype', y='BI', hue='Type')
# axis.set_ylabel('Branching Index (LR/cm)')
# plt.tight_layout()
# plt.savefig(path + '_BI.png')
# plt.savefig(path + '_BI.svg')
#
# fig, axis = plt.subplots(figsize=(4,3))
# sns.boxplot(data=Data_plot, x = 'Genotype', y='FW', hue='Type')
# axis.set_ylabel('Root Fresh Weight (mg)')
# plt.tight_layout()
# plt.savefig(path + '_FW.png')
# plt.savefig(path + '_FW.svg')
#
# fig, axis = plt.subplots(figsize=(4,3))
# sns.boxplot(data=Data_plot, x = 'Genotype', y='TL', hue='Type')
# axis.set_ylabel('Total Length (cm)')
# plt.tight_layout()
# plt.savefig(path + '_TL.png')
# plt.savefig(path + '_TL.svg')
#
# fig, axis = plt.subplots(figsize=(4,3))
# sns.boxplot(data=Data_plot, x = 'Genotype', y='LRn', hue='Type')
# axis.set_ylabel('Lateral Root Number')
# plt.tight_layout()
# plt.savefig(path + '_LRn.png')
# plt.savefig(path + '_LRn.svg')
#
#
# fig, axis = plt.subplots(figsize=(4,3))
# sns.boxplot(data=Data_plot, x = 'Genotype', y='PR', hue='Type')
# axis.set_ylabel('Primary Root Length (cm)')
# plt.tight_layout()
# plt.savefig(path + '_PR.png')
# plt.savefig(path + '_PR.svg')
#
# fig, axis = plt.subplots(figsize=(4,3))
# sns.boxplot(data=Data_plot, x = 'Genotype', y='LR', hue='Type')
# axis.set_ylabel('Lateral Root Length (cm)')
# plt.tight_layout()
# plt.savefig(path + '_LR.png')
# plt.savefig(path + '_LR.svg')

# fig, axis = plt.subplots(nrows=1, ncols =5, figsize=(16,8))
# sns.boxplot(data=Data_plot, x = 'Genotype', y='FW', hue='Type', ax=axis[0])
# sns.boxplot(data=Data_plot, x = 'Genotype', y='TL', hue='Type', ax=axis[1])
# sns.boxplot(data=Data_plot, x = 'Genotype', y='BI', hue='Type', ax=axis[2])
# sns.boxplot(data=Data_plot, x = 'Genotype', y='LRn',hue='Type', ax=axis[3])
# sns.boxplot(data=Data_plot, x = 'Genotype', y='PR', hue='Type', ax=axis[4])
# # plt.tight_layout()
DataWT_WT   =  Data[(Data['Genotype'] == 'WT/WT') ]
DataWT_bri  =  Data[(Data['Genotype'] == 'WT/bri') ]
Databri_WT  =  Data[(Data['Genotype'] == 'bri/WT') ]
Databri_bri =  Data[(Data['Genotype'] == 'bri/bri') ]
fig, ax = plt.subplots(ncols=3, figsize=(12,3))
geno=['WT/WT', 'WT/WT' ,'WT/WT' ,'WT/bri','WT/bri','WT/bri','bri/bri','bri/bri','bri/bri', 'bri/WT','bri/WT','bri/WT']
t= [7,14,21,7,14,21,7,14,21,7,14,21]
# PRwtwt =[np.mean(np.array(DataWT_WT['PR_d0']), axis=0 )/np.mean(np.array(DataWT_WT['PR_d0']), axis=0 ),np.mean(np.array(DataWT_WT['PR_d7']), axis=0 )/np.mean(np.array(DataWT_WT['PR_d7']), axis=0 ),np.mean(np.array(DataWT_WT['PR']), axis=0 )/np.mean(np.array(DataWT_WT['PR']), axis=0 )]
# PRwtbri =[np.mean(np.array(DataWT_bri['PR_d0']), axis=0 )/np.mean(np.array(DataWT_WT['PR_d0']), axis=0 ),np.mean(np.array(DataWT_bri['PR_d7'])/np.mean(np.array(DataWT_WT['PR_d7']), axis=0 ), axis=0 ),np.mean(np.array(DataWT_bri['PR']), axis=0 )/np.mean(np.array(DataWT_WT['PR']), axis=0 )]
# PRbriwt = [np.mean(np.array(Databri_WT['PR_d0']), axis=0 )/np.mean(np.array(DataWT_WT['PR_d0']), axis=0 ),np.mean(np.array(Databri_WT['PR_d7'])/np.mean(np.array(DataWT_WT['PR_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_WT['PR']), axis=0 )/np.mean(np.array(DataWT_WT['PR']), axis=0 )]
# PRbribri = [np.mean(np.array(Databri_bri['PR_d0']), axis=0 )/np.mean(np.array(DataWT_WT['PR_d0']), axis=0 ),np.mean(np.array(Databri_bri['PR_d7'])/np.mean(np.array(DataWT_WT['PR_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_bri['PR']), axis=0 )/np.mean(np.array(DataWT_WT['PR']), axis=0 )]
# PR = np.concatenate((PRwtwt,PRwtbri,PRbriwt,PRbribri))
# PRHeatmap = pd.DataFrame({'Genotype' : ['WT/WT', 'WT/WT' ,'WT/WT' ,'WT/bri','WT/bri','WT/bri', 'bri/WT','bri/WT','bri/WT','bri/bri','bri/bri','bri/bri'], 'Time': [0,7,14,0,7,14,0,7,14,0,7,14] , 'Value' : PR})
# PRHeatmap = PRHeatmap.pivot(index='Genotype', columns='Time', values='Value')
# sns.heatmap(PRHeatmap,cmap='crest_r',annot=True,ax = ax[0], vmin = 0 , vmax =1)
# ax[0].set_title('Primary Root Length')
LRNwtwt =[np.mean(np.array(DataWT_WT['LRn_d0']), axis=0 )/np.mean(np.array(DataWT_WT['LRn_d0']), axis=0 ),np.mean(np.array(DataWT_WT['LRn_d7']), axis=0 )/np.mean(np.array(DataWT_WT['LRn_d7']), axis=0 ),np.mean(np.array(DataWT_WT['LRn']), axis=0 )/np.mean(np.array(DataWT_WT['LRn']), axis=0 )]
LRNwtbri =[np.mean(np.array(DataWT_bri['LRn_d0']), axis=0 )/np.mean(np.array(DataWT_WT['LRn_d0']), axis=0 ),np.mean(np.array(DataWT_bri['LRn_d7'])/np.mean(np.array(DataWT_WT['LRn_d7']), axis=0 ), axis=0 ),np.mean(np.array(DataWT_bri['LRn']), axis=0 )/np.mean(np.array(DataWT_WT['LRn']), axis=0 )]
LRNbriwt = [np.mean(np.array(Databri_WT['LRn_d0']), axis=0 )/np.mean(np.array(DataWT_WT['LRn_d0']), axis=0 ),np.mean(np.array(Databri_WT['LRn_d7'])/np.mean(np.array(DataWT_WT['LRn_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_WT['LRn']), axis=0 )/np.mean(np.array(DataWT_WT['LRn']), axis=0 )]
LRNbribri = [np.mean(np.array(Databri_bri['LRn_d0']), axis=0 )/np.mean(np.array(DataWT_WT['LRn_d0']), axis=0 ),np.mean(np.array(Databri_bri['LRn_d7'])/np.mean(np.array(DataWT_WT['LRn_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_bri['LRn']), axis=0 )/np.mean(np.array(DataWT_WT['LRn']), axis=0 )]
LRn = np.concatenate((LRNwtwt,LRNwtbri,LRNbribri,LRNbriwt))
LRnHeatmap = pd.DataFrame({'Genotype' : geno, 'Time': t , 'Value' : LRn})
LRnHeatmap = LRnHeatmap.pivot(index='Genotype', columns='Time', values='Value')
sns.heatmap(LRnHeatmap,cmap='crest_r',annot=True,ax = ax[0], vmin = 0 , vmax =1)
ax[0].set_title('Lateral Root Number')

TLwtwt =[np.mean(np.array(DataWT_WT['TL_d0']), axis=0 )/np.mean(np.array(DataWT_WT['TL_d0']), axis=0 ),np.mean(np.array(DataWT_WT['TL_d7']), axis=0 )/np.mean(np.array(DataWT_WT['TL_d7']), axis=0 ),np.mean(np.array(DataWT_WT['TL']), axis=0 )/np.mean(np.array(DataWT_WT['TL']), axis=0 )]
TLwtbri =[np.mean(np.array(DataWT_bri['TL_d0']), axis=0 )/np.mean(np.array(DataWT_WT['TL_d0']), axis=0 ),np.mean(np.array(DataWT_bri['TL_d7'])/np.mean(np.array(DataWT_WT['TL_d7']), axis=0 ), axis=0 ),np.mean(np.array(DataWT_bri['TL']), axis=0 )/np.mean(np.array(DataWT_WT['TL']), axis=0 )]
TLbriwt = [np.mean(np.array(Databri_WT['TL_d0']), axis=0 )/np.mean(np.array(DataWT_WT['TL_d0']), axis=0 ),np.mean(np.array(Databri_WT['TL_d7'])/np.mean(np.array(DataWT_WT['TL_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_WT['TL']), axis=0 )/np.mean(np.array(DataWT_WT['TL']), axis=0 )]
TLbribri = [np.mean(np.array(Databri_bri['TL_d0']), axis=0 )/np.mean(np.array(DataWT_WT['TL_d0']), axis=0 ),np.mean(np.array(Databri_bri['TL_d7'])/np.mean(np.array(DataWT_WT['TL_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_bri['TL']), axis=0 )/np.mean(np.array(DataWT_WT['TL']), axis=0 )]

# TLwtwt =[np.mean(np.array(DataWT_WT['LR_d0']), axis=0 )/np.mean(np.array(DataWT_WT['LR_d0']), axis=0 ),np.mean(np.array(DataWT_WT['LR_d7']), axis=0 )/np.mean(np.array(DataWT_WT['LR_d7']), axis=0 ),np.mean(np.array(DataWT_WT['LR']), axis=0 )/np.mean(np.array(DataWT_WT['LR']), axis=0 )]
# TLwtbri =[np.mean(np.array(DataWT_bri['LR_d0']), axis=0 )/np.mean(np.array(DataWT_WT['LR_d0']), axis=0 ),np.mean(np.array(DataWT_bri['LR_d7'])/np.mean(np.array(DataWT_WT['LR_d7']), axis=0 ), axis=0 ),np.mean(np.array(DataWT_bri['LR']), axis=0 )/np.mean(np.array(DataWT_WT['LR']), axis=0 )]
# TLbriwt = [np.mean(np.array(Databri_WT['LR_d0']), axis=0 )/np.mean(np.array(DataWT_WT['LR_d0']), axis=0 ),np.mean(np.array(Databri_WT['LR_d7'])/np.mean(np.array(DataWT_WT['LR_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_WT['LR']), axis=0 )/np.mean(np.array(DataWT_WT['LR']), axis=0 )]
# TLbribri = [np.mean(np.array(Databri_bri['LR_d0']), axis=0 )/np.mean(np.array(DataWT_WT['LR_d0']), axis=0 ),np.mean(np.array(Databri_bri['LR_d7'])/np.mean(np.array(DataWT_WT['LR_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_bri['LR']), axis=0 )/np.mean(np.array(DataWT_WT['LR']), axis=0 )]

TL = np.concatenate((TLwtwt,TLwtbri,TLbribri,TLbriwt))
PRHeatmap = pd.DataFrame({'Genotype' : geno, 'Time': t , 'Value' : TL})
PRHeatmap = PRHeatmap.pivot(index='Genotype', columns='Time', values='Value')
sns.heatmap(PRHeatmap,cmap='crest_r',annot=True,ax = ax[1], vmin = 0 , vmax =1)
ax[1].set_title('Total Root Length')

TLwtwt =[np.mean(np.array(DataWT_WT['FW_d0']), axis=0 )/np.mean(np.array(DataWT_WT['FW_d0']), axis=0 ),np.mean(np.array(DataWT_WT['FW_d7']), axis=0 )/np.mean(np.array(DataWT_WT['FW_d7']), axis=0 ),np.mean(np.array(DataWT_WT['FW']), axis=0 )/np.mean(np.array(DataWT_WT['FW']), axis=0 )]
TLwtbri =[np.mean(np.array(DataWT_bri['FW_d0']), axis=0 )/np.mean(np.array(DataWT_WT['FW_d0']), axis=0 ),np.mean(np.array(DataWT_bri['FW_d7'])/np.mean(np.array(DataWT_WT['FW_d7']), axis=0 ), axis=0 ),np.mean(np.array(DataWT_bri['FW']), axis=0 )/np.mean(np.array(DataWT_WT['FW']), axis=0 )]
TLbriwt = [np.mean(np.array(Databri_WT['FW_d0']), axis=0 )/np.mean(np.array(DataWT_WT['FW_d0']), axis=0 ),np.mean(np.array(Databri_WT['FW_d7'])/np.mean(np.array(DataWT_WT['FW_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_WT['FW']), axis=0 )/np.mean(np.array(DataWT_WT['FW']), axis=0 )]
TLbribri = [np.mean(np.array(Databri_bri['FW_d0']), axis=0 )/np.mean(np.array(DataWT_WT['FW_d0']), axis=0 ),np.mean(np.array(Databri_bri['FW_d7'])/np.mean(np.array(DataWT_WT['FW_d7']), axis=0 ), axis=0 ),np.mean(np.array(Databri_bri['FW']), axis=0 )/np.mean(np.array(DataWT_WT['FW']), axis=0 )]
TL = np.concatenate((TLwtwt,TLwtbri,TLbribri,TLbriwt))
TLHeatmap = pd.DataFrame({'Genotype' : geno, 'Time': t , 'Value' : TL})
TLHeatmap = TLHeatmap.pivot(index='Genotype', columns='Time', values='Value')
sns.heatmap(TLHeatmap,cmap='crest_r',annot=True,ax = ax[2], vmin = 0 , vmax =1)
ax[2].set_title('Root Fresh Weight')

plt.tight_layout()
plt.savefig(path + '_Heatmap.png')
plt.savefig(path + '_Heatmap.svg')

# fig, ax = plt.subplots(ncols=3, figsize=(12,3))
#
# LRn =[ 1, 1, 1, 0.33,0.55,0.98,0.33,0.1,0.15,0.67,0.35,0.44]
# LRnHeatmap = pd.DataFrame({'Genotype' : geno, 'Time':t , 'Value' : LRn})
# LRnHeatmap = LRnHeatmap.pivot(index='Genotype', columns='Time', values='Value')
# sns.heatmap(LRnHeatmap,cmap='crest_r',annot=True,ax = ax[0], vmin = 0 , vmax =1)
# ax[0].set_title('Lateral Root Number')
#
# PR = [1,1,1,0.92,0.92,0.92,0.38,0.38,0.38,0.31,0.31,0.31]
# PRHeatmap = pd.DataFrame({'Genotype' :geno, 'Time': t , 'Value' : PR})
# PRHeatmap = PRHeatmap.pivot(index='Genotype', columns='Time', values='Value')
# sns.heatmap(PRHeatmap,cmap='crest_r',annot=True,ax = ax[2], vmin = 0 , vmax =1)
# ax[2].set_title('Root Fresh Weight')
#
# TL = [1,1,1,0.55,0.69,0.75,0.54,0.38,0.23,0.97,0.69,0.37]
# TLHeatmap = pd.DataFrame({'Genotype' :geno, 'Time': t , 'Value' : TL})
# TLHeatmap = TLHeatmap.pivot(index='Genotype', columns='Time', values='Value')
# sns.heatmap(TLHeatmap,cmap='crest_r',annot=True,ax = ax[1], vmin = 0 , vmax =1)
# ax[1].set_title('Total Root Length')
#
#
# # TL = [1,1,1,0.27,0.52,0.55,0.2,0.04,0.04,0.92,0.20,0.14]
# # TLHeatmap = pd.DataFrame({'Genotype' :geno, 'Time': t , 'Value' : TL})
# # TLHeatmap = TLHeatmap.pivot(index='Genotype', columns='Time', values='Value')
# # sns.heatmap(TLHeatmap,cmap='crest_r',annot=True,ax = ax[1], vmin = 0 , vmax =1)
# # ax[1].set_title('Lateral Root Length')
#
# plt.tight_layout()
# plt.savefig(path + '_Heatmap_Data.png')
# plt.savefig(path + '_Heatmap_Data.svg')
plt.show()
