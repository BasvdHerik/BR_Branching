from openalea.lpy import *
import matplotlib; matplotlib.use("TkAgg");import matplotlib.pyplot as plt;
import numpy as np;import random as ran;
import pandas as pd;import seaborn as sns;


exp_data = pd.read_excel('Tomato_graft_data_all.xlsx')

def Lsystem_calc(runs):
    sim = 'Model3.lpy';

    C_part = 0.66;
    bas1= ((ran.uniform(0.2, 0.25 )*C_part) + (ran.uniform(0.2, 0.25) *(1-C_part)))    #WT/WT
    bas2= ((ran.uniform(0.2, 0.25 )*C_part) + (ran.uniform(0.1, 0.125)*(1-C_part)))    #WT/bri
    bas3= ((ran.uniform(0.1, 0.125)*C_part) + (ran.uniform(0.2, 0.25) *(1-C_part)))    #bri/WT
    bas4= ((ran.uniform(0.1, 0.125)*C_part) + (ran.uniform(0.1, 0.125)*(1-C_part)))    #bri/bri
    #meristem size
    ml1=ran.randint(110,130)
    ml2=ran.randint(80, 100)


    em_age1 = 84
    em_age2 = 120
    alpha = ran.uniform(0.85, 1.15 )

    l = Lsystem(sim,  { 'stoch': alpha, 'em' :1, 'em_age1': em_age1, 'em_age2': em_age2, 'g_l' :0.4, 'basal1': bas1, 'basal2':bas2, 'ml1':ml1, 'ml2':ml2, 'basal3':bas3, 'basal4':bas4}); #, 'maxMeristemSize': ml})


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

Data = Data_runs;Data.reindex();Data_plot = pd.concat([exp_data,Data],ignore_index=True)

#Plots
fig, axis = plt.subplots(nrows=5, ncols=3, figsize=(8,11.5))
sns.boxplot(data=Data_plot, x='Genotype', y='FW_d0', hue='Type', ax=axis[4,0])
axis[4,0].set_ylabel('Fresh Weight')
sns.boxplot(data=Data_plot, x='Genotype', y='TL_d0', hue='Type', ax=axis[1,0])
axis[1,0].set_ylabel('Total Root Length')
sns.boxplot(data=Data_plot, x='Genotype', y='BI_d0', hue='Type', ax=axis[3,0])
axis[3,0].set_ylabel('Branching Index')
sns.boxplot(data=Data_plot, x='Genotype', y='LRn_d0', hue='Type', ax=axis[2,0])
axis[2,0].set_ylabel('Lateral Root Number')
sns.boxplot(data=Data_plot, x='Genotype', y='LR_d0', hue='Type', ax=axis[0,0])
axis[0,0].set_ylabel('Lateral Root Length')

sns.boxplot(data=Data_plot, x='Genotype', y='FW_d7', hue='Type', ax=axis[4,1])
axis[4,1].set_ylabel('Fresh Weight')
sns.boxplot(data=Data_plot, x='Genotype', y='TL_d7', hue='Type', ax=axis[1,1])
axis[1,1].set_ylabel('Total Root Length')
sns.boxplot(data=Data_plot, x='Genotype', y='BI_d7', hue='Type', ax=axis[3,1])
axis[3,1].set_ylabel('Branching Index')
sns.boxplot(data=Data_plot, x='Genotype', y='LRn_d7', hue='Type', ax=axis[2,1])
axis[2,1].set_ylabel('Lateral Root Number')
sns.boxplot(data=Data_plot, x='Genotype', y='LR_d7', hue='Type', ax=axis[0,1])
axis[0,1].set_ylabel('Lateral Root Length')

sns.boxplot(data=Data_plot, x='Genotype', y='FW', hue='Type', ax=axis[4,2])
axis[4,2].set_ylabel('Fresh Weight')
sns.boxplot(data=Data_plot, x='Genotype', y='TL', hue='Type', ax=axis[1,2])
axis[1,2].set_ylabel('Total Root Length')
sns.boxplot(data=Data_plot, x='Genotype', y='BI', hue='Type', ax=axis[3,2])
axis[3,2].set_ylabel('Branching Index')
sns.boxplot(data=Data_plot, x='Genotype', y='LRn', hue='Type', ax=axis[2,2])
axis[2,2].set_ylabel('Lateral Root Number')
sns.boxplot(data=Data_plot, x='Genotype', y='LR', hue='Type', ax=axis[0,2])
axis[0,2].set_ylabel('Lateral Root Length')
plt.tight_layout()
for ax in axis.flat:
    ax.legend_.remove()

plt.show()
