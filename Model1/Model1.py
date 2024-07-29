from openalea.lpy import *
import matplotlib; matplotlib.use("TkAgg");import matplotlib.pyplot as plt;
import numpy as np;import random as ran;
import pandas as pd;import seaborn as sns;

def Lsystem_calc(runs):
    sim = 'Model1.lpy';

    #Division Rate
    bas1 = ran.uniform(0.2, 0.25);    bas2= ran.uniform(0.1, 0.125)
    ml1= ran.randint(30,35);    ml2= ran.randint(20,25)

    l = Lsystem(sim,{  'em' :1/6, 'em_age': 48, 'g_l' :0.6, 'basal1': bas1, 'basal2':bas2, 'ml1':ml1, 'ml2':ml2}); #, 'maxMeristemSize': ml})
    lstring = l.derive()
    try :
        lat1_l = l.Roots1.secondary_length.item(); lat2_l = l.Roots2.secondary_length.item();#lat3_l = l.Roots3.secondary_length.item()
    except AttributeError:
        lat1_l = 0; lat2_l=0;  print('No lats formed');

    new_row_wt   = ['WT',   'Simulation' ,l.LRn1 , l.TL1 , l.Roots1.primary_length.item() ,    lat1_l, l.BI1 ]
    new_row_mut    = ['Mutant',    'Simulation' ,l.LRn2 , l.TL2 , l.Roots2.primary_length.item() ,    lat2_l, l.BI2 ]
    Data_multi = pd.DataFrame(columns = ['Genotype', 'Type','LRn','TL','PR','LR','BI']);
    Data_multi.loc[len(Data_multi)] =new_row_wt;
    Data_multi.loc[len(Data_multi)] =new_row_mut;

    return Data_multi

import multiprocessing
runs = np.arange(250)
with multiprocessing.Pool() as pool:
    Data_multi = pool.map(Lsystem_calc, runs)
    Data_runs = pd.concat(Data_multi)

Data = Data_runs;
Data.reindex();

exp_data = pd.read_excel('Arabidopsis_Data.xlsx')

path = 'test2_'
Data_plot = pd.concat([exp_data,Data],ignore_index=True)

#Save sim results
Data_plot.to_csv(path + '.csv')

#Plot Results
fig, axis = plt.subplots(figsize=(4,3))
sns.boxplot(data=Data_plot, x = 'Type', y='BI', hue='Genotype')
axis.set_ylabel('Branching Index (LR/cm)')
plt.tight_layout()
plt.savefig(path + '_BI.png')
plt.savefig(path + '_BI.svg')

fig, axis = plt.subplots(figsize=(4,3))
sns.boxplot(data=Data_plot, x = 'Type', y='TL', hue='Genotype')
axis.set_ylabel('Total Root Length (cm)')
plt.tight_layout()
plt.savefig(path + '_TL.png')
plt.savefig(path + '_TL.svg')

fig, axis = plt.subplots(figsize=(4,3))
sns.boxplot(data=Data_plot, x = 'Type', y='LRn', hue='Genotype')
axis.set_ylabel('Lateral Root Number')
plt.tight_layout()
plt.savefig(path + '_LRn.png')
plt.savefig(path + '_LRn.svg')

fig, axis = plt.subplots(figsize=(4,3))
sns.boxplot(data=Data_plot, x = 'Genotype', y='PR')
axis.set_ylabel('Primary Root Length (cm)')
plt.tight_layout()
plt.savefig(path + '_PR.png')
plt.savefig(path + '_PR.svg')

fig, axis = plt.subplots(figsize=(4,3))
sns.boxplot(data=Data_plot, x = 'Genotype', y='LR')
axis.set_ylabel('Lateral Root Length (cm)')
plt.tight_layout()
plt.savefig(path + '_LR.png')
plt.savefig(path + '_LR.svg')

plt.show()
