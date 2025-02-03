from openalea.lpy import *
import matplotlib; matplotlib.use("TkAgg"); import matplotlib.pyplot as plt
import numpy as np; import random as ran
import pandas as pd; import seaborn as sns
import multiprocessing

def Lsystem_calc(runs):
    sim = 'Model1.lpy'

    #Division Rate & Meristem Sizes
    bas1 = ran.uniform(0.2, 0.25);    bas2= ran.uniform(0.1, 0.125)
    ml1 = ran.randint(30,35);       ml2 = ran.randint(20,25)

    em = 1
    alpha = ran.uniform(0.85, 1.15 )

    l = Lsystem(sim,  { 'stoch': alpha,  'em' :em, 'em_age': 84, 'g_l' :0.66, 'basal1': bas1, 'basal2':bas2, 'ml1':ml1, 'ml2':ml2});
    lstring = l.derive()
    try :
        lat1_l = l.Roots1.secondary_length.item(); lat2_l = l.Roots2.secondary_length.item();
    except AttributeError:
        lat1_l = 0; lat2_l=0;  print('No lats formed');

    new_row_wt   = ['WT',   'Simulation' ,l.LRn1 , l.TL1 , l.Roots1.primary_length.item() ,    lat1_l, l.BI1 , ml1,bas1,em, l.PrimeTime1]
    new_row_mu    = ['Mutant',    'Simulation' ,l.LRn2 , l.TL2 , l.Roots2.primary_length.item() ,    lat2_l, l.BI2 , ml2,bas2,em, l.PrimeTime2]
    Data_multi = pd.DataFrame(columns = ['Genotype', 'Type','LRn','TL','PR','LR','BI','MZsize','DivRate', 'em','prime_dist']);
    Data_multi.loc[len(Data_multi)] =new_row_wt
    Data_multi.loc[len(Data_multi)] =new_row_mu

    return Data_multi

runs = np.arange(250)
with multiprocessing.Pool() as pool:
    Data_multi = pool.map(Lsystem_calc, runs)
    Data_runs = pd.concat(Data_multi)

Data = Data_runs;
Data.reindex();exp_data = pd.read_excel('ArabidopsisData.xlsx')
Data_plot = pd.concat([exp_data, Data], ignore_index=True)


# Plot Results
fig, ax = plt.subplots(1, 1, figsize=(12, 9))
sns.boxplot(data=Data_plot, x='Genotype', y='BI', hue='Genotype')
ax.set_ylabel('Branching Index (LR/cm)')
plt.tight_layout()
ax.set_ylim(0,4)

fig, ax = plt.subplots(2, 2, figsize=(12, 9))
sns.boxplot(data=Data_plot, x='Genotype', y='TL', hue='Genotype', ax=ax[0, 0])
ax[0, 0].set_ylabel('Total Root Length (cm)')
plt.tight_layout()


sns.boxplot(data=Data_plot, x='Genotype', y='LRn', hue='Genotype', ax=ax[0,1])
ax[0,1].set_ylabel('Lateral Root Number')
plt.tight_layout()


sns.boxplot(data=Data_plot, x='Genotype', y='PR', hue='Genotype', ax=ax[1, 0])
ax[1 ,0].set_ylabel('Primary Root Length (cm)')
plt.tight_layout()


sns.boxplot(data=Data_plot, x='Genotype', y='LR', hue='Genotype', ax=ax[1, 1])
ax[1, 1].set_ylabel('Lateral Root Length (cm)')
plt.tight_layout()

plt.show()