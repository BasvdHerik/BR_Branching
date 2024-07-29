from openalea.lpy import *
import matplotlib; matplotlib.use("TkAgg");import matplotlib.pyplot as plt;
import numpy as np;import random as ran;
import pandas as pd;import seaborn as sns;


exp_data = pd.read_csv('~/Projects/BR_Branching/Model_2/TMTcomb.csv')

def Lsystem_calc(runs):
    sim = 'Model2.lpy';
    #Division Rate
    bas1= ran.uniform(0.06, 0.075);    bas2= ran.uniform(0.03, 0.0375)
    #meristem size
    ml1=ran.randint(110,130);    ml2=ran.randint(80,100)

    em_age1 = 48;    em_age2 = 120

    l = Lsystem(sim,{ 'em' :1/2, 'em_age1': em_age1, 'em_age2': em_age2, 'g_l' :0.8, 'basal1': bas1, 'basal2':bas2, 'ml1':ml1, 'ml2':ml2}); #, 'maxMeristemSize': ml})
    lstring = l.derive()
    try :
        lat1_l = l.Roots1.secondary_length.item(); lat2_l = l.Roots2.secondary_length.item();
    except AttributeError:
        lat1_l = 0; lat2_l=0; lat3_l=0;lat4_l=0; print('No lats formed')
    new_row_wt_wt    = ['Simulation', bas1 ,  runs, 'WT',     l.LRn1_d10 ,l.LRn1_d20 , l.TL1_d10, l.TL1_d20 , l.BI1_d10,l.BI1_d20,l.LRn1_d14, l.TL1_d14,l.BI1_d14, ml1,l.Root1_store,l.LRn1_store,l.LRn1_store/l.Root1_store,l.time_store]
    new_row_cpd_cpd  = ['Simulation', bas2 ,  runs, 'Mutant', l.LRn2_d10 , l.LRn2_d20, l.TL2_d10, l.TL2_d20 , l.BI2_d10,l.BI2_d20,l.LRn2_d14, l.TL2_d14,l.BI2_d14,ml2,l.Root2_store,l.LRn2_store,l.LRn2_store/l.Root2_store,l.time_store]
    Data_multi = pd.DataFrame(columns = ['Type','TAdiv','ID', 'Genotype', 'LRn_d10','LRn_d20','TL_d10','TL_d20','BI_d10','BI_d20','LRn_d14','TL_d14','BI_d14','Meristem_Size','TL_time', 'LRn_time','BI_time', 'Time']);
    Data_multi.loc[len(Data_multi)] =new_row_wt_wt;
    Data_multi.loc[len(Data_multi)] =new_row_cpd_cpd;
    return Data_multi

import multiprocessing
runs = np.arange(250)
with multiprocessing.Pool() as pool:
    Data_multi = pool.map(Lsystem_calc, runs)
    Data_runs = pd.concat(Data_multi)

path = 'Output/';

#Plots
Data = Data_runs;
Data.reindex();
DataWT =  Data[(Data['Genotype'] == 'WT') ]
DataMut = Data[(Data['Genotype'] == 'Mutant') ]

Time = Data['Time'].mean()/(24)


WT_LRn_m = np.mean( np.array(DataWT['LRn_time']), axis=0 );WT_LRn_s =np.std( np.array(DataWT['LRn_time']), axis=0 )
WT_TL_m  = np.mean( np.array(DataWT['TL_time']), axis=0 ); WT_TL_s   =np.std( np.array(DataWT['TL_time']), axis=0 )
WT_BI_m   =np.mean( np.array(DataWT['BI_time']), axis=0 ); WT_BI_s  = np.std( np.array(DataWT['BI_time']), axis=0 )

Mut_LRn_m = np.mean( np.array(DataMut['LRn_time']), axis=0 ) ;Mut_LRn_s =np.std( np.array(DataMut['LRn_time']), axis=0 )
Mut_TL_m  = np.mean( np.array(DataMut['TL_time']), axis=0 )  ;Mut_TL_s  =np.std( np.array(DataMut['TL_time']), axis=0 )
Mut_BI_m  = np.mean( np.array(DataMut['BI_time']), axis=0 );  Mut_BI_s  =np.std( np.array(DataMut['BI_time']), axis=0 )

fig, axis = plt.subplots(nrows=1, figsize=(4,3))
axis.plot(Time, WT_BI_m,label = 'WT')
axis.fill_between(Time, WT_BI_m-WT_BI_s,WT_BI_m+WT_BI_s,alpha = 0.5)
axis.plot(Time, Mut_BI_m,label = 'Mutant');
axis.fill_between(Time, Mut_BI_m-Mut_BI_s,Mut_BI_m+Mut_BI_s,alpha = 0.5)
axis.set_xlabel('Time (day))'); axis.set_xlim(0,30)
plt.tight_layout();plt.legend()

sns.scatterplot(data=exp_data, x='dpg',y='BI', hue='Genotype', legend=False)
axis.set_ylabel('Branching Index (LR/cm)');
axis.set_ylim(0,2)

plt.savefig(path + 'timeplot_new.png')
plt.savefig(path + 'timeplot_new.svg')

fig, axis = plt.subplots(nrows=1, figsize=(8,6))
axis.plot(WT_TL_m, WT_BI_m,label = 'WT')
axis.fill_between(WT_TL_m, WT_BI_m-WT_BI_s,WT_BI_m+WT_BI_s,alpha = 0.5)
axis.plot(Mut_TL_m, Mut_BI_m,label = 'Mutant');
axis.fill_between(Mut_TL_m, Mut_BI_m-Mut_BI_s,Mut_BI_m+Mut_BI_s,alpha = 0.5)
axis.set_xlabel('Total Length (cm)'); axis.set_xlim(0,200)
plt.tight_layout();plt.legend(['WT','bri'])
sns.scatterplot(data=exp_data, x='Total_root_length',y='BI', hue='Genotype', legend=False, s=50)

predBI = 6.48/np.sqrt(WT_TL_m) - 5.6/WT_TL_m
axis.plot(WT_TL_m,predBI, color = 'grey',alpha = 0.95)

axis.set_ylabel('Branching Index (LR/cm)');
axis.set_ylim(0,2)
plt.savefig(path + 'lengthplot_new.svg')
plt.savefig(path + 'lengthplot_new.png')


exp_data = pd.read_csv('~/Projects/BR_Branching/Model_2/TomatoPlot.csv')
Data_plot = pd.concat([exp_data,Data])
fig, axis = plt.subplots(ncols=3, nrows=3, figsize=(8,6))
sns.boxplot(data=Data_plot, x = 'Type', y='LRn_d10', hue='Genotype', ax=axis[0,0])
sns.boxplot(data=Data_plot, x = 'Type', y='TL_d10', hue='Genotype', ax=axis[0,1])
sns.boxplot(data=Data_plot, x = 'Type', y='BI_d10', hue='Genotype', ax=axis[0,2])

sns.boxplot(data=Data_plot, x = 'Type', y='LRn_d14', hue='Genotype', ax=axis[1,0])
sns.boxplot(data=Data_plot, x = 'Type', y='TL_d14', hue='Genotype', ax=axis[1,1])
sns.boxplot(data=Data_plot, x = 'Type', y='BI_d14', hue='Genotype', ax=axis[1,2])


sns.boxplot(data=Data_plot, x = 'Type', y='LRn_d20', hue='Genotype', ax=axis[2,0])
sns.boxplot(data=Data_plot, x = 'Type', y='TL_d20', hue='Genotype', ax=axis[2,1])
sns.boxplot(data=Data_plot, x = 'Type', y='BI_d20', hue='Genotype', ax=axis[2,2])
plt.tight_layout()
plt.savefig(path + 'boxplot_new.png')
plt.savefig(path + 'boxplot_new.svg')

plt.show()
