from openalea.lpy import *
import matplotlib; matplotlib.use("TkAgg");import matplotlib.pyplot as plt;
import numpy as np;import random as ran;
import pandas as pd;import seaborn as sns;


def Lsystem_calc(runs):
    sim = 'BI_Sigal_Sugar_Radial.lpy';
    if runs < 50: bas2 = 0.075;                          ml2= ran.randint(20,25)
    elif runs >=50 and runs <100:  bas2 = 0.05+0.05;     ml2= ran.randint(20,25)
    elif runs >=100 and runs <150: bas2 = 0.075+0.05;    ml2= ran.randint(20,25)
    elif runs >=150 and runs <200: bas2 = 0.15;          ml2= ran.randint(20,25)
    elif runs >=200 and runs <250: bas2 = 0.125+0.05;    ml2= ran.randint(20,25)
    elif runs >=250 and runs <300: bas2 = 0.150+0.05;    ml2= ran.randint(20,25)

    #Division Rate
    bas1= 0.25;
    # bas1 = ran.uniform(0.2, 0.25);
    # bas2= ran.uniform(0.1, 0.125)
    ml1= ran.randint(30,32)
    ml2= ran.randint(22,25)

    l = Lsystem(sim,{'maxOrder' :3,  'em' :1/6, 'em_age': 36, 'g_l' :0.6, 'basal1': bas1, 'basal2':bas2, 'ml1':ml1, 'ml2':ml2}); #, 'maxMeristemSize': ml})
    lstring = l.derive()
    try :
        lat1_l = l.Roots1.secondary_length.item(); lat2_l = l.Roots2.secondary_length.item();lat3_l = l.Roots3.secondary_length.item()
    except AttributeError:
        lat1_l = 0; lat2_l=0; lat3_l=0; print('No lats formed')

    new_row_wt   = [l.FW1, bas1, runs,'WT',   'Simulation' , bas1,l.LRn1 , l.TL1 , l.Roots1.primary_length.item() ,    lat1_l, l.BI1 ]
    new_row_mut    = [l.FW2,bas2,runs, 'Mutant',    'Simulation' , bas2,l.LRn2 , l.TL2 , l.Roots2.primary_length.item() ,    lat2_l, l.BI2 ]
    Data_multi = pd.DataFrame(columns = ['FW','TAdiv','ID', 'Genotype', 'Type','Scenario','LRn','TL','PR','LR','BI']);
    Data_multi.loc[len(Data_multi)] =new_row_wt;
    Data_multi.loc[len(Data_multi)] =new_row_mut;

    return Data_multi

import multiprocessing
runs = np.arange(300)
with multiprocessing.Pool() as pool:
    Data_multi = pool.map(Lsystem_calc, runs)
    Data_runs = pd.concat(Data_multi)

Data = Data_runs;
# Data.reindex();

exp_data = pd.read_excel('~/Projects/Ara_FSPM/RootRatio/ModelSigal/Ara/AraPlot.xlsx')

path = 'SugarFigure/stable_meristem';
Data.to_csv(path + '.csv')

Data_plot = pd.concat([exp_data,Data])


#Plots
colors = ['Blue', 'Green', 'Purple','Orange']; label  = ['Data','Turing','Clock','Reflux']
ylabel = ['Lateral Root Number','Total Root Length (cm)','Primary Root Length (cm)','Lateral Root Length (cm)', 'Branching Index', 'Priming Site Production (h-1)', 'Priming Site Spacing (cm)']
ylim   = [100, 50, 7];
ylabel = ['Lateral Root Number','Total Root Length (cm)','Branching Index']
ydata  = ['LRn','TL','BI']


fig, axis = plt.subplots(nrows=2, ncols =4, figsize=(16,8))
sns.boxplot(data=Data, x = 'TAdiv', y='LRn' , hue='Genotype', ax=axis[0,0])
sns.boxplot(data=Data, x = 'TAdiv', y='TL' , hue='Genotype', ax=axis[0,1])
sns.boxplot(data=Data, x = 'TAdiv', y='BI' , hue='Genotype', ax=axis[0,2])
sns.boxplot(data=Data, x = 'TAdiv', y='FW' , hue='Genotype', ax=axis[0,3])

sns.boxplot(data=exp_data, x = 'TAdiv', y='LRn', hue='Genotype',  ax=axis[1,0])
sns.boxplot(data=exp_data, x = 'TAdiv', y='TL' ,hue='Genotype',  ax=axis[1,1])
sns.boxplot(data=exp_data, x = 'TAdiv', y='BI' ,hue='Genotype', ax=axis[1,2])
sns.boxplot(data=exp_data, x = 'TAdiv', y='FW' ,hue='Genotype', ax=axis[1,3])

axis[0,0].set_ylim(0,250);axis[1,0].set_ylim(0,250);
axis[0,1].set_ylim(0,125);axis[1,1].set_ylim(0,125);
axis[0,2].set_ylim(0,5);  axis[1,2].set_ylim(0,5);
axis[0,3].set_ylim(0,25); axis[1,3].set_ylim(0,25);

plt.tight_layout()
plt.savefig(path + '_boxplot.png')
plt.savefig(path + '_boxplot.svg')

fig, axis = plt.subplots(figsize=(5,4))
plt.scatter(x,y,s=size, c =color, cmap='Reds' ,edgecolors = 'grey')
plt.savefig(path + '_scatter.png')
plt.savefig(path + '_satter.svg')

fig, axis = plt.subplots(figsize=(5,4))
plt.scatter(Data['TAdiv'],Data['FW'],s=Data['TL'], c =Data['LRn'], cmap='Reds' ,edgecolors = 'grey')
plt.savefig(path + '_scatter_all.png')
plt.savefig(path + '_scatter_all.svg')

Data075 =  exp_data[(exp_data['TAdiv'] == 'CPD_0') ]; print(Data075)
Data100 =  exp_data[(exp_data['TAdiv'] == 'CPD_0_2')];
Data125 =  exp_data[(exp_data['TAdiv'] == 'CPD_0_5')]
Data150 =  exp_data[(exp_data['TAdiv'] == 'CPD_1')]
Data175 =  exp_data[(exp_data['TAdiv'] == 'CPD_1_5')]
Data200 =  exp_data[(exp_data['TAdiv'] == 'CPD_2')]
Data250 =  exp_data[(exp_data['TAdiv'] == 'WT_0')]

x_data = [0.075,0.075,0.100,0.125,0.150,0.175,0.200];
y_data = [np.mean( np.array(Data250['FW']), axis=0 ),
     np.mean( np.array(Data075['FW']), axis=0 ),
     np.mean( np.array(Data100['FW']), axis=0 ),
     np.mean( np.array(Data125['FW']), axis=0 ),
     np.mean( np.array(Data150['FW']), axis=0 ),
     np.mean( np.array(Data175['FW']), axis=0 ),
     np.mean( np.array(Data200['FW']), axis=0 )
 ]
color_data =[np.mean( np.array(Data250['LRn']), axis=0 ),
        np.mean( np.array(Data075['LRn']), axis=0 ),
      np.mean( np.array(Data100['LRn']), axis=0 ),
      np.mean( np.array(Data125['LRn']), axis=0 ),
      np.mean( np.array(Data150['LRn']), axis=0 ),
      np.mean( np.array(Data175['LRn']), axis=0 ),
      np.mean( np.array(Data200['LRn']), axis=0 )
]

size_data =[np.mean( np.array(Data250['TL']), axis=0 )*5,
        np.mean( np.array(Data075['TL']), axis=0 )*5,
      np.mean( np.array(Data100['TL']), axis=0 )*5,
      np.mean( np.array(Data125['TL']), axis=0 )*5,
      np.mean( np.array(Data150['TL']), axis=0 )*5,
      np.mean( np.array(Data175['TL']), axis=0 )*5,
      np.mean( np.array(Data200['TL']), axis=0 )*5
]
# WT_LRn_m = np.mean( np.array(Data_075['LRn_time']), axis=0 );WT_LRn_s =np.std( np.array(DataWT['LRn_time']), axis=0 )

fig, axis = plt.subplots(figsize=(4,3))
plt.scatter(x_data,y_data,s=size_data, c =color_data, cmap='Reds' ,edgecolors = 'grey', marker='s')
plt.scatter(x,y,s=size, c =color, cmap='Reds' ,edgecolors = 'grey')

plt.colorbar()
plt.savefig(path + '_scatter_data.png')
plt.savefig(path + '_scatter_data.svg')




plt.show()