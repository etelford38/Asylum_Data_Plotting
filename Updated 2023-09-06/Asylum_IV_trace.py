# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 15:51:16 2023

@author: Evan Telford (ejt2133@columbia.edu) and Maelle Kapfer (mak2294@columbia.edu)
"""
#%%
#import pertinent packages
import igor.binarywave as igor
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib as mt
from pathlib import Path
import shutil
#%%
def getMetaData(dat):
    note = dat['wave']['note'].decode('ISO-8859-1')
    metaData = {}
    name = []
    val = []
    for line in note.split('\r'):
        name.append(line.split(':')[0])
        try:
            val.append(line.split(':')[1])
        except IndexError:
            val.append('nan')
    metaData = dict(zip(name,val))
    
    return metaData

def getData_Line(dat):
    data = dat['wave']['wData']
    Label = dat['wave']['labels'][1]
    Label = Label[1:]
    dat = np.zeros((len(Label),len(data[:,0])))
    dict_res = {}
    
    for i in range(len(Label)):
        dat[i] = data[:,i]
        
    for l,d in zip(Label,dat):
        dict_res[l] = d

    return Label, dict_res

def getData2D(dat):
    data = dat['wave']['wData']
    Label = dat['wave']['labels'][2]
    Label = Label[1:]
    dat = np.zeros((len(Label),len(data),len(data)))
    dict_res = {}
    
    for i in range(len(Label)):
        dat[i] = data[:,:,i]
        
    for l,d in zip(Label,dat):
        dict_res[l] = d
    
    return Label, dict_res
#%%
key='*_IV_*'
xl='$V (V)$'
x=b'Bias'
#sets up the basic directory
dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path=Path(dir_path)
names=dir_path.glob('*.ibw')

#creates Raw data folder AND moves the files
name_path=Path(dir_path)/"Raw data"
if not os.path.exists(name_path):
   os.makedirs(name_path)

for m in names:
    base_name=os.path.basename(m)
    shutil.move(m,name_path/base_name)

#change directory to Raw data and load the files
os.chdir(name_path)
names=name_path.glob('*.ibw')

for i,j in enumerate(names):  
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 8}
    mt.rc('font', **font)
    mt.rcParams['pdf.fonttype'] = 42
    mt.rcParams['ps.fonttype'] = 42
    mt.rcParams['font.family'] = 'Arial'
    fig,axes = plt.subplots(nrows=1,ncols=2,figsize=(7,3),gridspec_kw={'wspace':0.4,'hspace':0.0},dpi=600)
    plt.rcParams['text.color'] = 'k'
    plt.rcParams['axes.labelcolor'] = 'k'
    plt.rcParams['xtick.color'] = 'k'
    plt.rcParams['ytick.color'] = 'k'
    
    for l in range(2):
        f=axes[l]
        g=f.twinx()
        g.set_zorder(1)
        f.set_zorder(0)
        f.spines['left'].set_color('k')
        f.spines['right'].set_color('k')
        f.spines['bottom'].set_color('k')
        f.spines['top'].set_color('k')
        f.spines['left'].set_linewidth(0.5)
        f.spines['right'].set_linewidth(0.5)
        f.spines['bottom'].set_linewidth(0.5)
        f.spines['top'].set_linewidth(0.5)
        f.tick_params(width=0.5)
        dat = igor.load(j)
        metaData = {}
        metaData = getMetaData(dat)
        label, dact= getData_Line(dat)
        
        x_temp=(dact[x])
        if l==0:
            I_temp=(dact[b'Cur']-dact[b'Cur'].mean())*1e9
        else:
            I_temp=(dact[b'Curr2']-dact[b'Curr2'].mean())*1e9
        F_temp=(dact[b'Defl']-np.mean(dact[b'Defl']))*1e9
        
        f.plot(x_temp, F_temp ,'ko-',markersize=1,linewidth=1)
        f.set_xlabel(xl,labelpad=0.0)
        f.set_ylabel('$\Delta F (nN)$',labelpad=0.0)
        f.set_ylim(-1,1)
        g.plot(x_temp,I_temp,'ro-',markersize=1,linewidth=1)
        g.set_xlabel(xl,labelpad=0.1)
        g.set_ylabel('$I (nA)$',labelpad=0.0)
        
        voltage=metaData['DeflectionSetpointVolts'].split('@')[0]
        voltagee=str(len(metaData['DeflectionSetpointVolts'].split('@'))-1)
        f.set_title('$F=$'+voltage+'$V$ - Errors: '+voltagee,fontsize=8,pad=3)

        f.patch.set_visible(False)
        f.tick_params(direction="in")
        g.tick_params(direction="in")
        f.xaxis.set_ticks_position(position='both')
        f.set_xlim(np.min(x_temp),np.max(x_temp))
    save_path = Path(dir_path)/"Plots uncorrected"
    fname=j
    fname=os.path.split(fname)[1]
    if not os.path.exists(save_path):
       os.makedirs(save_path)
    os.chdir(save_path)
    plt.savefig(fname[:-4]+'_figure_uncorrected.png', bbox_inches='tight',dpi=600)   
    plt.show()