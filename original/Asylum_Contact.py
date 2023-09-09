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

def line_correct(array,direction):
    me=np.mean(array)
    temp_array=array
    y, x = array.shape
    if direction=='vertical':
        for i in range(x):
            temp_x=array[:,i]
            temp_x=temp_x-np.mean(temp_x)
            temp_array[:,i]=temp_x
    elif direction=='horizontal':
        for i in range(y):
            temp_y=array[i,:]
            temp_y=temp_y-np.mean(temp_y)
            temp_array[i,:]=temp_y
    else:
        print('not a valid direction')
    return temp_array+me

def optimize_color(array):
    mo=np.median(array)
    sd=np.std(array)
    return mo, 1.5*sd
#%%
dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path=Path(dir_path)
name_path=dir_path/"Raw data"
names=name_path.glob('*.ibw')
for i,k in enumerate(names): 
    for l in range(2):
        data_file = k
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 6}
        mt.rc('font', **font)
        mt.rcParams['pdf.fonttype'] = 42
        mt.rcParams['ps.fonttype'] = 42
        mt.rcParams['font.family'] = 'Arial'
        fig, axes = plt.subplots(nrows=1,ncols=3,figsize=(7,4),gridspec_kw={'wspace':0.5,'hspace':0.0},dpi=300)
        plt.rcParams['text.color'] = 'k'
        plt.rcParams['axes.labelcolor'] = 'k'
        plt.rcParams['xtick.color'] = 'k'
        plt.rcParams['ytick.color'] = 'k'
        dat = igor.load(data_file)
        metaData = {}
        metaData = getMetaData(dat)
        label, dat= getData2D(dat)
        scansize=metaData['ScanSize'].split('@')[0]
        scane=str(len(metaData['ScanSize'].split('@'))-1)
        deflection=metaData['DeflectionSetpointVolts'].split('@')[0]
        deflectione=str(len(metaData['DeflectionSetpointVolts'].split('@'))-1)
        scanrate=metaData['ScanRate'].split('@')[0]
        scanratee=str(len(metaData['ScanRate'].split('@'))-1)
        title_text=[label[0],label[2],label[4]]
        title_text_t=['Height (nm) - Errors: '+scane+deflectione+scanratee,'Amplitude (nV)','Z Sensor (nm)']

        AFMZ=dat[title_text[0]]*1e9
        AFMD=dat[title_text[1]]*1e9
        AFMI=dat[title_text[2]]*1e9
        
        for i,f in enumerate(axes):
            f.spines['left'].set_color('k')
            f.spines['right'].set_color('k')
            f.spines['bottom'].set_color('k')
            f.spines['top'].set_color('k')
            f.spines['left'].set_linewidth(0.5)
            f.spines['right'].set_linewidth(0.5)
            f.spines['bottom'].set_linewidth(0.5)
            f.spines['top'].set_linewidth(0.5)
            f.tick_params(width=0.5,length=2.0)
            f.set_title(title_text_t[i],fontsize=6, pad=5)
            f.tick_params(direction="in")
            f.yaxis.set_ticks_position(position='both')
            f.xaxis.set_ticks_position(position='both')
            f.set_xlabel('$x (\mu m)$', labelpad=0.1)
            f.set_ylabel('$y (\mu m)$', labelpad=0.1)
        
        if l==0:
            pass
        if l==1:
            AFMZ=line_correct(AFMZ,'vertical')
            AFMD=line_correct(AFMD,'vertical')
            AFMI=line_correct(AFMI,'vertical')

        AFMZmo, AMFZsd = optimize_color(AFMZ)
        AFMDmo, AMFDsd = optimize_color(AFMD)
        AFMImo, AMFIsd = optimize_color(AFMI)
            
        f1=axes[0].imshow((AFMZ),extent = [0,float(scansize)*1e6,0,float(scansize)*1e6],cmap='Greys_r',vmin=AFMZmo-AMFZsd,vmax=AFMZmo+AMFZsd)
        c1=fig.colorbar(f1, ax=axes[0],fraction=0.046, pad=0.04)
        c1.outline.set_linewidth(0.5)
        f2=axes[1].imshow((AFMD),extent = [0,float(scansize)*1e6,0,float(scansize)*1e6], cmap='Blues_r',vmin=AFMDmo-AMFDsd,vmax=AFMDmo+AMFDsd)
        c2=fig.colorbar(f2, ax=axes[1],fraction=0.046, pad=0.04)
        c2.outline.set_linewidth(0.5)
        f3=axes[2].imshow((AFMI),extent = [0,float(scansize)*1e6,0,float(scansize)*1e6], cmap='hot',vmin=AFMImo-AMFIsd,vmax=AFMImo+AMFIsd)
        c3=fig.colorbar(f3, ax=axes[2],fraction=0.046, pad=0.04)
        c3.outline.set_linewidth(0.5)
        plt.tight_layout()
        
        if l==0:
            save_path = Path(dir_path)/"Images uncorrected"
            fname=k
            fname=os.path.split(fname)[1]
            if not os.path.exists(save_path):
               os.makedirs(save_path)
            os.chdir(save_path)
            plt.savefig(fname[:-4]+'_figure_uncorrected.png', bbox_inches='tight',dpi=600)   
        if l==1:
            save_path = Path(dir_path)/"Images line corrected"
            fname=k
            fname=os.path.split(fname)[1]
            if not os.path.exists(save_path):
               os.makedirs(save_path)
            os.chdir(save_path)
            plt.savefig(fname[:-4]+'_figure_line_corrected.png', bbox_inches='tight',dpi=600)  
        plt.show()
