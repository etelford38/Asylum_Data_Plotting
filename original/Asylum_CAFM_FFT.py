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
import cv2
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

def oneFFT(data):
    f = np.fft.fft2(data)
    fshift = np.fft.fftshift(f)
    absFFT = np.abs(fshift) + 1e-10
    magnitude_spectrum = np.log(absFFT)
    magnitude_spectrum = magnitude_spectrum - np.mean(magnitude_spectrum)
    
    # Remove the zero-peak and the vertical/horizontal noise lines
    magnitude_spectrum_size = magnitude_spectrum.shape[0]
    magnitude_spectrum[int(magnitude_spectrum_size//2),:] = 0#np.std(magnitude_spectrum)
    magnitude_spectrum[:,int(magnitude_spectrum_size//2)] = 0#np.std(magnitude_spectrum)
    
    magnitude_spectrum_norm = magnitude_spectrum - np.min(magnitude_spectrum)
    magnitude_spectrum_norm = (magnitude_spectrum_norm / np.max(magnitude_spectrum_norm)*255).astype(np.uint8)
    
    magnitude_spectrum_norm = cv2.medianBlur(magnitude_spectrum_norm,1)
    
    test = magnitude_spectrum_norm.copy()
    
    return test

def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 0)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value
#%%
dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path=Path(dir_path)
name_path=dir_path/"Raw data"
names=name_path.glob('*.ibw')
for i,k in enumerate(names): 
    for l in range(1):
        data_file = k
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 6}
        mt.rc('font', **font)
        mt.rcParams['pdf.fonttype'] = 42
        mt.rcParams['ps.fonttype'] = 42
        mt.rcParams['font.family'] = 'Arial'
        fig, axes = plt.subplots(nrows=1,ncols=2,figsize=(7,5),gridspec_kw={'wspace':0.2,'hspace':0.0},dpi=300)
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
        voltage=metaData['SurfaceVoltage'].split('@')[0]
        voltagee=str(len(metaData['SurfaceVoltage'].split('@'))-1)
        deflection=metaData['DeflectionSetpointVolts'].split('@')[0]
        deflectione=str(len(metaData['DeflectionSetpointVolts'].split('@'))-1)
        scanrate=metaData['ScanRate'].split('@')[0]
        scanratee=str(len(metaData['ScanRate'].split('@'))-1)
        title_text=[label[0],label[6],label[8]]
        title_text_t=['$Current (nA),V=$'+voltage+'$V$','$FFT[Current] (a.u.)$']
        AFMZ=dat[title_text[0]]*1e9
        AFMD=dat[title_text[1]]*1e9
        
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
            if i==1:
                f.set_xlabel('$x^{-1} (nm^{-1})$', labelpad=0.1)
                f.set_ylabel('$y^{-1} (nm^{-1})$', labelpad=0.1)
            else:
                f.set_xlabel('$x (\mu m)$', labelpad=0.1)
                f.set_ylabel('$y (\mu m)$', labelpad=0.1)
                
        if l==0:
            AFMZ=line_correct(AFMZ,'vertical')
            AFMD=line_correct(AFMD,'vertical')
        
        AFMD=AFMD-np.mean(AFMD)
        
        AFMZmo, AMFZsd = optimize_color(AFMZ)
        AFMDmo, AMFDsd = optimize_color(AFMD)
########################################################################################################################################################################
        #FFT stuff
        x_lines = float(metaData['PointsLines'])
        y_lines = float(metaData['ScanLines'])
        x_real = float(scansize)
        y_real=x_real
        nmPrPix = (x_real/1e-9)/(x_lines)
        FFTextent = 1/(nmPrPix*2)
        FFTpixSize = FFTextent / (x_lines*0.5)
########################################################################################################################################################################
        #FFT padding if desired
        data=AFMD
        padding = int(x_lines/2)
        data_pad = np.pad(data, padding, pad_with)
        
        #Recalculate sizes for padding
        x_lines_pad = x_lines + padding*2
        y_lines_pad = y_lines + padding*2
        
        x_real_pad = x_real + padding*2*nmPrPix
        y_real_pad = y_real + padding*2*nmPrPix
        
        nmPrPix_pad = x_real_pad/x_lines_pad
        FFTextent_pad = FFTextent
        FFTpixSize_pad = FFTextent_pad / (x_lines_pad*0.5)
        zoomFactor_pad = 8
        
        ymin_pad,ymax_pad = int(y_lines_pad*(zoomFactor_pad-1)/(zoomFactor_pad*2)),int(y_lines_pad*(zoomFactor_pad+1)/(zoomFactor_pad*2))
        xmin_pad,xmax_pad = int(x_lines_pad*(zoomFactor_pad-1)/(zoomFactor_pad*2)),int(x_lines_pad*(zoomFactor_pad+1)/(zoomFactor_pad*2))
        
        subFFT_pad = oneFFT(data_pad)[ymin_pad:ymax_pad,xmin_pad:xmax_pad]
########################################################################################################################################################################
        f2=axes[0].imshow((data),extent = [0,float(scansize)*1e6*(x_lines)/x_lines,0,float(scansize)*1e6*(x_lines)/x_lines], cmap='Blues_r',vmin=AFMDmo-AMFDsd,vmax=AFMDmo+AMFDsd)
        c2=fig.colorbar(f2, ax=axes[0],fraction=0.046, pad=0.04,ticks=[])
        c2.outline.set_linewidth(0.5)
        
        f3=axes[1].imshow((subFFT_pad),
                          extent = [-FFTextent_pad/zoomFactor_pad,FFTextent_pad/zoomFactor_pad,-FFTextent_pad/zoomFactor_pad,FFTextent_pad/zoomFactor_pad],
                          cmap='Greens',
                          vmin=np.max(subFFT_pad)*0.75,vmax=0.9*np.max(subFFT_pad))
        c3=fig.colorbar(f3, ax=axes[1],fraction=0.046, pad=0.04,ticks=[])
        c3.outline.set_linewidth(0.5)       
        plt.tight_layout()
        
        save_path = Path(dir_path)/"FFT analysis publication quality"
        fname=k
        fname=os.path.split(fname)[1]
        if not os.path.exists(save_path):
           os.makedirs(save_path)
        os.chdir(save_path)
        plt.savefig(fname[:-4]+'_figure_padded.png', bbox_inches='tight',dpi=600)  
        plt.show()
