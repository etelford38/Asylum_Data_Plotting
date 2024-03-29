# Asylum_Data_Plotting
Python codes to import, analyze, and save data generated by Asylum AFM instruments. Created by Evan Telford (ejt2133@columbia.edu) and Maelle Kapfer (mak2294@columbia.edu).

## All the included python files follow the general procedure:
1) Put all of the data to be analyzed in a folder called "Raw data".
2) Copy the corresponding python file into the folder containing "Raw data".
3) Run the python file. It will automatically create new folders containing plots of the raw and analyzed data.

NOTE: this section only applies to the original files. The updated files will create a "Raw data" folder on their own and organize the files accordingly.


NOTE: for the update folder (2023-11-02), all individual files have been consolidated into a single script. This single file can be placed inside a folder containing data files of various types (i.e. tapping, c-afm, contact, etc...) and, when run, will automatically sort the data based on the measurement type and will then analyze the data.

## To run all included python files, the following packages are required (links to anaconda):
* pathlib (https://anaconda.org/conda-forge/pathlib)
* os (https://anaconda.org/jmcmurray/os)
* igor (https://anaconda.org/conda-forge/igor)
* numpy (https://anaconda.org/anaconda/numpy)
* matplotlib (https://anaconda.org/conda-forge/matplotlib)
* opencv (https://anaconda.org/conda-forge/opencv)

## There are 7 python files in this project that perform the following tasks:
1) ["Asylum_CAFM.py"](https://github.com/etelford38/Asylum_Data_Plotting/blob/main/Asylum_CAFM.py)
* load, plot, and line correct conducting atomic force microscopy scans. 
* Both raw and analyzed plots are saved. 
* Both current preamplifier channels are plotted.
2) ["Asylum_CAFM_FFT.py"](https://github.com/etelford38/Asylum_Data_Plotting/blob/main/Asylum_CAFM_FFT.py)
* loads, plots, and fourier transforms conducting atomic force microscopy scans. 
3) ["Asylum_FD_trace.py"](https://github.com/etelford38/Asylum_Data_Plotting/blob/main/Asylum_FD_trace.py)
* loads and plots force-distance curves obtaining during approach scans with conducting tips. 
* Both current preamplifier channels are plotted.
4) ["Asylum_IV_trace.py"](https://github.com/etelford38/Asylum_Data_Plotting/blob/main/Asylum_IV_trace.py)
* loads and plots I-V curves. 
* Both current preamplifier channels are plotted.
5) ["Asylum_Tapping.py"](https://github.com/etelford38/Asylum_Data_Plotting/blob/main/Asylum_Tapping.py)
* loads, plots, and line corrects tapping scans. 
* Both raw and analyzed plots are saved.
* Topography, Amplitude Error, and Phase are plotted.
6) ["Asylum_Contact.py"](https://github.com/etelford38/Asylum_Data_Plotting/blob/main/original/Asylum_Contact.py)
* loads, plots, and line corrects contact-mode scans. 
* Both raw and analyzed plots are saved.
* Topography, Amplitude, and Z Sensor are plotted.
7) ["Asylum_LFM.py"](https://github.com/etelford38/Asylum_Data_Plotting/blob/main/original/Asylum_LFM.py)
* loads, plots, and line corrects lateral-force-mode scans. 
* Both raw and analyzed plots are saved.
* Topography, Amplitude, and Lateral Deflection are plotted.
