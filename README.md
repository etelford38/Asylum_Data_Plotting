# Asylum_Data_Plotting
Python codes to import, analyze, and save data generated by Asylum AFM instruments. Created by Evan Telford (ejt2133@columbia.edu) and Maelle Kapfer (mak2294@columbia.edu).

All the included python files follow the general procedure:
1) Put all of the data to be analyzed in a folder called "Raw data".
2) Copy the corresponding python file into the folder containing "Raw data".
3) Run the python file. It will automatically create new folders containing plots of the raw and analyzed data.

To run all included python files, the following packages are required (links are to anaconda installation information):
* pathlib (https://anaconda.org/conda-forge/pathlib)
* os (https://anaconda.org/jmcmurray/os)
* igor (https://anaconda.org/conda-forge/igor)
* numpy (https://anaconda.org/anaconda/numpy)
* matplotlib (https://anaconda.org/conda-forge/matplotlib)
* opencv (https://anaconda.org/conda-forge/opencv)

There are 5 python files in this project that perform the following tasks:
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
