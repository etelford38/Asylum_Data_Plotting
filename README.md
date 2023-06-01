# Asylum_Data_Plotting
Python codes to import, analyze, and save data generated by Asylum AFM instruments. Created by Evan Telford and Maelle Kapfer.

All the included python files follow the format:
1) Put all of the data to be analyzed in a folder called "Raw data".
2) Copy the corresponding python file into the folder containing "Raw data".
3) Run the python file. It will automatically create new folders containing plots of the raw and analyzed data.

To run all included python files, the following packages are required (links are to anaconda installation information):
1) pathlib (https://anaconda.org/conda-forge/pathlib)
2) os (https://anaconda.org/jmcmurray/os)
3) igor (https://anaconda.org/conda-forge/igor)
4) numpy (https://anaconda.org/anaconda/numpy)
5) matplotlib (https://anaconda.org/conda-forge/matplotlib)
6) opencv (https://anaconda.org/conda-forge/opencv)

There are 5 python files in this project that perform the following tasks:
1) "Asylum_CAFM.py"
        load, plot, and line correct conducting atomic force microscopy scans. Both raw and analyzed plots are saved. Both current preamplifier channels are plotted.
3) "Asylum_CAFM_FFT.py": load, plot, and fourier transforms conducting atomic force microscopy scans. 
4) "Asylum_FD_trace.py": loads and plots force-distance curves obtaining during approach scans with conducting tips. Both current preamplifier channels are plotted.
5) "Asylum_IV_trace.py": loads and plots I-V curves. Both current preamplifier channels are plotted.
6) "Asylum_Tapping.py": loads, plots, and line corrects tapping scans. Both raw and analyzed plots are saved.
