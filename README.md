# MatlabOpenData
**************************************** MATLAB INSTRUCTIONS ****************************************

1. Download MuonAnalysisMaster.m
2. Download muons100k.mat
    a. Or use txtFileReader.m to make your own workspace/.mat file
        i. Use NanoAOD2Arrays.py (in PythonAnalysis folder) to create custom txt file
        ii. Make sure to save the txt file
        iii. Run txtFileReader.m for your txt file to create a custom .mat file
3. Download the function pollsf.m from https://github.com/AlejGarcia/NM4P/tree/master/Matlab
4. Make sure the .mat file and pollsf.m are in the Matlab path
5. Load workspace with command "load ('workspaceName.mat')"
6. Run

**************************************** PYTHON INSTRUCTIONS ****************************************

1. Download MuonAnalysisMaster_python2020.py, DoubleMuParked_100k.py, and pollsf.py from the github
2. Install these packages for python: matplotlib, pickle, numpy, and scipy
    Links: matplotlib: https://matplotlib.org/
           pickle: https://docs.python.org/3/library/pickle.html
           numpy: https://numpy.org/
           scipy: https://www.scipy.org/
3. Make sure MuonAnalysisMaster_python2020.py, DoubleMuParked_100k.py, and pollsf.py are in the same file
4. Run MuonAnalysisMaster_python2020.py

OPTIONAL: If you would like to generate your own .pkl file of a different size, download and run NanoAOD2Arrays.py.
