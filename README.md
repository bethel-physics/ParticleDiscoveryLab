[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bethel-physics/ParticleDiscoveryLab/HEAD)

# Particle Discovery lab (for students!)

The particle discovery lab uses [CMS dimuon data from 2012](http://doi.org/10.7483/OPENDATA.CMS.LVG5.QT81) published via the [CERN Open Data Portal](http://opendata.cern.ch/). 
We have developed an undergraduate intermediate-level lab exercise to complement the many high school-level exercises available via the Open Data Portal.
Solutions and student code are available in both MATLAB and Python, but do not require ROOT or Open Data Virtual Machines for students or instructors.

The goal of this exercise is for students to reconstruct decays of unknown particle X (initial state) to 2 muons (final state). They will use histograms to display their calculated mass for particle X, and learn about fitting and subtracting background contributions from data. Uncertainty propagation concepts are included through each step of the analysis. After isolating the signal distribution they will determine which particle they have discovered and compare their observed properties (mass and width) to the known properties. 

![J/psi](images/MuonLab_JpsiSigBkg.png)

## Visualize the data
To help introduce this exerience, there are event displays available from the DoubleMuParked dataset that we are studying. The [ISpy Event Display](http://opendata.cern.ch/visualise/events/cms#) is a web-based tool to study events interactively, so students could do this quickly on their laptops. To get
to our dataset:
 * Click on the "folder" icon in the upper left corner
 * Select "Open files from the web"
 * Single-click on "Run2012B/"
 * Scroll down and single-click on "DoubleMuParked_0.ig"
 * Click on any of the individual events that appear in the right-hand column
 * Click "Load"
Click and drag to rotate the image around! The yellow cylinder represents a detector element for scale, and the red boxes show which muon detector elements were hit by the 2 muons in each event. More detector or physics elements can be showed by making selections in the left-hand menu

![ISpy event display](images/eventDisplay.PNG)

## MATLAB (MatlabAnalysis folder)
**Software setup**: MATLAB should be installed on whatever computers students will use for the exercise. The "curve-fitting toolbox" is the only package
used beyond the basic MATLAB install. This package might be part of your MATLAB license, or trials are available for 30 days. https://www.mathworks.com/products/matlab.html

### Understanding the workspace
The workspace stored in `DoubleMuParked_100k.mat` contains the energy and 3-momentum for "muon 1" and "muon 2". When you load the workspace you'll see 4 variables: 
arrays of size (2 x nEvents) for `E`, `px`, `py`, and `pz`. The default unit for these numerical values is the GeV.

Values for "muon 1" are stored in the first index, and values for "muon 2" are stored in the second index. Energies can be accessed like this:
```
i = 294;   # Let's look at event number 294
E_muon1 = E(i,1);  # MATLAB counts indices from 1, not 0!
E_muon2 = E(i,2);
```

### To start from the 100k event workspace
Students can start by opening `MuonAnalysis_student.mlx` with MATLAB. This contains the instruction of how to set their file path and load the workspace. 


### To reprocess the data and make a new input workspace
 * `NanoAOD2Arrays.py` can be run (`python -u NanoAOD2Arrays.py`) to produce a text file with as many events as you want -- up to the ~61M stored in the NanoAOD ROOT file. See the comments in that file for accessing the NanoAOD file either via the web or by downloading a local file. You will need the ROOT program installed, either [directly](https://root.cern.ch/) or by setting up the [Open Data Software](http://opendata.cern.ch/docs/cms-virtual-machine-2011).
 * `txtFileReader.m` can be run in MATLAB to create a workspace from the text file. The variables will appear in the Workspace panel to the right -- you can see the size/shape of each array and double-click on a variable to see its contents. 
 * In MATLAB, click on the Workspace panel and type CTRL-s. Save a .mat file and post it for your students to download. 

## Python (PythonAnalysis folder)

**Binder**: You can run these Jupyter notebooks on the web! https://mybinder.org/v2/gh/bethel-physics/ParticleDiscoveryLab/HEAD

**Non-binder software setup**: If not using Binder, Python should be installed on whatever computers students will use for the exercise. My experience is using [Cygwin](https://www.cygwin.com/) on Windows, which provides a unix-based terminal in which python can be run. There are many other methods -- if your students' computational courses use Python via a certain program, I recommend following that protocol. Contact me (j-hogan@bethel.edu) to brainstorm python solutions if needed. 

Install these packages for python:
 * [matplotlib](https://matplotlib.org/)
 * [pickle](https://docs.python.org/3/library/pickle.html)
 * [numpy](https://numpy.org/)
 * [scipy](https://www.scipy.org/)
 * [jupyter](https://jupyter.org/) optional, for using the notebooks. 


If you fork this github repository to make it your own, you can create your own Binder link. Go to [mybinder.org](mybinder.org) and enter the link to your own github repository, such as: https://github.com/<yourUserNameHere>/ParticleDiscoveryLab.

### Understanding the workspace
The data is stored in `DoubleMuParked_100k.pkl` as a "list of lists". Each item in the list contains 8 numerical values in this order: `E1`, `E2`, `px1`, `px2`, `py1`, `py2`, `pz1`, `pz2` (where "1" refers to "muon 1" and "2" refers to "muon 2") The default unit for these numerical values is the GeV.

Energies for each muon can be accessed like this:
```
data = pickle.load(open('DoubleMuParked_100k.pkl','rb'))
i = 294;   # Let's look at event number 294
E_muon1 = data[i][0];  # Python counts indices from 0, different from MATLAB
E_muon2 = data[i][1];
```

### To start from the 100k event workspace
Students can start directly from `MuonAnalysis_student.py`, which is already set to import packages and load the pickle file. The simplest method is to click on the Binder link and navigate to `MuonAnalysis_student.ipynb`. Or, using python locally: 
```
$ cd /path/to/files/they/downloaded/
$ vi MuonAnalysis_student.py  # or their favorite text editor
```
To use jupyter outside of Binder, students should launch a notebook, wait for the web browser to launch (or paste the link), and then click on the MuonAnalysis_student.ipynb file.
```
$ cd /path/to/files/they/downloaded/
$ jupyter notebook
[I 10:53:47.697 NotebookApp] Serving notebooks from local directory: /home/...somepath.../MatlabOpenData
[I 10:53:47.697 NotebookApp] 0 active kernels 
[I 10:53:47.697 NotebookApp] The Jupyter Notebook is running at: http://localhost:8888/?token=lonnnnnggggstringofrandomletters
[I 10:53:47.697 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 10:53:47.697 NotebookApp] 
    
    Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://localhost:8888/?token=lonnnggggggstringofrandomletters
```

### To reprocess the data and make a new input workspace
 * `NanoAOD2Arrays.py` can be run to produce a new pickle file with as many events as you want -- up to the ~61M stored in the NanoAOD ROOT file. See the comments in that file for accessing the NanoAOD file either via the web or by downloading a local file. You will need the ROOT program installed, either [directly](https://root.cern.ch/) or by setting up the [Open Data Software](http://opendata.cern.ch/docs/cms-virtual-machine-2011).
 * Execute `python -u NanoAOD2Arrays.py` and provide students with the new pickle file.  
