[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bethel-physics/ParticleDiscoveryLab/HEAD) [![DOI](https://zenodo.org/badge/132026337.svg)](https://zenodo.org/doi/10.5281/zenodo.13787598)

# Particle Discovery lab (for instructors!)

The particle discovery lab uses [CMS dimuon data from 2016](http://doi.org/10.7483/OPENDATA.CMS.UZD7.Z50M) published via the [CERN Open Data Portal](http://opendata.cern.ch/). 
We have developed an undergraduate intermediate-level lab exercise to complement the many high school-level exercises available via the Open Data Portal.
Solutions and student code are available in both MATLAB and Python, and do not require special CMS Open Data software for students or instructors. Instructions for using either software option are in separate sections below.

The goal of this exercise is for students to reconstruct decays of unknown particle X (initial state) to 2 muons (final state). They will use histograms to display their calculated mass for particle X, and learn about fitting and subtracting background contributions from data. Uncertainty propagation concepts are included through each step of the analysis. After isolating the signal distribution they will determine which particle they have discovered and compare their observed properties (mass and width) to the known properties. 

![J/psi](images/MuonLab_JpsiSigBkg.png)

## Visualize the data
To help introduce this exerience, there are event displays available from the "DoubleMuParked" 2012 dataset that is conceptially identical to the dataset we are studying. The [ISpy Event Display](http://opendata.cern.ch/visualise/events/cms#) is a web-based tool to study events interactively, so students could do this quickly on their laptops. To get
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
**Software setup**: [MATLAB](https://www.mathworks.com/products/matlab.html) should be installed on whatever computers students will use for the exercise. The "statistics and machine learning toolbox" is the only package used beyond the basic MATLAB install. This package might be part of your MATLAB license. Write to Julie Hogan for options to use the "curve-fitting toolbox" instead. This exercise is written using MATLAB "live scripts", which can display text and instructions with typical web formatting and provide smaller code boxes for students to use.

**Exercise materials**:

 * Download this Github code!
 * Post the student blank (`MuonAnalysis_student.mlx`) and the workspace (`DoubleMuon_2016H_200k.mat`) on your LMS.
 * Download the `pollsf.m` function file from Alejandro Garcia's [*Numerical Methods for Physics*](https://github.com/AlejGarcia/NM4P/tree/master/MatlabRevised) repository and post it for students on the LMS.

### Understanding the workspace
The workspace stored in `DoubleMuon_2016H_200k.mat` contains the energy and 3-momentum for "muon 1" and "muon 2". When you load the workspace you'll see 4 variables: 
arrays of size (2 x nEvents) for `E`, `px`, `py`, and `pz`. The default unit for these numerical values is the GeV.

Values for "muon 1" are stored in the first index, and values for "muon 2" are stored in the second index. Energies can be accessed like this:
```
i = 294;   # Let's look at event number 294
E_muon1 = E(i,1);  # MATLAB counts indices from 1, not 0!
E_muon2 = E(i,2);
```

### To start from the 200k event workspace
Students can start by opening `MuonAnalysis_student.mlx` with MATLAB. This contains the instruction of how to set their file path and load the workspace. 
Your reference is `MuonAnalysis_solutions.mlx`. This is a fully-coded solutions guide that you can run in MATLAB to get a set of solutions plots. 
You will want to **provide students a mass window** in which to search. Tested, working suggestions include:
 * J/psi meson: 2.8 GeV -- 3.5 GeV
 * Upsilon mesons: 8 -- 12 GeV (there are 3 peaks here, so I give this window to the strongest group)
 * Z boson: 60 -- 120 GeV (or for a more challenging fit, assign 30 -- 150)
 * [More ideas](https://github.com/cms-opendata-analyses/DimuonSpectrumNanoAODOutreachAnalysis)
 

### To reprocess the data and make a new input workspace
 * `NanoAOD2Arrays.py` can be run to produce a text file with as many events as you want! See the comments in that file for accessing the NanoAOD file either via the web or by downloading a local file.
   * You will need the ROOT program installed to make a new workspace -- check out the [CMS Open Data ROOT docker container](https://opendata.cern.ch/docs/cms-guide-docker) for a simple work environment with all the required software.
   * In the docker container, execute `python NanoAOD2Arrays.py` to create a new workspace.
   * You need the output file called `DoubleMuon_2016H_200k.txt` (or whatever new name you choose!) for the next step.
 * The MATLAB script `txtFileReader.m` can be run in MATLAB to create a workspace from the text file. The variables will appear in the Workspace panel to the right -- you can see the size/shape of each array and double-click on a variable to see its contents. 
   * In MATLAB, click on the Workspace panel and type `CTRL-s`. Save a `.mat` file and post it for your students to download. 

## Python (PythonAnalysis folder)

**Software setup**: Python 3 should be installed on whatever computers students will use for the exercise. There are many methods for students to access Python:
 * If your students' computational courses use Python via a certain program, such as Anaconda, that will almost certainly be the simplest thing to use.
 * Binder can be used to run the Juypter notebooks on the web! See below for instructions
 * The [CMS Open Data Python docker container](https://opendata.cern.ch/docs/cms-guide-docker) can be installed on lab or individual computers and used for the exercises.
 * Contact Julie Hogan to brainstorm python solutions if needed. 

The following python packages are needed for the exercise:
 * [matplotlib](https://matplotlib.org/)
 * [PyQt6](https://pypi.org/project/PyQt6/)
 * [pickle](https://docs.python.org/3/library/pickle.html)
 * [numpy](https://numpy.org/)
 * [scipy](https://www.scipy.org/)
 * [jupyter](https://jupyter.org/) -- optional, if using the notebooks.
 * Use `pip install <packagename>` to install any of these packages

**Binder**: You can run these Jupyter notebooks on the web! https://mybinder.org/v2/gh/bethel-physics/ParticleDiscoveryLab/HEAD

If you fork this github repository to make it your own, you can provide students your own Binder link. Go to [mybinder.org](mybinder.org) and enter the link to your own github repository, such as: https://github.com/<yourUserNameHere>/ParticleDiscoveryLab.

**Exercise Materials**:

 * Download this Github code!
 * Post the student blank (`MuonAnalysis_student.py`), the pickle file (`DoubleMuon_2016H_200k.pkl`), and `pollsf.py` on your LMS. The `pollsf.py` function file is adapted directly from Alejandro Garcia's [*Numerical Methods for Physics*](https://github.com/AlejGarcia/NM4P/tree/master/MatlabRevised) repository.
 * If your students will be able to use Jupyter notebooks, you can provide the `MuonAnalysis_student.ipynb` file instead. 

### Understanding the workspace
The data is stored in `DoubleMuon_2016H_200k.pkl` as a "list of lists". Each item in the list contains 8 numerical values in this order: `E1`, `E2`, `px1`, `px2`, `py1`, `py2`, `pz1`, `pz2` (where "1" refers to "muon 1" and "2" refers to "muon 2") The default unit for these numerical values is the GeV.

Energies for each muon can be accessed like this:
```
data = pickle.load(open('DoubleMuon_2016H_200k.pkl','rb'))
i = 294;   # Let's look at event number 294
E_muon1 = data[i][0];  # Python counts indices from 0, different from MATLAB
E_muon2 = data[i][1];
```

### To start from the 200k event workspace
Students can start directly from `MuonAnalysis_student.py` (or the corresponding notebook), which is already set to import packages and load the pickle file:
```
$ cd /path/to/files/they/downloaded/
$ vi MuonAnalysis_student.py  # or their favorite text editor
```

To use jupyter, students should launch a notebook, wait for the web browser to launch (or paste the link), and then click on the MuonAnalysis_student.ipynb file. The instructions to launch a jupyter notebook vary based on the python platform used -- an example is provided here, but adaptation might be necessary. Please test this with your systems during exercise prep.
```
$ cd /path/to/files/they/downloaded/
$ jupyter notebook  # in the Open Data docker container use: jupyter-lab --ip=0.0.0.0 --no-browser
[I 10:53:47.697 NotebookApp] Serving notebooks from local directory: /home/...somepath.../MatlabOpenData
[I 10:53:47.697 NotebookApp] 0 active kernels 
[I 10:53:47.697 NotebookApp] The Jupyter Notebook is running at: http://localhost:8888/?token=lonnnnnggggstringofrandomletters
[I 10:53:47.697 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 10:53:47.697 NotebookApp] 
    
    Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://localhost:8888/?token=lonnnggggggstringofrandomletters
```

Your reference is `MuonAnalysis_solutions.py` (again with a corresponding notebook), which can be opened in similar ways. The jupyter notebook can be viewed directly on github! This is a fully-coded solutions guide that you can run (`python -u MuonAnalysis_solutions.py`) to get a set of solutions plots. 
You will want to **provide students a mass window** in which to search. Tested, working suggestions include:
 * J/psi meson: 2.8 GeV -- 3.5 GeV
 * Upsilon mesons: 8 -- 12 GeV (there are 3 peaks here, so I give this window to the strongest group)
 * Z boson: 60 -- 120 GeV (or for a more challenging fit, assign 30 -- 150)
 * [More ideas](https://github.com/cms-opendata-analyses/DimuonSpectrumNanoAODOutreachAnalysis)
 

### To reprocess the data and make a new input workspace
 * `NanoAOD2Arrays.py` can be run to produce a text file with as many events as you want! See the comments in that file for accessing the NanoAOD file either via the web or by downloading a local file.
   * You will need the ROOT program installed to make a new workspace -- check out the [CMS Open Data ROOT docker container](https://opendata.cern.ch/docs/cms-guide-docker) for a simple work environment with all the required software.
   * In the docker container, execute `python NanoAOD2Arrays.py` to create a new workspace.
   * You need the output file called `DoubleMuon_2016H_200k.pkl` (or whatever new name you choose!) to provide to your students.

