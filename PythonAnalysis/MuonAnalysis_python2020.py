######################################################
#
#  June 2020
#
#  Script to read the Open Data parked dimuon sample
#  and produce a "data" object with the 4-vector
#  that is saved in a pickle file. 
#
#  Julie Hogan, j-hogan@bethel.edu
#
######################################################

import math
import matplotlib.pyplot as plt
import pickle
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import pollsf

#attempt with pickle
#data = pickle.load(open('DoubleMuParked_100k.pkl','rb'))
data = pickle.load(open('DoubleMuParked_1M.pkl','rb'))

# This program begins by loading a Txt file that contains the  
# momentum and energy data from about 1M events. It produces a 
# graph that has mass on the x-axis and the number of muon pairs 
# with that mass on the y-axis. 
# The user can choose the range of GeV values shown. 

## DAY 1 = RECONSTRUCTION

# Choose how many events to process
#Ntoprocess = input("How many events to process? ")
Ntoprocess=1000000
Ntoprocess = int(Ntoprocess)

# Initialize a vector that will hold 1 invariant mass per event
Masses = []
KineticEnergy = []

# User-chosen min/max values and resolution
#Min = input("Type in your min (in GeV): ")
Min=2.8
Min = float(Min)
#Max = input("Type in your max (in GeV): ")
Max=3.5
Max = float(Max)
#n = input("Type in the number of bins: ")
n=50
n = int(n)
BinWidth = (Max - Min)/n

E=0
xMom=0
M=0
Px=0
Py=0
Pz=0
xK=0

# Loop over the number of events with at least 2 muons
print("Looping over ", Ntoprocess, " events...", sep=" ")
for i in range(Ntoprocess):
    # COMPUTE the mass of particle X -> mu mu
    
    E = data[i][0] + data[i][1]
    Px = data[i][2] + data[i][3]
    Py = data[i][4] + data[i][5]
    Pz = data[i][6] + data[i][7]
    M = math.sqrt(E**2 - Px**2 - Py**2 - Pz**2)
    Masses.append(M.real)
    


    # THINK: Is this mass in your window from Min to Max? 
    #        What should you do if it"s outside the window?


    
    # Calculate the Kinetic Energy of particle X. 
    # Store KE and mass values to plot later
    # Tip: make sure you mass value if "real" by using real(massvalue)
    if M> Min and M<Max:
        M=(M).real
        xK=E-M
        KineticEnergy.append(xK)


#end of the for loop

# HISTOGRAMMING -- create mass and KE histograms
# THINK: What do you expect your kinetic energy histogram to look like?
#        What do you expect your mass histogram to look like?
#        Make a quick sketch of what you expect for both plots
#
# Vocab: imagine plot with 3 bins on x-axis: 0-10, 10-20, 20-30
# "Bin edges": 0, 10, 20, 30
# "Bin centers": 5, 15, 25 (want dots on plot to be here!)
# "Bin width": 10



# Draw a HISTOGRAM of counts versus mass
# HISTOGRAM needs: Mass values list, Mass "bin edges" list 
plt.figure()
x, bins, patches=plt.hist(x=Masses, bins=n, range=(Min, Max), histtype='step',alpha=.5)
center = (Max - Min)/n
centers = []
for i in range(len(bins)-1):
    centers.append(bins[i]+center*.5)

err=[[],[]]
for i in range(len(x)):
    if x[i]==0:
        err[1].append(1)
    else:
        err[1].append( math.sqrt(x[i]))
    err[0].append( math.sqrt(x[i]))

plt.errorbar(centers, x, yerr=err, fmt='none', ecolor='red')
plt.xlabel('GeV') 
plt.ylabel('Number of muon pairs') 
  
plt.title('MuonLab GeV/Counts\n\n', fontweight ="bold") 
#plt.show()

counts=x
binEdg=bins
error=err
midpoints=centers



# THINK: What should the ERROR BARS be for each bin? 
#        What should you do if the bin has ZERO entries?
# Tools: 
# HISTCOUNTS: gives counts, needs bin edges + input values (mass)
# ERRORBAR: draws dots+bars, needs bin centers, y values, down 
#           uncert list, up uncert list 




# Draw another HISTOGRAM of counts vs kinetic energy
# Add ERROR BARS
plt.figure()
x, bins, patches=plt.hist(x=KineticEnergy, bins=n, range=(0, 800),histtype='step', log=True,alpha=.5)
center = (800)/n
centers=[]
for i in range(len(bins)-1):
    centers.append(bins[i]+center*.5)

err=[[],[]]
for i in range(len(x)):
    if x[i]==0:
        err[1].append(1)
    else:
        err[1].append( math.sqrt(x[i]))
    err[0].append( math.sqrt(x[i]))

plt.errorbar(centers, x, yerr=err, fmt='none', ecolor='red')
plt.xlabel('KE') 
plt.ylabel('Number of muon pairs') 
  
plt.title('MuonLab KE/Counts\n\n', fontweight ="bold") 
#plt.show()



### Great work! SAVE these plots to represent your RAW DATA.
###             SAVE a DAY1 workspace (delete large E, px, py, pz)

## DAY 2 = FITTING -- fit background on either side of the y
# LOAD your DAY1 workspace

# Vocab: imagine a mass plot with a bump in the middle
# "y window": region along x-axis under the y
# "background": smoothly falling slope of random events, 
#               including some of the events in the y window
# "signal": events in the y window minus the background

# CHOOSE mass values in GeV for where the y lies. 



# REMOVE the y window completely from your list of: 
# bin edges, bin centers, counts, uncertainties. 
# This forms your BACKGROUND dataset


bkgCount=[]
bkgBin=[]
bkgErr=[[],[]]
bkgMids=[]
for i in range(17):
    bkgCount.append(counts[i])
    bkgBin.append(binEdg[i])
    bkgErr[0].append(error[0][i])
    bkgErr[1].append(error[1][i])
    bkgMids.append(midpoints[i])
for i in range(17):
    bkgCount.append(counts[i+33])
    bkgBin.append(binEdg[i+33])
    bkgErr[0].append(error[0][i+33])
    bkgErr[1].append(error[1][i+33])
    bkgMids.append(midpoints[i+33])
    
bkgBin.append(binEdg[50])
#print(len(bkgBin))
#print(bkgMids)
#print(bkgCount)
#print("{")
#for i in range(len(bkgMids)-1):
#    print("{", bkgMids[i], ",", bkgCount[i],"},", sep=" ")
#print("{", bkgMids[44], ",", bkgCount[44],"} }", sep=" ")


# PERFORM a polynominal fit to the background
# THINK: Which type of curve do you expect will match your data best?
#        Imagine a curve connecting the two sides under your y.
# Tool: POLLSF gives fit params, uncerts, y-values, chi^2 value
#       needs bin centers, counts, uncertainties, N params
#       This is a least-squares fitter that uses uncerts!
m=2
fitparams, uncerts, yvalues, chisq = pollsf.pollsf(bkgMids, bkgCount, bkgErr[1], m)
print(chisq/(50-m))


# EVALUATE your fit by chi^2 and plotting
# -- Plotting: does the shape make any sense? Make a helpful plot
# -- Chi^2 (or "SSE") is defined in Eq. 29. It describes the difference 
#    between the points and the fitted curve. LARGER chi^2 tends to mean 
#    more difference or scatter of points.
#    OPTIMALLY, Chi^2 / (# points - # parameters) is around 1
# REPEAT fitting until you are satisfied
plt.figure()
x, bins, patches=plt.hist(x=Masses, bins=n, range=(Min, Max), histtype='step',alpha=.5)
center = (Max - Min)/n
centers = []
for i in range(len(bins)-1):
    centers.append(bins[i]+center*.5)

err=[[],[]]
for i in range(len(x)):
    if x[i]==0:
        err[1].append(1)
    else:
        err[1].append( math.sqrt(x[i]))
    err[0].append( math.sqrt(x[i]))

plt.errorbar(centers, x, yerr=err, fmt='none', ecolor='red')
plt.xlabel('GeV') 
plt.ylabel('Number of muon pairs') 
  
plt.title('Fit Function\n\n', fontweight ="bold") 
plt.plot(bkgMids,yvalues)
#plt.show()



## SUBTRACTION -- now you will subtract that background from data
# THINK: How will you estimate background in the signal y window?
#        What do you expect the curve to look like after bkg subtraction?

# CALCULATE background = yourFit(bin center) for all bins
noBkg=np.zeros(len(midpoints))

for j in range(m):
    for i in range(len(midpoints)):
        noBkg[i]=noBkg[i]+(midpoints[i]**(j))*fitparams[j]
for i in range(len(midpoints)):
    noBkg[i]=abs(counts[i]-noBkg[i])
        

plt.figure()
plt.step(x=midpoints, y=noBkg ,where='mid', alpha=.3 )
 
err=[[],[]]
for i in range(len(noBkg)):
    if math.sqrt(noBkg[i])<1:
        err[1].append(1)
    else:
        err[1].append( math.sqrt(noBkg[i]))
    err[0].append( math.sqrt(noBkg[i]))

plt.errorbar(midpoints, noBkg, yerr=err, fmt='none', ecolor='red')
plt.xlabel('GeV') 
plt.ylabel('Number of muon pairs') 

  
plt.title('Signal\n\n', fontweight ="bold") 
#plt.show()


# PLOT the background curve on top of your mass histogram (save it!)
# THINK: Are your estimated bkg values at all uncertain? 
##figure()


# EVALUATE signal = data - background
# THINK: What should you do if the background estimate is > data?
#        How could you find the uncertainty in data - background?




# PLOT the signal-only y with ERROR BARS
##figure()


# Great work! Save the data+background and signal-only plots as ANALYSIS
#             Save a DAY2 workspace!

## DAY 3 = CHARACTERIZATION of your signal 
# LOAD your DAY2 workspace

# EXTRACT the characteristics (mean, width, uncerts) of your signal y
# THINK: Which statistical distribution describes your signal y?
# TOOL: Curve-fitting app has a GUI to walk you through fitting, with 
#       automatic plotting and quality prints! Select x and y data and the 
#       function you"d like. (use Fit -> Save to Workspace) 
# SAVE a screencapture of your fit and its parameters.
y = np.delete(noBkg,[33,47])
y = np.delete(y,[0,16])
x = np.delete(midpoints,[33,47])
x = np.delete(x,[0,16])


# weighted arithmetic mean (corrected - check the section below)
mean = sum(x * y) / sum(y)
sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

def Gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

popt,pcov = curve_fit(Gaus,x,y,p0=[max(y),mean,sigma])
print(popt)
perr = np.sqrt(np.diag(pcov))
print(perr)
x=np.linspace(30,150,501).tolist()

plt.plot(x,Gaus(x,*popt),'r',label='fit')
plt.show()

# COMPARE: NSignal in signal y to NBackground under the y region
# THINK: how can you find the number of events in the signal y?
#        how can you find the number of bkg events under the y?
# PRINT: these values along with their uncertainties



# Almost done!
# THINK: Can you statistically distinguish signal from background?
#        Can you find this particle with a web search for you mass?
#        Research this particle (pdg.gov), find its width (capital Gamma)
#        Do your mass & width agree with the known values? Find percent
#        differences and also discrepancy/significance.
#        If your width is *much* larger than accepted, why might this be 



