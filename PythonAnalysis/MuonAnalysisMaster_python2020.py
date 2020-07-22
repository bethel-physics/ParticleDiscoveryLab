import math
import matplotlib.pyplot as plt
import pickle
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import pollsf

#Loads in data from a pickle file
data = pickle.load(open('DoubleMuParked_1M.pkl','rb'))

# This program begins by loading a pickle file that contains the  
# momentum and energy data in the form (E1, E2, Px1, Px2, Py1, Py2, Pz1, Pz2)
# from about 1M events. It produces a graph that has mass on the x-axis 
# and the number of muon pairs with that mass on the y-axis. 
# The user can choose the range of GeV values shown. 

## DAY 1 = RECONSTRUCTION

# Choose how many events to process
Ntoprocess=input("How many events would you like to process? ")
Ntoprocess = int(Ntoprocess)

# Initialize a vector that will hold 1 invariant mass per event
Masses = []
KineticEnergy = []

# User-chosen min/max values and resolution
Min=input("What is the min GeV? ")
Min = float(Min)
Max=input("What is the max GeV? ")
Max = float(Max)
n = input("Type in the number of bins: ")
n = int(n)
BinWidth = (Max - Min)/n

E=0
xMom=0
M=0
Px=0
Py=0
Pz=0
xK=0

# Loop over the number of events 
print("Looping over ", Ntoprocess, " events...", sep=" ")
for i in range(Ntoprocess):
   
    # COMPUTE the mass of particle X -> mu mu
    #E=E1+E2
    E = data[i][0] + data[i][1]

    #p=p1+p2 vector sum
    Px = data[i][2] + data[i][3]
    Py = data[i][4] + data[i][5]
    Pz = data[i][6] + data[i][7]

    # mc2 = sqrt( E^2 - p^2 ), where E and p are sums of the muons
    # GOAL for session 1 is implementing this correctly after intros
    M = math.sqrt(E**2 - Px**2 - Py**2 - Pz**2)
    Masses.append(M.real)

    # Calculate the Kinetic Energy of particle X. 
    # Store KE and mass values to plot later
    # Tip: make sure you mass value if "real" by using (massvalue).real
    if M> Min and M<Max:
        M=(M).real
        xK=E-M
        KineticEnergy.append(xK)


#end of the for loop

# HISTOGRAMMING -- create mass and KE histograms
# Draw a HISTOGRAM of counts versus mass
# HISTOGRAM needs: a list of values, the number of bins, the range, histtype and alpha are for visuals 
# plt.hist will return the y values of counts, the x values of the bin edges
# we won't use patches, but it gets mad if we don't include it
plt.figure()
massCounts, binEdg, patches =plt.hist(x=Masses, bins=n, range=(Min, Max), histtype='step',alpha=.5)

#Set up an array with the midpoints of the bins
center = (Max - Min)/n
midpoints = []
for i in range(len(binEdg)-1):
    midpoints.append(binEdg[i]+center*.5)

#Create a 2d array of error bar values one d for x, and one d for y
error=[[],[]]
for i in range(len(massCounts)):
    if massCounts[i]==0:
        error[1].append(1)
    else:
        error[1].append( math.sqrt(massCounts[i]))
    error[0].append( math.sqrt(massCounts[i]))

#Plot errorbars on top of histogram
#Needs: the center of the bins, the y values from the plt.hist function, the 2d error array
plt.errorbar(midpoints, massCounts, yerr=error, fmt='none', ecolor='red')

#Throw some labels on
plt.xlabel('GeV') 
plt.ylabel('Number of muon pairs') 
plt.title('MuonLab GeV/Counts\n\n', fontweight ="bold") 

#Show the plot
plt.show()


# Draw another HISTOGRAM of counts vs kinetic energy
plt.figure()
y, bins, patches =plt.hist(x=KineticEnergy, bins=n, range=(0, 800),histtype='step', log=True,alpha=.5)
center = (800)/n
centers=[]
for i in range(len(bins)-1):
    centers.append(bins[i]+center*.5)

error=[[],[]]
for i in range(len(y)):
    if y[i]==0:
        error[1].append(1)
    else:
        error[1].append( math.sqrt(y[i]))
    error[0].append( math.sqrt(y[i]))

plt.errorbar(centers, y, yerr=error, fmt='none', ecolor='red')

plt.xlabel('KE') 
plt.ylabel('Number of muon pairs') 
plt.title('MuonLab KE/Counts\n\n', fontweight ="bold") 
plt.show()



### Great work! SAVE these plots to represent your RAW DATA.

## DAY 2 = FITTING -- fit background on either side of the y

# CHOOSE mass values in GeV for where the y lies. 
minCut=input("What is the GeV to start the cut? ")
minCut=float(minCut)
maxCut=input("What is the GeV to end the cut? ")
maxCut=float(maxCut)


# REMOVE the y window completely from your list of: 
# bin edges, bin centers, counts, uncertainties. 
# This forms your BACKGROUND dataset
iterLow=0
iterHigh=0
for i in range(len(binEdg)):
    if binEdg[i]<minCut:
        iterLow=i
    if binEdg[i]<maxCut:
        iterHigh=len(binEdg)-i
add = len(massCounts) - iterHigh

# Set up lists to contain the background data
bkgCount=[]
bkgBin=[]
bkgErr=[[],[]]
bkgMids=[]

# Trim the data
for i in range(iterLow):
    bkgCount.append(massCounts[i])
    bkgBin.append(binEdg[i])
    bkgErr[1].append(error[1][i])
    bkgMids.append(midpoints[i])
for i in range(iterHigh):
    bkgCount.append(massCounts[i+add])
    bkgBin.append(binEdg[i+add])
    bkgErr[1].append(error[1][i+add])
    bkgMids.append(midpoints[i+add])
    
bkgBin.append(binEdg[len(binEdg)-1])

# Perform the polynominal fit
# pollsf needs: x-value centers, y-values, y uncertainties, N params
# pollsf returns: parameters, their uncerts, fitted y-values, chi^2
m=input("Enter the degree of polynomial: ")
m=int(m)
fitparams, uncerts, yvalues, chisq = pollsf.pollsf(bkgMids, bkgCount, bkgErr[1], m)

# printing chi2/degrees-of-freedom as a check (should be ~1, < 2)
print(chisq/(len(bkgMids)-m))


# Replot the histogram with the fit from pollsf on the same figure
plt.figure()

x, bins, patches=plt.hist(x=Masses, bins=n, range=(Min, Max), histtype='step',alpha=.5)
center = (Max - Min)/n
centers = []
for i in range(len(bins)-1):
    centers.append(bins[i]+center*.5)

error=[[],[]]
for i in range(len(x)):
    if x[i]==0:
        error[1].append(1)
    else:
        error[1].append( math.sqrt(x[i]))
    error[0].append( math.sqrt(x[i]))

plt.errorbar(centers, x, yerr=error, fmt='none', ecolor='red')
plt.xlabel('GeV') 
plt.ylabel('Number of muon pairs') 
plt.title('Fit Function\n\n', fontweight ="bold") 

# Plot the bkgMids and the yvalues from pollsf
plt.plot(bkgMids,yvalues)
plt.show()

# subtract the background from data

# set up an empty array to hold the signal data
noBkg=np.zeros(len(midpoints))

# loop through the number of exponents and calculate the y value of the fit for each x value in midpoints
for j in range(m):
    for i in range(len(midpoints)):
        noBkg[i]=noBkg[i]+(midpoints[i]**(j))*fitparams[j]
for i in range(len(midpoints)):
    # subtract the array of y values for the fit by the array of y values for the count
    noBkg[i]=abs(massCounts[i]-noBkg[i])
        

plt.figure()
# use a step plot to plot this signal
plt.step(x=midpoints, y=noBkg ,where='mid', alpha=.3 )
 
# plot error bars on top
error=[[],[]]
for i in range(len(noBkg)):
    if math.sqrt(noBkg[i])<1:
        error[1].append(1)
    else:
        error[1].append( math.sqrt(noBkg[i]))
    error[0].append( math.sqrt(noBkg[i]))

plt.errorbar(midpoints, noBkg, yerr=error, fmt='none', ecolor='red')
plt.xlabel('GeV') 
plt.ylabel('Number of muon pairs') 
plt.title('Signal\n\n', fontweight ="bold") 


# Great work! Save the data+background and signal-only plots as ANALYSIS

## DAY 3 = CHARACTERIZATION of your signal 

# set up arrays that only include the signal
# if there is more than one peak, make sure to include only the highest peak in these arrays
y = np.delete(noBkg,[add,len(noBkg)-1])
y = np.delete(y,[0,iterLow-1])
x = np.delete(midpoints,[add,len(noBkg)-1])
x = np.delete(x,[0,iterLow-1])


# calculate mean and sigma for the curve_fit function
mean = sum(x * y) / sum(y)
sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

# define the gausian function for curve_fit
def Gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

# call curve_fit, it takes in the definition of the funtion you would like to use,
# the x and y arrays of only the signal, and an array that includes the max value
# of y, and the mean and sigma values calculated earlier
popt,pcov = curve_fit(Gaus,x,y,p0=[max(y),mean,sigma])

# set up alot of xvales to plot the gausian on
xPlot=np.linspace(Min,Max,501).tolist()

# plot the gausian fit over the signal histogram
plt.plot(xPlot,Gaus(xPlot,*popt),'r',label='fit')
plt.show()

