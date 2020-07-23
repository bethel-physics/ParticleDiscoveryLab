import math, pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import pollsf


# Open the data file provided by your instructor
data = pickle.load(open('DoubleMuParked_100k.pkl','rb'))


# ## Day 1 : Reconstruction
# 
# Your first task is to load the CMS data file!
# Each data element has 8 pieces of information:
# 
# `E1`, `E2`, `px1`, `px2`, `py1`, `py2`, `pz1`, `pz2`
# 
# First choose a number of events to process, and the boundaries of your analysis window:

Ntoprocess = int(raw_input("How many events to process? "))
Min = float(raw_input("Type in your min (in GeV): "))
Max = float(raw_input("Type in your max (in GeV): "))
n = int(raw_input("How many x-axis bins would you like? "))

# Use these to compute a bin width
BinWidth = (Max - Min)/n

# let's get some empty objects ready for later
Masses = []
KineticEnergy = []


# Now we're ready to loop over the events in the data file and calculate the invariant mass of particle X. 
# 
# ### Think: 
#  * How will you use the 8 pieces of information to calculate the mass of X?
#  * How can you save only the events with a mass value inside your window?
#  * How can you calculate the relativistic kinetic energy of particle X? 
#  
# Write code to calculate the mass and KE of particle X. Store the results in Masses and KineticEnergy if the event has a mass inside your window.

print "Looping over",Ntoprocess,"events..."
for i in range(Ntoprocess):
    
    ## COMPUTE the mass of particle X that decays to 2 muons
    
    E = data[i][0] + data[i][1]  ## conserve E
    Px = data[i][2] + data[i][3] ## conserve px
    Py = data[i][4] + data[i][5] ## conserve py
    Pz = data[i][6] + data[i][7] ## conserve pz
    
    ## Invariant mass from E^2 = M^2 + |p|^2   (using natural units with c = 1!)
    M = math.sqrt(E**2 - Px**2 - Py**2 - Pz**2) 
    
    if M > Min and M < Max:
        Masses.append(M.real)
        
        KE = E - M  ## relativistic KE the easy way!
        KineticEnergy.append(KE)
        
print "Done!"


# ### HISTOGRAMMING -- create mass and KE histograms              
# 
# #### THINK: 
#  * What do you expect your kinetic energy histogram to look like?                   
#  * What do you expect your mass histogram to look like?                             
# 
# Make a quick sketch of what you expect for both plots                            
#                                                                                          
# #### Vocab: imagine plot with 3 bins on x-axis: 0-10, 10-20, 20-30                           
#  * "Bin edges": 0, 10, 20, 30                                                              
#  * "Bin centers": 5, 15, 25 (want dots on plot to be here!)                                
#  * "Bin width": 10  (you already have this for mass)
#  
# #### Tools: plt.hist
# plt.hist creates histograms when given a list of data, number of bins, and x-axis range. Look up its arguments and outputs!
# 
# Create a MASS histogram -- do you see a bump?

plt.figure(1)
massCounts, massEdges, patches = plt.hist(Masses,n,(Min,Max),histtype='step') ## want to keep the y-axis values 
plt.xlabel('dimuon mass (GeV)')
plt.ylabel('number of events')
plt.title('Mass validation')
plt.show()


# #### THINK: 
#  * What should the ERROR BARS be for each bin?                                      
#  * What should you do if the bin has ZERO entries?  
#  
# #### Tools:   plt.errorbar: 
# 
# plt.errorbar draws dots+bars when given x-axis bin centers, y-axis values, and up/down uncertainties. Look up its drawing options!                                                                                      

# Calculate the uncertainties
error = [[],[]]
error[0] = np.sqrt(massCounts)
error[1] = np.maximum(np.sqrt(massCounts),np.ones(len(massCounts)))

# Define an array of bin centers
massCenters = massEdges+BinWidth*0.5
massCenters = massCenters[:-1] # cut off the extra at the end

# Draw the new plot with error bars (k is black)
plt.errorbar(massCenters, massCounts, yerr=error, fmt='.k',ecolor='k')
plt.xlabel('dimuon mass (GeV)')
plt.ylabel('number of muon pairs')
plt.title('Mass')
plt.show()


#  #### Draw another HISTOGRAM with error bars of counts vs kinetic energy

plt.figure(2)
keCounts, keEdges, patches = plt.hist(KineticEnergy,n,(0,800),histtype='step',log=True) ## want to keep the y-axis values 

# Calculate the uncertainties
keerror = [[],[]]
keerror[0] = np.sqrt(keCounts)
keerror[1] = np.maximum(np.sqrt(keCounts),np.ones(len(keCounts)))

# Define an array of bin centers
keCenters = keEdges+(800/n)*0.5
keCenters = keCenters[:-1] # cut off the extra at the end

# Draw the new plot with error bars (k is black)
plt.errorbar(keCenters, keCounts, yerr=keerror, fmt='.k',ecolor='k')
plt.xlabel('kinetic energy (GeV)')
plt.ylabel('number of muon pairs')
plt.title('Kinetic Energy')
plt.show()


# #### Great work! 
# Save these plots to represent your raw data in your report. If you're using a jupyter notebook, save and checkpoint the notebook here. 

# ## Day 2: Fitting
# Fit the background on either side of the signal peak in your mass distribution. 
# 
# #### Vocab: imagine a mass plot with a bump in the middle
#  * "Peak window": region along x-axis under the peak
#  * "background": smoothly falling slope of random events, including some of the events in the peak window
#  * "signal": events in the peak window minus the background

# Choose mass values or bin numbers for where the peak lies
peakmin = float(input('Enter your peak minimum (in GeV)'))
peakmax = float(input('Enter your peak maximum (in GeV) '))

# Convert these mass values to bin numbers
iMin = int(round((peakmin-Min)/BinWidth))
iMax = int(round((peakmax-Min)/BinWidth))
print iMin,iMax

# REMOVE the peak window completely from your list of: 
# mass edges, mass centers, mass counts, mass uncertainties. 
# This forms your BACKGROUND dataset
bkgEdges = massEdges[0:iMin].tolist() + massEdges[iMax:-1].tolist()
bkgCounts = massCounts[0:iMin].tolist() + massCounts[iMax:-1].tolist()
bkgCenters = massCenters[0:iMin].tolist() + massCenters[iMax:-1].tolist()
bkgError = [[],[]]
bkgError[0] = error[0][0:iMin].tolist() + error[0][iMax:-1].tolist()
bkgError[1] = error[1][0:iMin].tolist() + error[1][iMax:-1].tolist()

# Check: edges should have 1 more than the others
print len(bkgEdges),len(bkgCounts),len(bkgCenters),len(bkgError[0])


# #### PERFORM a polynominal fit to the background
# #### THINK: 
# Which type of curve do you expect will match your data best? Imagine a curve connecting the two sides under your peak.
# 
# *SOLUTION: Probably a line, or 2rd/3rd order poly, likely not much higher*
# 
# #### Tool: 
# The function *pollsf* gives fit params, uncerts, y-values, chi^2 value. Needs bin centers, counts, uncertainties, M parameters in the polynominal.
# 
# #### EVALUATE your fit:
#  * Plotting: does the shape make any sense? 
#  * Chi^2 is defined in Eq. 29. It describes the difference between the points and the fitted curve. LARGER chi^2 tends to mean more difference or scatter of points.
#  * OPTIMALLY, Chi^2 / (# points - # parameters) is around 1
# 
# #### REPEAT fitting until you are satisfied

# Use pollsf to fit a polynomial
numpars = int(input('How many polynomial parameters? 1 (flat), 2 (line), etc: '))
params,paramerrs,fityvals,chisq = pollsf.pollsf(bkgCenters,bkgCounts,bkgError[1],numpars)

# Print the chi2 metric
print chisq/(len(bkgCenters) - numpars)

# Plot the fit on top of the background points
plt.figure(3)
plt.errorbar(bkgCenters, bkgCounts, yerr=bkgError, fmt='.k', ecolor='k')
plt.plot(bkgCenters,fityvals,'b-')
plt.xlabel('dimuon mass (GeV)')
plt.ylabel('number of muon pairs')
plt.title('Fit validation')
plt.show()


# ### SUBTRACTION -- now you will subtract that background from data
# #### THINK: 
# How will you estimate background in the signal peak window? What do you expect the curve to look like after bkg subtraction?
# 
# *SOLUTION: evaluate the function at x-values inside the peak window. After subtraction should look like a ~Gaussian peak*

# Draw your background estimate on top of your full mass distribution
paramlist = params.tolist()
paramlist.reverse() # polyval wants a flipped order of parameters
fittedCounts = np.polyval(paramlist,massCenters)

plt.figure(4)
plt.errorbar(massCenters, massCounts, yerr=error, fmt='.k', ecolor='k')
plt.plot(massCenters,fittedCounts,'b-')
plt.xlabel('dimuon mass (GeV)')
plt.ylabel('number of muon pairs')
plt.title('Data and background estimate')
plt.show()


# #### THINK: 
# Are your estimated bkg values at all uncertain? How could you evaluate an uncertainty on the number of background events in each bin?
# 
# *SOLUTION: Yes, of course! But we have not discussed covariance, so we won't attempt to use the fitted parameters' uncertainties. You could estimate an uncertainty of 0, or use the Poisson sqrt(N) formula as shown here.*  
# 
# #### EVALUATE signal = data - background
# #### THINK: 
# What should you do if the background estimate is > data? How could you find the uncertainty in data - background?
# 
# *SOLUTION: values less than 0 are not unphysical anymore, since this is an "estimate". Uncertainty must be propagated: sqrt(errData^2 + errBkg^2)*

# Subtract background and plot the resulting signal with error bars
signalCounts = massCounts - fittedCounts
signalErrors = [[],[]]
signalErrors[0] = np.sqrt(signalCounts + fittedCounts)
signalErrors[1] = np.maximum(np.sqrt(signalCounts + fittedCounts),np.ones(len(signalCounts)))

plt.errorbar(massCenters,signalCounts,yerr=signalErrors,fmt='.r',ecolor='r')
plt.xlabel('dimuon mass (GeV)')
plt.ylabel('number of events')
plt.title('Signal = Data - Background')
plt.show()


# #### Great work!
# Save the data+background and signal-only plots for the analysis section of your report. 

# ## Day 3: Characterization
# Determine which particle you've discovered and use a fit to find its properties. 
# 
# #### EXTRACT the characteristics of your signal peak
# #### THINK: 
# Which statistical distribution describes your signal peak?
# 
# #### Tools: 
#  * A Gaussian function *Gaus* has been defined below. It takes x-axis values, an amplitude, a mean, and a width.
#  * The *curve_fit* function returns lists of fitted parameters and uncertainties when given a fit function, x and y-axis values, and initial conditions for the function's parameters.
#  
# *SOLUTION: they should do well with a Gaussian shape. The initial conditions can be basic like p0=[1, (eyeballed peak center), 0.5]*

def Gaus(x,amplitude,mean,sigma):
    return amplitude*np.exp(-(x-mean)**2/(2*sigma**2))


# Use curve_fit to fit your signal peak using Gaus
gausParams,gausUncerts = curve_fit(Gaus,massCenters,signalCounts,p0=[1,3,0.1])
print gausParams


# Plot the fitted function on top of your signal distribution 
# (tip: use many more x-axis values for smoothness!)
xGaus = np.linspace(Min,Max,501).tolist()

plt.figure(5)
plt.errorbar(massCenters, signalCounts, yerr=signalErrors, fmt='r.',ecolor='r')
plt.plot(xGaus, Gaus(xGaus,*gausParams), 'b')
plt.xlabel('dimuon mass (GeV)')
plt.ylabel('number of events')
plt.title('Fitted signal')
plt.show()

# Print out the mean and width of your curve with uncertainties
print "Mean =",gausParams[1],"+/-",gausUncerts[1][1]
print "Width =",abs(gausParams[2]),"+/-",gausUncerts[2][2]


# #### COMPARE: NSignal in signal peak to NBackground under the peak region
# #### THINK: 
# How can you find the number of events in the signal peak? How can you find the number of bkg events under the peak?
# 
# #### PRINT: these values along with their uncertainties
# 
# *SOLUTION: NSignal = sum up counts from "sig counts". NBackground = sum up counts from "fittedcounts" (iMin to iMax). Of course, the groups are free to integrate their 2 fitted functions from peakmin to peakmax! They should be able to show that sqrt(sumtotal) is the proper error propagation across the sum over bins.*

bkginpeak = sum(fittedCounts[iMin:iMax])
siginpeak = sum(signalCounts[iMin:iMax])
print 'NBkg =',bkginpeak,'+/-',math.sqrt(bkginpeak)
print 'NSig =',siginpeak,'+/-',math.sqrt(siginpeak)


# #### Almost done!
# #### THINK: 
#  * Can you statistically distinguish signal from background?
#  * Can you find this particle with a web search for you mass?
# 
# Research this particle (pdg.gov), find its width (capital Gamma). 
#  * Do your mass & width agree with the known values? 
#  * Find percent differences and also discrepancy/significance. 
#  * If your width is *much* larger than accepted, why might this be?
#  
# *SOLUTION: for the Z boson their values should agree well -- mass should agree at least within +/- width! For light mesons they will find that the width they measure is much much too large, because of limited detector resolution. The lifetimes of these mesons can't be determined from this mass plot.*
