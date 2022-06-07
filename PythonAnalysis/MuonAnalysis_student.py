### Welcome to the Particle Discovery Lab!

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
    



    ## Store mass and KE for events with mass inside your window




        
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
# Create a MASS histogram:

# Draw your mass histogram. Use plt.show() to draw your plot. 
# Be sure to save your y-axis values! 




# #### THINK: 
#  * What should the ERROR BARS be for each bin?                                      
#  * What should you do if the bin has ZERO entries?  
#  
# #### Tools:   plt.errorbar: 
# 
# plt.errorbar draws dots+bars when given x-axis bin centers, y-axis values, and up/down uncertainties. 
# Look up its drawing options: https://matplotlib.org/stable/gallery/statistics/errorbar_features                                                                                       

# Calculate lists of uncertainty values for plt.errorbar




# Define an array of bin centers




# Draw the new plot with error bars






#  #### Draw another HISTOGRAM with error bars of counts vs kinetic energy
# Get the y-axis values by drawing a new KE histogram




# Calculate the uncertainties





# Define an array of bin centers





# Draw the new plot with error bars






# #### Great work! 
# Save these plots to represent your raw data in your report. If you're using a jupyter notebook, save and checkpoint the notebook here. 

# ## Day 2 : Fitting
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






# REMOVE the peak window completely from your list of: 
# mass bin centers, mass counts, and mass uncertainties. 
# This forms your BACKGROUND dataset







# #### PERFORM a polynominal fit to the background
# #### THINK: 
# Which type of curve do you expect will match your data best? Imagine a curve connecting the two sides under your peak.
# 
# #### Tool: 
# The function *pollsf* is defined locally in pollsf.py.  Read pollsf.py to find information on the input and output parameters.
# 
# #### EVALUATE your fit:
#  * Plotting: does the shape make any sense? 
#  * Chi^2 is defined in "Place Holder". It describes the difference between the points and the fitted curve. LARGER chi^2 tends to mean more difference or scatter of points.
#  * OPTIMALLY, Chi^2 / (# points - # parameters) is around 1
# 
# #### REPEAT fitting until you are satisfied with both of these metrics

# Use pollsf to fit a polynomial
numpars = int(input('How many polynomial parameters? 1 (flat), 2 (line), etc: '))




# Print the chi2 metric described above





# Plot the fit on top of the background points
# Avoid using connecting lines between the points in your background fit




# ### SUBTRACTION -- now you will subtract that background from data
#
# In order to subtract the background contribution from your orginal data, you will need to estimate the background in your signal peak window.
#
# #### THINK: 
# How will you estimate background in the signal peak window? 

# Draw your background estimate on top of your full mass distribution




# #### THINK: 
# Are your estimated bkg values at all uncertain? 
# 
# How could you evaluate an uncertainty on the number of background events in each bin?
# 
# What do you expect the curve to look like after background subtraction?
#
# #### EVALUATE signal = data - background
# #### THINK: 
# Do you have any bins where the background estimate is larger than the data? What do you think about this situation? 
# 
# How could you find the uncertainty in data - background?

# Subtract background and plot the resulting signal with error bars





# #### Great work!
# Save the data+background and signal-only plots for the analysis section of your report. 

# ## Day 3 : Characterization
# Determine which particle you've discovered and use a fit to find its properties. 
# 
# #### EXTRACT the characteristics of your signal peak
# #### THINK: 
# Which statistical distribution describes your signal peak?
# 
# #### Tools: 
#  * A Gaussian function *Gaus* has been defined below. It takes x-axis values, an amplitude, a mean, and a width.
#  * The *curve_fit* function returns lists of fitted parameters and uncertainties when given a fit function, x and y-axis values, and initial conditions for the function's parameters. 
#  * Read about how to use this function at https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html

def Gaus(x,amplitude,mean,sigma):
    return amplitude*np.exp(-(x-mean)**2/(2*sigma**2))


* A Gaussian function *Gaus* has been defined below. It takes x-axis values, an amplitude, a mean, and a width.
#  * The *curve_fit* function returns lists of fitted parameters and uncertainties when given a fit function, x and y-axis values, and initial conditions for the function's parameters. 
#  * Read about how to use this function at https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
#  
# *SOLUTION: they should do well with a Gaussian shape. The initial conditions can be basic like p0=[1, (eyeballed peak center), 0.5]*

def Gaus(x,amplitude,mean,sigma):
    return amplitude*np.exp(-(x-mean)**2/(2*sigma**2))


# Use curve_fit to fit your signal peak using Gaus as the fit function




# Plot the fitted function on top of your signal distribution 
# xGaus below gives you lots of x-axis points to plot a smooth curve
xGaus = np.linspace(Min,Max,501).tolist()





# Print out the mean and width of your curve with uncertainties






# #### COMPARE: the number of signal events in signal peak window to the number of background events under the peak window.
# #### THINK: 
# How can you find the number of events in the signal peak? 
# 
# How can you find the number of bkg events under the peak?
# 
# #### PRINT: these values along with their uncertainties

# Print signal and background counts with uncertainties





# #### Almost done!
# #### THINK: 
#  * Can you statistically distinguish signal from background?
#  * Can you find this particle with a web search for you mass?
# 
# Research this particle (https://pdg.lbl.gov), find its width (capital Gamma). 
#  * Do your mass & width agree with the known values? 
#  * Find percent differences and also discrepancy/significance. 
#  * If your width is *much* larger than accepted, why might this be?
