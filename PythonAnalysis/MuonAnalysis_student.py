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

# Loading data from the .pkl file
data = pickle.load(open('DoubleMuParked_100k.pkl','rb'))

# This program begins by loading a pkl file that contains the  
# momentum and energy data from about 1M events. It produces a 
# graph that has mass on the x-axis and the number of muon pairs 
# with that mass on the y-axis. 
# The user can choose the range of GeV values shown. 

## DAY 1 = RECONSTRUCTION

# Choose how many events to process
Ntoprocess = input("How many events to process? ")
Ntoprocess = int(Ntoprocess)

# Initialize a vector that will hold 1 invariant mass per event
Masses = []
KineticEnergy = []

# User-chosen min/max values and resolution
Min = input("Type in your min (in GeV): ")
Min = float(Min)
Max = input("Type in your max (in GeV): ")
Max = float(Max)
n = input("Type in the number of bins: ")
n = int(n)
BinWidth = (Max - Min)/n

# Loop over the number of events with at least 2 muons
print("Looping over ", Ntoprocess, " events...", sep=" ")
for i in range(Ntoprocess):
   
    # COMPUTE the mass of particle X -> mu mu
    

    # THINK: Is this mass in your window from Min to Max? 
    #        What should you do if it"s outside the window?


    
    # Calculate the Kinetic Energy of particle X. 
    # Store KE and mass values to plot later
    # Tip: make sure you mass value if "real" by using real(massvalue)
    
    
    


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
# HISTOGRAM needs: a list of values, the number of bins, and the range
# plt.hist will return the y values of counts, the x values of the bin edges, and a silent list of individual patches
plt.figure()



# THINK: What should the ERROR BARS be for each bin? 
#        What should you do if the bin has ZERO entries?
# Tools: 
# plt.errorbars: draws dots+bars, needs bin centers, y values, 2d array of down 
#           uncert list/up uncert list 


plt.show()

# Draw another HISTOGRAM of counts vs kinetic energy
# Add ERROR BARS
plt.figure()


plt.show()



### Great work! SAVE these plots to represent your RAW DATA.

## DAY 2 = FITTING -- fit background on either side of the y

# Vocab: imagine a mass plot with a bump in the middle
# "Peak window": region along x-axis under the y
# "background": smoothly falling slope of random events, 
#               including some of the events in the y window
# "signal": events in the y window minus the background

# CHOOSE mass values in GeV for where the y lies. 



# REMOVE the peak window completely from your list of: 
# bin edges, bin centers, counts, uncertainties. 
# This forms your BACKGROUND dataset




# PERFORM a polynominal fit to the background
# THINK: Which type of curve do you expect will match your data best?
#        Imagine a curve connecting the two sides under your y.
# Tool: POLLSF gives fit params, uncerts, y-values, chi^2 value
#       needs bin centers, counts, up uncertainties, N params
#       This is a least-squares fitter that uses uncerts!



# EVALUATE your fit by chi^2 and plotting
# -- Plotting: does the shape make any sense? Make a helpful plot
# -- Chi^2 (or "SSE") is defined in Eq. 29. It describes the difference 
#    between the points and the fitted curve. LARGER chi^2 tends to mean 
#    more difference or scatter of points.
#    OPTIMALLY, Chi^2 / (# points - # parameters) is around 1
# REPEAT fitting until you are satisfied
plt.figure()

plt.show()



## SUBTRACTION -- now you will subtract that background from data
# THINK: How will you estimate background in the signal y window?
#        What do you expect the curve to look like after bkg subtraction?

# CALCULATE background = yourFit(bin center) for all bins
        




# PLOT the background curve on top of your mass histogram (save it!)
# THINK: Are your estimated bkg values at all uncertain? 
plt.figure()

plt.show()


# EVALUATE signal = data - background
# THINK: What should you do if the background estimate is > data?
#        How could you find the uncertainty in data - background?




# PLOT the signal-only y with ERROR BARS
plt.figure()

plt.show()


# Great work! Save the data+background and signal-only plots as ANALYSIS
#             Save a DAY2 workspace!

## DAY 3 = CHARACTERIZATION of your signal 
# LOAD your DAY2 workspace

# EXTRACT the characteristics (mean, width, uncerts) of your signal y
# THINK: Which statistical distribution describes your signal y?
# TOOL: curve_fit is a function that takes in the function you want to fit, x and y data,
#       and an array in the form p0=[max(y),mean,uncerts]
# SAVE a screencapture of your fit and its parameters.
plt.figure()

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



