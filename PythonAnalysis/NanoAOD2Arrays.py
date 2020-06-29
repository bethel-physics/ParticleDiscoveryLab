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

import os,sys
from array import array
import pickle
import math

## You will need ROOT to read the Open Data file
## Information at: https://root.cern.ch/  and http://opendata.cern.ch/docs/cms-getting-started-2011
from ROOT import TFile, TLorentzVector, TH1D, TCanvas

## Information at: http://opendata.cern.ch/record/12341,   DOI:10.7483/OPENDATA.CMS.LVG5.QT81
## You can also DOWNLOAD the ROOT file from that webpage and open the local version!
NanoMuons = TFile.Open("root://eospublic.cern.ch//eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/Run2012BC_DoubleMuParked_Muons.root")

## Histogram for testing -- currently set to show Upsilon mesons
## Example resonance ranges: 2.8 -- 3.5 GeV (J/Psi meson), 8 -- 12 GeV (Upsilon meson), 30 -- 150 GeV (Z boson)
hist = TH1D("hist",";dimuon mass (GeV)",40,8,12)

## Read the TTree object from the file
t = NanoMuons.Get("Events")

## Initializations
posMuon1 = TLorentzVector()
negMuon1 = TLorentzVector()
data = []
isaved = 0

## Loop over the events in the TTree
for ievent in xrange(t.GetEntries()):

    ## Cap the number of events in your output (smaller file sizes! 100k is enough for the exercise)
    if isaved >= 100000: continue

    if ievent % 100000 == 0: print 'Processed',ievent,' / ',t.GetEntries()

    ## This line grabs one event's specific information
    t.GetEntry(ievent)

    negMu1 = -1
    posMu1 = -1

    ## Save the highest pT positive and negative muons (list goes in pT-decreasing order)
    for imu in range(t.nMuon):
        if t.Muon_charge[imu] < 0:  
            if negMu1 < 0: negMu1 = imu
        else:
            if posMu1 < 0: posMu1 = imu

    ## Skip events without a +/- charge pair
    if posMu1 < 0 or negMu1 < 0: continue 

    ## Fill these chosen muons into the TLorentzVectors
    posMuon1.SetPtEtaPhiM(t.Muon_pt[posMu1],t.Muon_eta[posMu1],t.Muon_phi[posMu1],t.Muon_mass[posMu1])
    negMuon1.SetPtEtaPhiM(t.Muon_pt[negMu1],t.Muon_eta[negMu1],t.Muon_phi[negMu1],t.Muon_mass[negMu1])

    ## Create a list of physics quantities
    allmuons = [posMuon1.E(), negMuon1.E(), posMuon1.Px(), negMuon1.Px(), posMuon1.Py(), negMuon1.Py(), posMuon1.Pz(), negMuon1.Pz()]

    ## Test it! This is just for you
    Z = posMuon1+negMuon1
    hist.Fill(Z.M())

    ## Store this event's information into the data object
    data.append(allmuons)
    isaved += 1    

## Write the data object into a pickle file for students to use
pickle.dump(data,open('DoubleMuParked_100k.pkl','wb'))

## Test drawing your histogram
hist.Draw()

## Test a mini-replication of the student's future code
## The matplotlib plot *should* match your ROOT histogram
masses = []
for i in range(len(data)):
    E = data[i][0] + data[i][1]
    Px = data[i][2] + data[i][3]
    Py = data[i][4] + data[i][5]
    Pz = data[i][6] + data[i][7]
    M = math.sqrt(E**2 - Px**2 - Py**2 - Pz**2)
    masses.append(M.real)

import matplotlib.pyplot as plt
from numpy import linspace

edges = linspace(8,12,41).tolist()

n = plt.hist(masses,edges)
plt.show()
    
        



    
