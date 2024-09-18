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

## You will need ROOT to read the Open Data file -- we recommend the ROOT Docker container
## Information at: https://opendata.cern.ch/docs/cms-guide-docker
## Guide for analyzing NanoAOD is: https://opendata.cern.ch/docs/cms-getting-started-nanoaod
from ROOT import TFile, TLorentzVector, TH1D, TCanvas, TChain

## Prepare an object that will hold the tree called "Events" from the input files
t = TChain("Events")

## Information at: https://opendata.cern.ch/record/30555, 10.7483/OPENDATA.CMS.UZD7.Z50M
## You can also DOWNLOAD ROOT files from that webpage and open the local version!
filelist = [
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/127C2975-1B1C-A046-AABF-62B77E757A86.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/183BFB78-7B5E-734F-BBF5-174A73020F89.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/1BE226A3-7A8D-1B43-AADC-201B563F3319.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/1DE780E2-BCC2-DC48-815D-9A97B2A4A2CD.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/21DA4CE5-4E50-024F-9CE1-50C77254DD4E.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/2C6A0345-8E2E-9B41-BB51-DB56DFDFB89A.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/3676E287-A650-8F44-BBCB-3B8556966406.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/411A019C-7058-FD42-AD50-DE74433E6859.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/46A8960A-E58F-4648-9C12-2708FE7C12FB.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/4F0B53A7-6440-924B-AF48-B5B61D3CE23F.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/790F8A75-8256-3B46-8209-850DE0BE3C77.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/7F53D1DE-439E-AD48-871E-D3458DABA798.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/8A696857-C147-B04A-905A-F85FB76EDA23.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/8B253755-51F2-CB49-A4B6-C79637CAE23F.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/9528EA75-1C0B-9047-A9A3-6A47564F7A98.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/A6605227-0B58-864E-8422-B8990D18F622.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/B2DC29E0-8679-1D4F-A5AE-E7D0284A20D4.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/B450B2B3-BEF8-8C43-82BF-7AD0EF2EA7EA.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/B7AA7F04-5D5F-514A-83A6-9A275198852C.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/B93B57BF-4239-A049-9531-4C542C370185.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/C4558F81-9F2C-1349-B528-6B9DD6838D6D.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/C8CFC890-D4B8-8A4F-8699-C6ACCDF1620A.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/CAA285FF-7A12-F945-9183-DC7042178535.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/CD267D88-E57D-3B44-AC45-0712E2E12B87.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/E7C51551-7A75-5C41-B468-46FB922F36A9.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/EBC200F4-C06F-CE45-BAAA-7CAECDD3076F.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/EEB2FE3F-7CF3-BF4A-9F70-3F89FACE698E.root',
    'root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/2510000/F5E234F9-1E9C-0042-B395-AB6407E4A336.root'
]

## Add the tree from each file to the "chain"
for ifile in filelist:
    t.AddFile(ifile)

## Histogram for testing -- currently set to show Upsilon mesons
## Example resonance ranges:
##    2.8 -- 3.5 GeV (J/Psi meson),
##    8 -- 12 GeV (Upsilon meson),
##    30 -- 150 GeV (Z boson)
hist = TH1D("hist",";dimuon mass (GeV)",40,8,12)

## Initializations
posMuon1 = TLorentzVector()
negMuon1 = TLorentzVector()
data = []
isaved = 0
textfile = open("DoubleMuParked_100k.txt", "w")

## Loop over the events in the tree chain
for ievent in range(t.GetEntries()):

    ## Cap the number of events in your output (smaller file sizes! 100k - 200k is enough for the exercise)
    if isaved >= 200000: continue

    if ievent % 200000 == 0: print('Processed '+str(ievent)+' / '+str(t.GetEntries()))

    ## This line grabs one event's specific information
    t.GetEntry(ievent)

    ## To begin with, we have not identified either muon that we need
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

    ## Store this event's information into the text file for MATLAB analysis
    textfile.write(str(allmuons)[1:-2]+'\n')

    ## Store this event's information into the data object
    data.append(allmuons)
    isaved += 1    

## Write the data object into a text file (MATLAB analysis) and a pickle file (Python analysis)
textfile.close()
pickle.dump(data,open('DoubleMuon_2016H_200k.pkl','wb'))

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
    
        



    
