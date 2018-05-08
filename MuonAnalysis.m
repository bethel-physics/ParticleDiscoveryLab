% This program begins by loading a workspace that contains the momentum 
% and energy data from about 1M events. It produces a graph that has
% mass on the x-axis and the number of muon pairs with that mass on the
% y-axis. The user can choose the range of GeV values shown. 

%% Loop over events and calculate invariant mass, kinetic energy
% FIRST LOAD THE WORKSPACE on the command window

% Initialize a vector that will hold 1 invariant mass per event
Masses(1:length(px)) = 0;
KineticEnergy(1:length(px)) = 0;

% User-chosen min/max values and resolution
Min = input('Type in your min (in GeV): ');
Max = input('Type in your max (in GeV): ');
n = input('Type in the number of bins: ');
BinWidth = (Max - Min)/n;

% Loop over the number of events with at least 2 muons
disp(['Looping over ' num2str(length(px)) ' events...']);
for i = 1:length(px)
    if MomMarkers2(i) == 0
        continue; % this is not a good event
    end

    % Compute mass and kinetic energy
    
    % Store mass and kinetic energy for events in the chosen window

end

%% Histogramming -- create mass and KE histograms

% make vectors of x-value bin edges for mass, KE[Min Min+BinWidth ... Max]


% make vectors of center x-values, shift low edge up by width/2


% Draw the histogram of counts vs kinetic energy
% histogram needs the KE values list and KE x-value edges list
figure(1)


% Draw the histogram of counts versus mass 
figure(2)


% Get the y-values (counts) in each bin to find the uncertainty
% histcounts needs the mass value list and the x-value edges list


% Find uncertainties to use in the errorbar method
% Positive uncertainty should not be 0, negative uncertainty can be


% Draw the error bars on top of the histogram
% errorbar needs: x-value centers, y-values, down uncert, up uncert



%% Fit the background on either side of the peak

% choose mass values in GeV for where the peak lies. 


% trim the lists of x-value edges and centers for only background bins


% trim the list of y-values (counts) and find uncerts in each bkg bin


% Perform the polynominal fit
% pollsf needs: x-value centers, y-values, y uncertainties, N params
% pollsf returns: parameters, their uncerts, fitted y-values, chi^2


% Print chi2/degrees-of-freedom as a check (should be ~1, < 2)


%% Subtract fitted background from the data to find signal

% Get counts from the fitted function over the full range


% Plot the fit on top of the histogram


% Signal = Data - background (not allowed to be negative)


% Uncertainty on signal count -- sqrt(NTot + NBkg)


% Plot the signal-only peak on a new graph with errorbars
figure(3)



%% Find the number of signal and background under the peak

% extract the peak region of the background counts


% extract the peak region of the signal counts (label it siginpeak)


% Print N background & N signal under peak with uncertainties



%% Statistics of the peak

% No quick fix for mean value after we do Data - Bkg...
% We'll loop over the peak bins and make a list of their contents 
peakcenters = 0; %make a list here of the bin centers
sigmasses = [];
for ibin = 1:length(peakcenters)
    % if the bin from 60-62 GeV has 5 counts we want a list like
    % [61 61 61 61 61]
    temp = peakcenters(ibin)*ones(1,round(siginpeak(ibin)));
    
    % add items to the list each time
    sigmasses = [sigmasses temp];
end
    
% Print the mean and standard deviation of sigmasses

