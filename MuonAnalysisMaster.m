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
    
    % E = E1 + E2
    E = Energy{i}{MomMarkers1(i)} + Energy{i}{MomMarkers2(i)};
    
    % p = p1 + p2 vector sum
    Px = px{i}{MomMarkers1(i)} + px{i}{MomMarkers2(i)};
    Py = py{i}{MomMarkers1(i)} + py{i}{MomMarkers2(i)};
    Pz = pz{i}{MomMarkers1(i)} + pz{i}{MomMarkers2(i)};
        
    % mc2 = sqrt( E^2 - p^2 ), where E and p are sums of the muons
    % GOAL for session 1 is implementing this correctly after intros
    mc2 = sqrt( E^2 - Px^2 - Py^2 - Pz^2 );
    
    % Store mass and kinetic energy for events in the chosen window
    if(mc2 > Min & mc2 < Max)
        Masses(i) = mc2;
        KineticEnergy(i) = E - mc2;
    end
end

%% Histogramming -- create mass and KE histograms

% create a vector of x-value bin edges [Min Min+BinWidth ... Max]
edges = Min:BinWidth:Max;  % 0, 10, 20, 30...400
KEedges = 0:10:800; 
% make a vector of center x-values, shift low edge up by width/2
centers = edges(1:length(edges)-1) + 0.5*BinWidth; 

% Draw the histogram of kinetic energy
% histogram needs the mass values list and x-value edges list
figure(1);
histogram(KineticEnergy,KEedges); 
set(gca,'YScale','log'); %log scale on the y axis
xlabel('KE [GeV]'); 
ylabel('Count');

% Draw the histogram of counts versus mass
% real( ) is to protect against it yelling about a complex mass value...
figure(2);
histogram(real(Masses),edges,'DisplayStyle','stairs');
xlabel('Mass [GeV]');
ylabel('Number of muon pairs');
hold on; % this lets you draw on top of this plot later
axis manual; % this keeps the axes fixed to where they are now

% Get the y-values (counts) in each bin to find the uncertainty
% histcounts needs the mass value list and the x-value edges list
[counts, xlow] = histcounts(real(Masses),edges);
disp(['Peak bin content: ' num2str(max(counts))]);

% Positive uncertainty should not be 0, negative uncertainty can be
counterrpos = max(sqrt(counts),1);
counterrneg = sqrt(counts);

% Draw the error bars on top of the histogram
% errorbar needs: x-value centers, y-values, down uncert, up uncert
errorbar(centers,counts,counterrneg,counterrpos,'b.'); % b. = blue dots

%% Fit the background on either side of the peak -- SESSION 2

% choose mass values in GeV for where the peak lies. 
peakmin = input('Type in your peak min boundary (in GeV): ');
peakmax = input('Type in your peak max boundary (in GeV): ');

% make a vector of background masses "less than min OR greater than max"
%BkgMass = Masses(find(Masses < peakmin | Masses > peakmax));

% make a new list of x-value edges for only background bins
% (peak-Min)/BinWidth find the index of the cutoff (ex: bin 15 of 50)
peakminIndex = round((peakmin-Min)/BinWidth);
peakmaxIndex = round((peakmax-Min)/BinWidth) + 1;
bkgedges = [edges(1:peakminIndex) edges(peakmaxIndex:length(edges))];
% Shift by width/2 to get centers
bkgcenters = bkgedges(1:length(bkgedges)-1) + 0.5*BinWidth;

% trim the list of y-values (counts) as well
bkgcounts = [counts(1:peakminIndex) counts(peakmaxIndex:length(counts))];
% Uncertainty is at least 1
bkgerr = max(sqrt(bkgcounts),1);

% Perform the polynominal fit
% pollsf needs: x-value centers, y-values, y uncertainties, N params
% pollsf returns: parameters, their uncerts, fitted y-values, chi^2
numpars = input('How many polynomial orders? 0, 1, 2, etc: ');
[params,paramerrs,fityvals,chisq] = pollsf(bkgcenters,bkgcounts,bkgerr,numpars);

% Printing chi2/degrees-of-freedom as a check (should be ~1, < 2)
disp(['Chi2/dof = ' num2str(chisq/(length(bkgcenters)-numpars))])

%% Subtract fitted background from the data to find signal

% Get counts from the fitted function over the full range
% pollsf gives a column, fliplr(params') gives a row in the opposite order
params = fliplr(params');
fittedcounts = polyval(params,centers);

% Plot the fit on top of the histogram as a red line
plot(centers,fittedcounts,'r-')

% Signal = Data - background (not allowed to be negative)
sigcount = max(counts-fittedcounts,0);

% Uncertainty on signal count -- sqrt(NTot + NBkg)
sigcounterrpos = max(sqrt(counts+fittedcounts),1);
sigcounterrneg = sqrt(counts+fittedcounts);
% Stop the uncertainty from going below zero counts
for ibin = 1:length(sigcounterrneg)
    if sigcounterrneg(ibin) > sigcount(ibin)
        sigcounterrneg(ibin) = sigcount(ibin);
    end
end

% Plot the signal-only peak on top of all the others as green circles
figure(3);
plot(centers,sigcount,'b.'); %plot with just dots, draw errors below
xlabel('Mass [GeV]');
ylabel('Number of muon pairs');
hold on; axis manual; 
errorbar(centers,sigcount,sigcounterrneg,sigcounterrpos,'b.'); %errorbars

%% Find the number of signal and background under the peak -- SESSION 3?

% extract the peak region of the background counts
bkginpeak = fittedcounts(peakminIndex:peakmaxIndex-1);
disp(['Number of background events in peak = ' num2str(sum(bkginpeak))...
     ' +/- ' num2str(sqrt(sum(bkginpeak)))]);

% extract the peak region of the full counts 
peakcounts = counts(peakminIndex:peakmaxIndex-1);

% Signal = data - background
siginpeak = peakcounts - bkginpeak;

% Uncertainty = sqrt(Ntotal + Nbkg)
disp(['Number of signal events in peak = ' num2str(sum(siginpeak))...
     ' +/- ' num2str(sqrt(sum(peakcounts+bkginpeak)))]);

%% Statistics of the peak

% No quick fix for mean value after we do Data - Bkg...
% We'll loop over the peak bins and make a list of their contents 
peakcenters = centers(peakminIndex:peakmaxIndex-1);
sigmasses = [];
for ibin = 1:length(peakcenters)
    % if the bin from 60-62 GeV has 5 counts we want a list like
    % [61 61 61 61 61]
    temp = peakcenters(ibin)*ones(1,round(peakcounts(ibin)));
    
    % add items to the list each time
    sigmasses = [sigmasses temp];
end
disp(['Mean value = ' num2str(mean(sigmasses))]);
disp(['Std dev = ' num2str(std(sigmasses))]);
    
%% DONE! (well, didn't print uncertainties on yields ;)
