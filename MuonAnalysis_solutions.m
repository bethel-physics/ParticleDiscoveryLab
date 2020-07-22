%% SETUP
% Use the RUN SECTION button to advance through the code in chunks
% LOAD THE WORKSPACE on the command window
load('DoubleMuParked_100k.pkl')

%% DAY 1 = RECONSTRUCTION

% Choose how many events to process
Ntoprocess = input('How many events to process? ');

% Initialize a vector that will hold 1 invariant mass per event
Masses = [];
KineticEnergy = [];

% User-chosen min/max values and resolution
Min = input('Type in your min (in GeV): ');
Max = input('Type in your max (in GeV): ');
n = input('Type in the number of bins: ');
BinWidth = (Max - Min)/n;

% Loop over the number of events with at least 2 muons
disp(['Looping over ' num2str(Ntoprocess) ' events...']);
for i = 1:Ntoprocess
    if MomMarkers2(i) == 0
        continue; % this is not a good event
    end

    % COMPUTE the mass of particle X -> mu mu
    % SOLUTION: E = E1 + E2, p = p1 + p2 vector sum
    %           mc2 = sqrt( E^2 - p^2 )
    E = Energy{i}{MomMarkers1(i)} + Energy{i}{MomMarkers2(i)};
    Px = px{i}{MomMarkers1(i)} + px{i}{MomMarkers2(i)};
    Py = py{i}{MomMarkers1(i)} + py{i}{MomMarkers2(i)};
    Pz = pz{i}{MomMarkers1(i)} + pz{i}{MomMarkers2(i)};
    mc2 = sqrt( E^2 - Px^2 - Py^2 - Pz^2 );

    % THINK: Is this mass in your window from Min to Max? 
    %        What should you do if it's outside the window?
    % SOLUTION: discard those events

    if(mc2 > Min && mc2 < Max)
        
        % Calculate the Kinetic Energy of particle X. 
        % Store KE and mass values to plot later  
        % Tip: make sure you mass value if "real" by using real(massvalue)
        Masses = [Masses real(mc2)];  
        KineticEnergy = [KineticEnergy E - real(mc2)]; 
        % they'll try to solve for v! Redirect to Etotal - Erest = K
    end

end

%% HISTOGRAMMING -- create mass and KE histograms
% THINK: What do you expect your kinetic energy histogram to look like?
%        What do you expect your mass histogram to look like?
% SOLUTION: higher energies are always less probable, so falling from 0
%           mass is similar: falling from low -> high, but with a bump
%
% Vocab: imagine plot with 3 bins on x-axis: 0-10, 10-20, 20-30
% "Bin edges": 0, 10, 20, 30
% "Bin low edges": 0, 10, 20
% "Bin centers": 5, 15, 25 (want dots on plot to be here!)
% "Bin width": 10
MassEdges = Min:BinWidth:Max;  
MassCenters = MassEdges(1:length(MassEdges)-1) + 0.5*BinWidth; 
KEedges = 0:10:800; %chosen kind of arbitrarily, should fall from 0
KEcenters = KEedges(1:length(KEedges)-1) + 0.5*10; 

% Draw a HISTOGRAM of counts versus mass
% HISTOGRAM needs: Mass values list, Mass "bin edges" list 
figure(1)
histogram(real(Masses),MassEdges,'DisplayStyle','stairs');
xlabel('Mass [GeV]');
ylabel('Number of muon pairs');
hold on; % this lets you draw on top of this plot later
axis manual; % this keeps the axes fixed to where they are now

% THINK: What should the ERROR BARS be for each bin? 
%        What should you do if the bin has ZERO entries?
% SOLUTION: Error on N = sqrt(N)
%           Zero is not exact! Just lack of data. Use error_up = 1
%           But error bars can't dip below 0, that would be unphysical
% Tools: 
% HISTCOUNTS: gives counts, needs bin edges + input values (mass)
% ERRORBAR: draws dots+bars, needs bin centers, y values, down 
%           uncert list, up uncert list 
%[counts, xlow] = histcounts(Masses,MassEdges);
counts = histcounts(Masses, MassEdges);
disp(['Peak bin content: ' num2str(max(counts))]); % just to check

error_up = max(sqrt(counts),1); % avoid error of 0 going up
error_down = sqrt(counts);      % but error going down can be 0
errorbar(MassCenters,counts,error_down,error_up,'b.'); % b. = blue dots
hold off;

% Draw another HISTOGRAM of counts vs kinetic energy
% Add ERROR BARS
figure(2)
histogram(KineticEnergy,KEedges,'DisplayStyle','stairs'); 
set(gca,'YScale','log'); %log scale on the y axis
xlabel('KE [GeV]'); 
ylabel('Number of muon pairs');
hold on; % this lets you draw on top of this plot later
axis manual; % this keeps the axes fixed to where they are now

KEcounts = histcounts(KineticEnergy, KEedges);
KEerror_up = max(sqrt(KEcounts),1); % avoid error of 0 going up
KEerror_down = sqrt(KEcounts);      % but error going down can be 0
errorbar(KEcenters,KEcounts,KEerror_down,KEerror_up,'b.'); % b. = blue dots
hold off;

%%% Great work! SAVE these plots to represent your RAW DATA.
%%%             SAVE a DAY1 workspace (delete large E, px, py, pz)

%% DAY 2 = FITTING -- fit background on either side of the peak
% LOAD your DAY1 workspace

% Vocab: imagine a mass plot with a bump in the middle
% "Peak window": region along x-axis under the peak
% "background": smoothly falling slope of random events, 
%               including some of the events in the peak window
% "signal": events in the peak window minus the background

% CHOOSE mass values or bin numbers for where the peak lies. 
peakmin = input('Type in your peak min boundary (in GeV): ');
peakmax = input('Type in your peak max boundary (in GeV): ');
iMin = round((peakmin-Min)/BinWidth);
iMax = round((peakmax-Min)/BinWidth) + 1;

% REMOVE the peak window completely from your list of: 
% mass bin edges, mass bin centers, mass counts, mass uncertainties. 
% This forms your BACKGROUND dataset
bkgedges = [edges(1:iMin) edges(iMax:length(edges))];
bkgcounts = [counts(1:iMin) counts(iMax:length(counts))];
bkgcenters = [centers(1:iMin) centers(iMax:length(centers))];
bkgerror_up = [error_up(1:iMin) error_up(iMax:length(error_up))];
bkgerror_down = [error_down(1:iMin) error_down(iMax:length(error_down))];

% PERFORM a polynominal fit to the background
% THINK: Which type of curve do you expect will match your data best?
%        Imagine a curve connecting the two sides under your peak.
% SOLUTION: Probably a line, or 2rd/3rd order poly, likely not much higher
% Tools:
% POLLSF: gives fit params, uncerts, y-values, chi^2 value
%         needs bin centers, counts, uncertainties, N params
% CURVE FITTING APP: GUI to walk you through fitting, with automatic
%         plotting and quality prints! (use Fit -> Save to Workspace) 
%         Same as the "fit" command line function

% SOLUTION: pollsf method (actually curve-fitting app is simpler):
numpars = input('How many polynomial orders? 0, 1, 2, etc: ');
[params,paramerrs,fityvals,chisq] = pollsf(bkgcenters,bkgcounts,bkgerr,numpars);

% SOLUTION: curve-fitting app method (easier):
% Run their program up to this point
% Click on Apps in the menu and choose Curve Fitting
% x-axis data = bkgcenters, y-axis data = bkgcounts, choose "Polynomial"
% Choose an order of polynominal -- should fit and plot automatically
% "SSE" in the printout on the left of the GUI is the Chi^2 value 
% Fit->Save to Workspace will save the fit object 

% EVALUATE your fit by chi^2 and plotting
% -- Plotting: does the shape make any sense? 
% -- Chi^2 (or "SSE") is defined in Eq. 29. It describes the difference 
%    between the points and the fitted curve. LARGER chi^2 tends to mean 
%    more difference or scatter of points.
%    OPTIMALLY, Chi^2 / (# points - # parameters) is around 1
% REPEAT fitting until you are satisfied

% SOLUTION: pollsf method
disp(['Chi2/dof = ' num2str(chisq/(length(bkgcenters)-numpars))])
params = fliplr(params');
fittedbkgcounts = polyval(params,bkgcenters);

figure(3)
errorbar(bkgcenters,bkgcounts,bkgerror_down,bkgerror_up,'b.');
xlabel('Mass [GeV]');
ylabel('Number of muon pairs');
hold on;
plot(bkgcenters,fittedbkgcounts,'r-')
hold off;

%%% SUBTRACTION -- now you will subtract that background from data
%%% THINK: How will you estimate background in the signal peak window?
%%%        What do you expect the curve to look like after bkg subtraction?
% SOLUTION: evaluate the function at x-values inside the peak window
%           After subtraction should look like a ~Gaussian peak

% CALCULATE background = yourFit(bin center) for all bins
fittedcounts = polyval(params,centers); % more x-axis bins than before

% PLOT the background curve on top of your mass histogram (save it!)
% THINK: Are your estimated bkg values at all uncertain?
% SOLUTION: Yes, of course! But we have not discussed covariance and will
%           make the ~safe assumption that our background uncert is small
figure(4)
errorbar(MassCenters,counts,error_down,error_up,'b.');
xlabel('Mass [GeV]');
ylabel('Number of muon pairs');
hold on;
plot(MassCenters,fittedcounts,'r-')
hold off;

% EVALUATE signal = data - background
% THINK: What should you do if the background estimate is > data?
%        How could you find the uncertainty in data - background?
% SOLUTION: we cannot allow values < 0, so "floor" the subtraction
%           Uncert is like radioactivity: err = sqrt(errData^2 + errBkg^2)
sigcount = max(counts-fittedcounts,0);
sigerror_up = max(sqrt(counts+fittedcounts),1); % again, up error not 0
sigerror_down = sqrt(counts+fittedcounts);

% PLOT the signal-only peak with ERROR BARS
figure(3)
%plot(centers,sigcount,'b.'); %plot with just dots, draw errors below
errorbar(MassCenters,sigcount,sigerror_neg,sigerror_pos,'b.'); 
xlabel('Mass [GeV]');
ylabel('Number of muon pairs');

% Great work! Save the data+background and signal-only plots as ANALYSIS
%             Save a DAY2 workspace!

%% DAY 3 = CHARACTERIZATION of your signal 
% LOAD your DAY2 workspace

% EXTRACT the characteristics of your signal peak
% THINK: Which statistical distribution describes your signal peak?
% CURVE FITTING app: select x and y data and the function you'd like.
%                    Find MEAN and WIDTH with uncertainties in the output
% SOLUTION: they should do well with a Gaussian shape
%           They can do this without typing anything here, but should use
%           Fit->Save to Workspace. Example here names it "signalGaus"
gausparams = coeffvalues(signalGaus);
intervals = confint(signalGaus);
mu = gausparams(2);
mu_uncert = intervals(2,2) - mu;
sigma = gausparams(3);
sigma_uncert = intervals(3,2) - sigma;
disp(['Gaussian mean = ' num2str(mu) ' +/- ' num2str(mu_uncert)]);
disp(['Gaussian mean = ' num2str(sigma) ' +/- ' num2str(sigma_uncert)]);

% COMPARE: NSignal in signal peak to NBackground under the peak region
% THINK: how can you find the number of events in the signal peak?
%        how can you find the number of bkg events under the peak?
% PRINT: these values along with their uncertainties
% SOLUTION: NSignal = sum up counts from "sig counts"
%           NBackground = sum up counts from "fittedcounts" (iMin to iMax)
%           Of course, the groups are free to integrate their 2 fitted
%           functions from peakmin to peakmax!
bkginpeak = fittedcounts(iMin:iMax-1);
peakcounts = counts(iMin:iMax-1);
siginpeak = sigcount(iMin:iMax-1);

disp(['Number of background events in peak = ' num2str(sum(bkginpeak))...
     ' +/- ' num2str(sqrt(sum(bkginpeak)))]);
disp(['Number of signal events in peak = ' num2str(sum(siginpeak))...
     ' +/- ' num2str(sqrt(sum(peakcounts+bkginpeak)))]);
