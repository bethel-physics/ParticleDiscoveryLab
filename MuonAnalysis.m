% This program begins by loading a workspace that contains the  
% momentum and energy data from about 1M events. It produces a 
% graph that has mass on the x-axis and the number of muon pairs 
% with that mass on the y-axis. 
% The user can choose the range of GeV values shown. 

%% DAY 1 = RECONSTRUCTION
% LOAD THE WORKSPACE on the command window
% Use the RUN SECTION button to advance through the code in chunks

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

    % COMPUTE the mass of particle X -> mu mu


    % THINK: Is this mass in your window from Min to Max? 
    %        What should you do if it's outside the window?


    
    % Calculate the Kinetic Energy of particle X. 
    % Store KE and mass values to plot later
    % Tip: make sure you mass value if "real" by using real(massvalue)



end

% HISTOGRAMMING -- create mass and KE histograms
% THINK: What do you expect your kinetic energy histogram to look like?
%        What do you expect your mass histogram to look like?
%        Make a quick sketch of what you expect for both plots
%
% Vocab: imagine plot with 3 bins on x-axis: 0-10, 10-20, 20-30
% "Bin edges": 0, 10, 20, 30
% "Bin centers": 5, 15, 25 (want dots on plot to be here!)
% "Bin width": 10



% Draw a HISTOGRAM of counts versus mass
% HISTOGRAM needs: Mass values list, Mass "bin edges" list 
figure()


% THINK: What should the ERROR BARS be for each bin? 
%        What should you do if the bin has ZERO entries?
% Tools: 
% HISTCOUNTS: gives counts, needs bin edges + input values (mass)
% ERRORBAR: draws dots+bars, needs bin centers, y values, down 
%           uncert list, up uncert list 




% Draw another HISTOGRAM of counts vs kinetic energy
% Add ERROR BARS
figure()



%%% Great work! SAVE these plots to represent your RAW DATA.
%%%             SAVE a DAY1 workspace (delete large E, px, py, pz)

%% DAY 2 = FITTING -- fit background on either side of the peak
% LOAD your DAY1 workspace

% Vocab: imagine a mass plot with a bump in the middle
% "Peak window": region along x-axis under the peak
% "background": smoothly falling slope of random events, 
%               including some of the events in the peak window
% "signal": events in the peak window minus the background

% CHOOSE mass values in GeV for where the peak lies. 



% REMOVE the peak window completely from your list of: 
% bin edges, bin centers, counts, uncertainties. 
% This forms your BACKGROUND dataset




% PERFORM a polynominal fit to the background
% THINK: Which type of curve do you expect will match your data best?
%        Imagine a curve connecting the two sides under your peak.
% Tool: POLLSF gives fit params, uncerts, y-values, chi^2 value
%       needs bin centers, counts, uncertainties, N params
%       This is a least-squares fitter that uses uncerts!



% EVALUATE your fit by chi^2 and plotting
% -- Plotting: does the shape make any sense? Make a helpful plot
% -- Chi^2 (or "SSE") is defined in Eq. 29. It describes the difference 
%    between the points and the fitted curve. LARGER chi^2 tends to mean 
%    more difference or scatter of points.
%    OPTIMALLY, Chi^2 / (# points - # parameters) is around 1
% REPEAT fitting until you are satisfied
figure()



%% SUBTRACTION -- now you will subtract that background from data
% THINK: How will you estimate background in the signal peak window?
%        What do you expect the curve to look like after bkg subtraction?

% CALCULATE background = yourFit(bin center) for all bins



% PLOT the background curve on top of your mass histogram (save it!)
% THINK: Are your estimated bkg values at all uncertain? 
figure()


% EVALUATE signal = data - background
% THINK: What should you do if the background estimate is > data?
%        How could you find the uncertainty in data - background?




% PLOT the signal-only peak with ERROR BARS
figure()


% Great work! Save the data+background and signal-only plots as ANALYSIS
%             Save a DAY2 workspace!

%% DAY 3 = CHARACTERIZATION of your signal 
% LOAD your DAY2 workspace

% EXTRACT the characteristics (mean, width, uncerts) of your signal peak
% THINK: Which statistical distribution describes your signal peak?
% TOOL: Curve-fitting app has a GUI to walk you through fitting, with 
%       automatic plotting and quality prints! Select x and y data and the 
%       function you'd like. (use Fit -> Save to Workspace) 
% SAVE a screencapture of your fit and its parameters.


% COMPARE: NSignal in signal peak to NBackground under the peak region
% THINK: how can you find the number of events in the signal peak?
%        how can you find the number of bkg events under the peak?
% PRINT: these values along with their uncertainties



% Almost done!
% THINK: Can you statistically distinguish signal from background?
%        Can you find this particle with a web search for you mass?
%        Research this particle (pdg.gov), find its width (capital Gamma)
%        Do your mass & width agree with the known values? Find percent
%        differences and also discrepancy/significance.
%        If your width is *much* larger than accepted, why might this be?