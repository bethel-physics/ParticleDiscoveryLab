%% Particle Discovery Lab
% SETUP
% Download the MATLAB workspace for this exercise. Go to the Home menu, Environment 
% section and select "Set Path" -- choose the location of your workspace file. 
% Then run:

load('DoubleMuParked_100k.mat')
%% DAY 1 = RECONSTRUCTION
% Use energy and momentum conservation to reconstruct a decay of particle X 
% -> mu + mu

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
%% HISTOGRAMMING -- create mass and KE histograms
% THINK: What do you expect your kinetic energy histogram to look like? What 
% do you expect your mass histogram to look like? 
% 
% Vocab: imagine plot with 3 bins on x-axis: 0-10, 10-20, 20-30 
% "Bin edges": 0, 10, 20, 30 
% "Bin low edges": 0, 10, 20
% "Bin centers": 5, 15, 25 (want dots on plot to be here!) 
% "Bin width": 10
% 
% Draw a |histogram| of counts versus mass. The |histogram| function needs: 
% Mass values list, Mass "bin edges" list 

figure(1)



hold on % this lets you draw on top of this plot later


 
% THINK: What should the ERROR BARS be for each bin? What should you do if the 
% bin has ZERO entries?
% 
% Tools:  
% * |histcounts|: gives counts, needs bin edges + input values (mass)
% * |errorbar|: draws dots+bars, needs bin centers, y values, down uncert list, 
% up uncert list 




hold off;


%% Draw another |histogram| of counts vs kinetic energy. Add error bars.
figure(2)



hold on; 



hold off;


% Great work! 
% 
% * *SAVE* these plots to represent your raw data in your report.
% * *SAVE* a DAY1 workspace to start from next time (Click in the workspace 
% area, then CTRL-S)
%
%
%% DAY 2 = FITTING -- fit background on either side of the peak
% LOAD your DAY1 workspace
% 
% Vocab: imagine a mass plot with a bump in the middle 
% * "Peak window": region along x-axis under the peak
% * "background": smoothly falling slope of random events, including some of 
% the events in the peak window
% * "signal": events in the peak window minus the background
% 
% CHOOSE mass values or bin numbers for where the peak lies. 
peakmin = input('Type in your peak min boundary (in GeV): ');
peakmax = input('Type in your peak max boundary (in GeV): ');

% Convert these GeV boundaries into bin numbers:



 
% REMOVE the peak window completely from your list of: mass bin edges, mass 
% bin centers, mass counts, mass uncertainties. This forms your BACKGROUND dataset






 
% PERFORM a polynominal fit to the background
% 
% THINK: Which type of curve do you expect will match your data best? Imagine 
% a curve connecting the two sides under your peak.
% 
% Tool: |pollsf| gives fit params, uncerts, y-values, chi^2 value. Needs bin 
% centers, counts, uncertainties, N params

numpars = input('How many params in your polynominal? 1 (flat), 2 (line), etc: ');




 
% EVALUATE your fit by chi^2 and plotting 
% * Plotting: does the shape make any sense? 
% * Chi^2 is defined in Eq. 29. It describes the difference between the points 
% and the fitted curve. LARGER chi^2 tends to mean more difference or scatter 
% of points.
% * OPTIMALLY, Chi^2 / (# points - # parameters) is around 1
% 
% REPEAT fitting until you are satisfied with both of these metrics





figure(3)



hold on;



hold off;


%% SUBTRACTION -- now you will subtract that background from data
% THINK: How will you estimate background in the signal peak window?
% 
% What do you expect the curve to look like after background subtraction?






 
% PLOT the background curve on top of your mass histogram (save it!)
% 
% THINK: Are your estimated bkg values at all uncertain?

figure(4)



hold on;



hold off;

%% 
% EVALUATE signal = data - background
% 
% THINK: What should you do if the background estimate is > data? How could 
% you find the uncertainty in data - background?





 
% PLOT the signal-only peak with ERROR BARS

figure(5)




% Great work! 
% 
% * *SAVE* the data+background and signal-only plots as analysis in your report.
% * *SAVE* a DAY2 workspace!


%% DAY 3 = CHARACTERIZATION of your signal
% LOAD your DAY2 workspace
% 
% EXTRACT the characteristics of your signal peak
% THINK: Which statistical distribution describes your signal peak?
% 
% TOOL nlinfit: fits a given nonlinear function to data 
%   * outputs: [paramsGaus, R, J, covariance]. We need 
%       * paramsGaus (the fitted function parameters)
%       * covariance (uncert^2 for each param runs down the diagonal)
%   * inputs: nlinfit(x-axis data, y-axis data, function handle, 
%     initial parameter guess)
%   * usage example: 
%     [paramsGaus, R, J, covariance] = nlinfit(X, Y, @functionName, guess)
% 
% LOOK at your signal peak to create a vector of initial guesses for the
% Gaussian parameters. THINK: what do the parameters represent?
% 
% CALL "nlinfit" on your signal peak using the function "gausfit".  


 

% EXTRACT the MEAN and WIDTH of your fitted curve. Print these values
% with their uncertainties.
 


%% 
% PLOT the Gaussian fit on top of your signal peak as a final figure. 
% Include a LEGEND. 

gaus_x = Min:(Max-Min)/100:Max; % a fine-division x-axis for a smooth curve
gaus_y = gausfit(paramsGaus,gaus_x); % evaluate your fit at all x-axis values
 
figure(6)
% Repeat your plotting instructions from figure(5) first




%% Almost done!
% THINK: 
% * Can you statistically distinguish signal from background?
% * Can you find this particle with a web search for you mass?
% * Research this particle (pdg.lbl.gov), find its width (capital Gamma)
% * Do your mass & width agree with the known values within their uncertainties? 
%   Find percent differences and also discrepancy/significance.
% * If your width is *much* larger than accepted, why might this be?
%
%
%
% Great work! 
% Save plots and all the numerical values (mean, width, Nsignal, Nbackground, 
% all with uncertainties, percent differences, discrepancy significances) for 
% the Results of your report.