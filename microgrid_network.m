%% Define parameters

% Network topology 
M = 5; % Number of microgrids
L = ones(M, M) - eye(M, M); % Links between microgrids

% General microgrid parameters 
longitudes = []; % Microgrid longitudes 
latitudes =  []; % Microgrid latitudes 
B_char = []; % Charging efficiency
B_dis = [];  % Discharging efficiency 

% Wind modeling parameters (these are vectors because they differ depending on
% the microgrid)
vc = 3 * ones(1,M); % Cut-in speed (3-4 m/s)
vr = 12 * ones(1,M); % Rated speed (12-13 m/s)
vf = 25 * ones(1,M); % Cut-out speed (25 m/s)
Pr = [500 500 500 500 500]; % Rated power (kW) (100 kW - 1 MW for microgrids)
k = [2.537 3.463 2.971 2.963 2.385];  % Weibull shape parameter - download from https://globalwindatlas.info/en/
c = [8.29 9.01 9.38 10.39 9.19];  % Weibull scale parameter - download from https://globalwindatlas.info/en/

% Solar modeling parameters 
Spv = [];   % Solar cell area (m^2)
Pf = [];    % Packing factor
epv = [];   % Module reference efficiency
epc = [];   % Power conditioning efficiency 
phi = [];   % Parameter of Beta distribution
theta = []; % Parameter of Beta distribution

% Control parameters
Nc = 24; % Control horizon (hours)
Np = 24; % Prediction horizon (hours)
dt = 1;  % Time step (hours)
psi_dno = []; % Interest of purchasing power from DNO
psi_mg = [];  % Interest of purchasing power from MG
phi_dno = []; % Interest of selling power to DNO
phi_mg = [];  % Interest of selling power to MG

%% Plot the Weibull PDF 

figure(1)
for i=1:M 
    X = [0:0.2:40];
    Y = wblpdf(X, c(i), k(i));
    plot(X,Y)
    xlabel("Windspeed (m/s)")
    ylabel("Probability")
    hold on
end
hold off