%% Define parameters

% Control parameters
Nc = 24; % Control horizon (hours)
Np = 24; % Prediction horizon (hours)
dt = 1;  % Time step (hours)

% Microgrid network topology 
M = 5; % Number of microgrids
L = ones(M, M) - eye(M, M); % Links between microgrids

% General microgrid parameters 
B_char = []; % Charging efficiency
B_dis = [];  % Discharging efficiency


%% Some sample code to test it out
month = 5;
expected_wind_power = zeros(M,24,1);
expected_solar_power = zeros(M,24,1);

mg1 = define_microgrid(51.93, 4.5, 500, 3, 12, 25, 500, 0.8, 0.2, 0.8, 1, 1, 1, 1);
mg2 = define_microgrid(51.93, 4.5, 400, 3, 12, 25, 600, 0.8, 0.2, 0.8, 1, 1, 1, 1);
mg3 = define_microgrid(51.93, 4.5, 300, 3, 12, 25, 700, 0.8, 0.2, 0.8, 1, 1, 1, 1);
mg4 = define_microgrid(51.93, 4.5, 200, 3, 12, 25, 800, 0.8, 0.2, 0.8, 1, 1, 1, 1);
mg5 = define_microgrid(51.93, 4.5, 100, 3, 12, 25, 900, 0.8, 0.2, 0.8, 1, 1, 1, 1);

mgs = [mg1, mg2, mg3, mg4, mg5];

for mg_idx = 1:M
    for hour = 1:24
        mg = mgs(mg_idx);

        [beta_params, weibull_params] = load_params([mg.latitude],[mg.longitude]);
        weibull_params_mg = get_weibull_params(weibull_params, [mg.latitude], [mg.longitude], month, hour);
        beta_params_mg = get_beta_params(beta_params, [mg.latitude], [mg.longitude], month, hour);
        
        v = linspace(0, 60, 1000);
        wbl_pdf = wblpdf(v, weibull_params_mg(1), weibull_params_mg(2));
        i = linspace(0.001, beta_params_mg(3), 1000);
    
        if (beta_params_mg(3) == 0)
            beta_pdf = zeros(1, 1000);
        else
            beta_pdf = betapdf(i, beta_params_mg(1), beta_params_mg(2));
            beta_pdf = beta_pdf*beta_params_mg(3);
        end
    
        % Wind power function
        wp = zeros(size(v));
        wp(v >= mg.vc & v < mg.vr) = mg.Pr * (v(v >= mg.vc & v < mg.vr) - mg.vc) / (mg.vr - mg.vc);
        wp(v >= mg.vr & v <= mg.vf) = mg.Pr;
    
        % Solar power function
        sp = mg.Pf*mg.Spv*mg.epc*mg.epv.*i;
        
        expected_wind_power(mg_idx, hour) = trapz(v, wp .* wbl_pdf);
        expected_solar_power(mg_idx, hour) = trapz(i, sp .* beta_pdf);
    
    end
end

figure(1)
plot(0:23, expected_solar_power, "o-")
xlabel("Hour")
ylabel("Expected solar power (kW)")
title("Expected hourly solar power generation for one day in month ", month)
legend(["MG1", "MG2", "MG3", "MG4", "MG5"])

figure(2)
plot(0:23, cumsum(expected_solar_power,2), "o-")
xlabel("Hour")
ylabel("Expected solar power (kW)")
title("Expected cumulative solar power for one day in month ", month)
legend(["MG1", "MG2", "MG3", "MG4", "MG5"])

figure(3)
plot(0:23, expected_wind_power, "o-")
xlabel("Hour")
ylabel("Expected wind power (kW)")
title("Expected hourly wind power generation for one day in month ", month)
legend(["MG1", "MG2", "MG3", "MG4", "MG5"])

figure(4)
plot(0:23, cumsum(expected_wind_power,2), "o-")
xlabel("Hour")
ylabel("Expected wind power (kW)")
title("Expected cumulative wind power for one day in month ", month)
legend(["MG1", "MG2", "MG3", "MG4", "MG5"])


%% Functions
% Create a microgrid struct with the given parameters
function mg = define_microgrid(latitude, longitude, Pr, vc, vr, vf, Spv, Pf, epv, epc, psi_dno, psi_mg, phi_dno, phi_mg)
    % Microgrid parameters
    mg.latitude = latitude;   % Latitude of the MG
    mg.longitude = longitude; % Longitude of the MG
    mg.Pr = Pr;               % Wind power rating (kW) (100 kW - 1 MW for microgrids)
    mg.vc = vc;               % Cut-in speed (3-4 m/s)
    mg.vr = vr;               % Rated speed (12-13 m/s)
    mg.vf = vf;               % Cut-out speed (25 m/s)
    mg.Spv = Spv;             % Solar cell area (m^2) (500 - 800 m^2 for microgrids)
    mg.Pf = Pf;               % Packing factor (30 - 50%)
    mg.epv = epv;             % Module reference efficiency (10 - 23% depending on the type of panel)
    mg.epc = epc;             % Power conditioning efficiency (need to figure this out, just put 1 for now)

    % Optimization parameters
    mg.psi_dno = psi_dno;     % Interest of purchasing power from DNO
    mg.psi_mg = psi_mg;       % Interest of purchasing power from MG
    mg.phi_dno = phi_dno;     % Interest of selling power to DNO
    mg.phi_mg = phi_mg;       % Interest of selling power to MG
end

function [beta_params, weibull_params] = load_params(latitudes, longitudes)
    beta_params = containers.Map(); 
    weibull_params = containers.Map(); 

    for i = 1:length(latitudes)
        coordinates = mat2str([latitudes(i), longitudes(i)]); 
        beta_params(coordinates) = load(latitudes(i) + "_" + longitudes(i) + "/beta_params.mat");
        weibull_params(coordinates) = load(latitudes(i) + "_" + longitudes(i) + "/weibull_params.mat");
    end

end

function params = get_weibull_params(weibull_params, latitude, longitude, month, hour)
    coordinates = mat2str([latitude, longitude]); 
    data = weibull_params(coordinates);
    params = [data.weibull_params_a{month, hour}, data.weibull_params_b{month, hour}];
end

function params = get_beta_params(beta_params, latitude, longitude, month, hour)
    coordinates = mat2str([latitude, longitude]); 
    data = beta_params(coordinates);
    params = [data.beta_params_a{month, hour}, data.beta_params_b{month, hour}, data.beta_params_scaling{month, hour}];
end
