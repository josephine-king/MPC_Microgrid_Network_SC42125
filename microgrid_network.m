% Define parameters
% Control parameters
Nc = 24; % Control horizon (hours)
Np = 24; % Prediction horizon (hours)
dt = 1;  % Time step (hours)
n = 3;   % Number of state variables
m = 24;  % Number of input variables
dt = 1;  % Time step in hours
% Microgrid network topology 
M = 3; % Number of microgrids
L = ones(M, M) - eye(M, M); % Links between microgrids
% General microgrid parameters 
beta_c = 1; % Charging efficiency
beta_d = 1;  % Discharging efficiency

% Initialize microgrids
mg1 = define_microgrid(51.93, 4.5, 0, 3, 12, 25, 5000, 0.8, 0.2, 0.8, 1000, beta_c, beta_d); % Solar only 
mg2 = define_microgrid(51.93, 4.5, 800, 3, 12, 25, 0, 0.8, 0.2, 0.8, 1000, beta_c, beta_d);  % Wind only 
mg3 = define_microgrid(51.93, 4.5, 400, 3, 12, 25, 2500, 0.8, 0.2, 0.8, 1000, beta_c, beta_d); % Hybrid
mgs = [mg1, mg2, mg3];

% Initialize MPC config
mpc_config = get_mpc_config(n, m, Nc, Np, dt);
demand = load_demand(2016);

%% Solve using our own optimization loop with CVX
num_time_steps = 72;
num_hours = num_time_steps + Np;

% Load energy and demand data
[wt, pv] = get_energy_data(mgs, 5, 1, num_hours);
mg1_demand = get_demand_data(demand.demand_data, 5, 15, 1, num_hours, "DE_KN_industrial1_grid_import");
mg2_demand = 20.*get_demand_data(demand.demand_data, 5, 15, 1, num_hours, "DE_KN_industrial2_grid_import");
mg3_demand = 20.*get_demand_data(demand.demand_data, 5, 15, 1, num_hours, "DE_KN_residential1_grid_import");
D = [mg1_demand; mg2_demand; mg3_demand];

% Initialize state and inputs
x = zeros(n, num_time_steps+1);
u = zeros(m, num_time_steps+1);
for i = 1:n
    x(i,1) = mgs(i).cap./2;
end
u(:,1) = zeros(m,1);

for k = 1:num_time_steps
    % Get energy generation and demand data
    wt_k = wt(:, k:k+Np);
    pv_k = pv(:, k:k+Np);
    D_k = D(:, k:k+Np);
    % Get optimal input from MPC controller
    u(:,k+1) = solve_mpc(x(:,k), mpc_config, mgs, wt_k, pv_k, D_k);
    % Apply state update
    x(:,k+1) = reshape(state_function(x(:,k), u(:,k+1), Nc, M, beta_c, beta_d, wt_k, pv_k, D_k), n, 1);
end

function first_u = solve_mpc(state, mpc_config, mgs, wt, pv, D)

    N = mpc_config.Np;
    m = mpc_config.m;
    n = mpc_config.n;

    charge_idx = 1;
    discharge_idx = 2;
    grid_sell_idx = 3;
    mg_sell_idx_start = 4;
    mg_sell_idx_end = 5;
    grid_buy_idx = 6;
    mg_buy_idx_start = 7;
    mg_buy_idx_end = 8;
    num_mg_inputs = 8;

    cap = [];
    for mg = mgs
        cap = [cap; mg.cap];
    end

    cvx_begin 
        variable u(m,N); % Control sequence
        variable x(n,N); % State sequence
        
        % Define cost function
        J = 0;
        for i = 1:N
            for mg = 1:n
                mg_offset = num_mg_inputs*(mg-1);
                J = J + u(grid_buy_idx + mg_offset, i) + u(grid_sell_idx + mg_offset, i); 
                %for offset = mg_sell_idx_start:mg_sell_idx_end
                %    J = J + u(offset + mg_offset, i);
                %end
            end
        end
        minimize(J)

        % Constraints
        subject to

            % Initial state
            x(:,1) == state;  
            % State update
            for i = 1:N-1
                for mg = 1:n
                    mg_offset = num_mg_inputs*(mg-1);
                    beta_c = mgs(mg).beta_c;
                    beta_d = mgs(mg).beta_d;
                    % State update equation
                    x(mg,i+1) == x(mg,i) + beta_c*u(charge_idx + mg_offset, i) - beta_d*u(discharge_idx + mg_offset, i);
                end
            end

            for i = 1:N
                % Energy balance constraints 
                u_i = u(:,i);
                %u_bal_lhs = zeros(n, 1);
                %u_bal_rhs = zeros(n, 1);
                for mg = 1:n
                    % MG variables
                    mg_offset = num_mg_inputs*(mg-1);
                    beta_c = mgs(mg).beta_c;
                    beta_d = mgs(mg).beta_d;
                    % Left side of equation 
                    ubal_lhs(mg) = wt(mg,i) + pv(mg,i) - D(mg,i);
                    % Right side of equation
                    energy_sold = sum(u_i(grid_sell_idx + mg_offset : mg_sell_idx_end + mg_offset));
                    energy_bought = sum(u_i(grid_buy_idx + mg_offset : mg_buy_idx_end + mg_offset));
                    u_char = u_i(charge_idx + mg_offset);
                    u_dischar = u_i(discharge_idx + mg_offset);
                    ubal_rhs(mg) = energy_sold - energy_bought + beta_c*u_char - beta_d*u_dischar;
                end
                ubal_lhs - ubal_rhs == 0;
               
                % Energy sold from one grid must match energy bought from other grid
                u_i(4)  - u_i(15) == 0;
                u_i(5)  - u_i(23) == 0;
                u_i(13) - u_i(24) == 0;
                u_i(12) - u_i(7)  == 0;
                u_i(20) - u_i(8)  == 0;
                u_i(21) - u_i(16) == 0;

                % Constraints on states and inputs
                0.2.*cap <= x(:,i) <= 0.8.*cap;
                0 <= u_i;
            end

    cvx_end

    % Get first control input
    first_u = u(:,1);

end

%% Plotting
figure(1)
plot(x(1,:))
hold on
plot(x(2,:))
hold on
plot(x(3,:))
legend(["MG1", "MG2", "MG3"])
title("Energy stored")
ylabel("Energy stored (kWh)")
xlabel("Time step (hour)")
%ylim([0,max(cap)+50])

figure(2)
plot(u(4,:))
hold on
plot(u(5,:))
hold on
plot(u(12,:))
hold on
plot(u(13,:))
hold on
plot(u(20,:))
hold on
plot(u(21,:))
hold on
legend(["MG1 to MG2", "MG1 to MG3", "MG2 to MG1", "MG2 to MG3", "MG3 to MG1", "MG3 to MG2"])
title("Power sold to MGs")
ylabel("Power sold (kW)")
xlabel("Time step (hour)")

figure(3)
plot(u(7,:))
hold on
plot(u(8,:))
hold on
plot(u(15,:))
hold on
plot(u(16,:))
hold on
plot(u(23,:))
hold on
plot(u(24,:))
hold on
legend(["MG1 to MG2", "MG1 to MG3", "MG2 to MG1", "MG2 to MG3", "MG3 to MG1", "MG3 to MG2"])
title("Power purchased from MGs")
ylabel("Power purchased (kW)")
xlabel("Time step (hour)")

figure(4)
plot(u(3,:))
hold on
plot(u(11,:))
hold on
plot(u(19,:))
legend(["MG1", "MG2", "MG3"])
title("Power sold to DNO")
ylabel("Power sold (kW)")
xlabel("Time step (hour)")

figure(5)
plot(u(6,:))
hold on
plot(u(14,:))
hold on
plot(u(22,:))
legend(["MG1", "MG2", "MG3"])
title("Power purchased from DNO")
ylabel("Power purchased (kW)")
xlabel("Time step (hour)")

figure(6)
plot(D(1,:))
hold on
plot(D(2,:))
hold on
plot(D(3,:))
legend(["MG1", "MG2", "MG3"])
title("Power demand")
ylabel("Power demand (kW)")
xlabel("Time step (hour)")

%% Testing and plotting energy modeling code

month = 5;
expected_wind_power = zeros(M,24,1);
expected_solar_power = zeros(M,24,1);

mgs = [mg1, mg2, mg3];
num_hours = 72;
[expected_wind_power, expected_solar_power] = get_energy_data(mgs, month, 1, num_hours);

function [wt, pv] = get_energy_data(mgs, month, start_hour, num_hours)
    
    wt = zeros(length(mgs), num_hours);
    pv = zeros(length(mgs), num_hours);
    for m = 1:length(mgs)
        mg = mgs(m);

        for hour = start_hour:start_hour+num_hours-1

            hour_in_day = mod(hour, 24) + 1;
    
            [beta_params, weibull_params] = load_params([mg.latitude],[mg.longitude]);
            weibull_params_mg = get_weibull_params(weibull_params, [mg.latitude], [mg.longitude], month, hour_in_day);
            beta_params_mg = get_beta_params(beta_params, [mg.latitude], [mg.longitude], month, hour_in_day);
            
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
            
            wt(m, hour) = trapz(v, wp .* wbl_pdf);
            pv(m, hour) = trapz(i, sp .* beta_pdf);
        
        end
    end
end

figure(1)
plot(0:num_hours-1, expected_solar_power, "o-")
xlabel("Hour")
ylabel("Expected solar power (kW)")
title("Expected hourly solar power generation for one day in month ", month)
legend(["MG1", "MG2", "MG3", "MG4", "MG5"])

figure(2)
plot(0:num_hours-1, cumsum(expected_solar_power,2), "o-")
xlabel("Hour")
ylabel("Expected solar power (kW)")
title("Expected cumulative solar power for one day in month ", month)
legend(["MG1", "MG2", "MG3", "MG4", "MG5"])

figure(3)
plot(0:num_hours-1, expected_wind_power, "o-")
xlabel("Hour")
ylabel("Expected wind power (kW)")
title("Expected hourly wind power generation for one day in month ", month)
legend(["MG1", "MG2", "MG3", "MG4", "MG5"])

figure(4)
plot(0:num_hours-1, cumsum(expected_wind_power,2), "o-")
xlabel("Hour")
ylabel("Expected wind power (kW)")
title("Expected cumulative wind power for one day in month ", month)
legend(["MG1", "MG2", "MG3", "MG4", "MG5"])

figure(5)
plot(0:num_hours-1, cumsum(expected_wind_power,2)+cumsum(expected_solar_power,2), "o-")
xlabel("Hour")
ylabel("Expected wind power (kW)")
title("Expected cumulative power generation one day in month ", month)
legend(["MG1", "MG2", "MG3", "MG4", "MG5"])



%% Functions
% Create a microgrid struct with the given parameters
function mg = define_microgrid(latitude, longitude, Pr, vc, vr, vf, Spv, Pf, epv, epc, cap, beta_c, beta_d)
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
    mg.cap = cap;             % Storage capacity, in kWh
    mg.beta_c = beta_c;
    mg.beta_d = beta_d;
    % Optimization parameters
    %mg.psi_dno = psi_dno;     % Interest of purchasing power from DNO
    %mg.psi_mg = psi_mg;       % Interest of purchasing power from MG
    %mg.phi_dno = phi_dno;     % Interest of selling power to DNO
    %mg.phi_mg = phi_mg;       % Interest of selling power to MG
end

function mpc_config = get_mpc_config(n, m, Nc, Np, dt)
    mpc_config.n = n;
    mpc_config.m = m;
    mpc_config.Nc = Nc;
    mpc_config.Np = Np;
    mpc_config.dt = dt;
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

function demand = load_demand(selected_year)
    demand = load(sprintf("energy_demand_%d.mat", selected_year));
end

function demand_data = get_demand_data(demand_table, month, day, hour, num_hours, param_name)
    demand_data = [];
    % Find the row index that matches the given year, month, day, and hour
    idx = find(demand_table.("month") == month & demand_table.("day") == day & demand_table.("hour") == hour, 1);
    % Check if index was found
    if isempty(idx)
        error('No matching timestamp found in the data.');
    end
    % Extract the requested data for the specified number of hours
    end_idx = min(idx + num_hours - 1, height(demand_table)); % Ensure it doesn't exceed table size
    demand_data = table2array(demand_table(idx:end_idx, param_name))'; % Return the subset of data
end


%% Use MPC Controller
% Define MPC Controller
controller = nlmpc(3, 3, 24);
controller.PredictionHorizon = Np;
controller.ControlHorizon = Nc;
controller.Model.StateFcn = 'state_function';
controller.Model.OutputFcn = 'output_function';
controller.Model.IsContinuousTime = false;
controller.Model.NumberOfParameters = 8;
controller.Optimization.CustomCostFcn = 'cost_function';
controller.Optimization.CustomEqConFcn = 'eq_constraints';
controller.Optimization.CustomIneqConFcn = 'ineq_constraints';

num_time_steps = 24;
options = nlmpcmoveopt;

cap = [mg1.cap; mg2.cap; mg3.cap];

% Load energy and demand data
num_hours = num_time_steps + Np;
[wt, pv] = get_energy_data([mg1, mg2, mg3], 5, 1, num_hours);
D = ones(1,num_hours).*(wt(:,1) + pv(:,1) + 30);

% Initial state
x = zeros(n, num_time_steps+1);
u = zeros(m, num_time_steps+1);
x(:,1) = cap/2;
u(:,1) = zeros(m,1);
%u(2) = 30;
%u(10) = 30;
%u(18) = 30;

for k = 1:num_time_steps
    wt_k = wt(:, k:k+Np);
    pv_k = pv(:, k:k+Np);
    D_k = D(:, k:k+Np);
    options.Parameters = {Nc, M, beta_c, beta_d, cap, wt_k, pv_k, D_k};
    u(:,k+1) = nlmpcmove(controller, x(:,k), u(:,k), [], [], options);
    x(:,k+1) = reshape(state_function(x(:,k), u(:,k+1), Nc, M, beta_c, beta_d, wt_k, pv_k, D_k), n, 1);
end