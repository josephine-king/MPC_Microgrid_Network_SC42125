%% Initialize MPC + Microgrid parameters

% Control parameters
N = 3; % Control and prediction horizon (hours)
dt = 1;  % Time step (hours)
n = 3;   % Number of state variables (number of microgrids)
m = 21;  % Number of input variables
dt = 1;  % Time step in hours
nd = 3;  % Number of disturbances
% Microgrid network topology 
L = ones(n, n) - eye(n, n); % Links between microgrids

% General microgrid parameters 
rte = sqrt(0.87); % Round trip efficiency, square root because we only look at one direction (charge or discharge)
% Choose cost function
cost_function = 1;
if (cost_function == 1)
    state_weight = 16;
    buy_weight = 2;
    sell_weight = 1;
elseif (cost_function == 2)
    state_weight = 1;
    buy_weight = 4;
    sell_weight = 4;
else
    state_weight = 0;
    buy_weight = 1;
    sell_weight = 1;
end
max_grid_buy = 500;
max_grid_sell = 500;
max_mg_buy = 500;
max_mg_sell = 500;
max_charge = 1000;
min_charge = -1000;

% Initialize microgrids
% Solar only, residential
mg1 = define_microgrid(1, L, 51.93, 4.5, 0, 3, 12, 25, 6000, 0.7, 0.2, 2000, rte, ...
    max_grid_buy, max_grid_sell, max_mg_buy, max_mg_sell, max_charge, min_charge, ...
    state_weight, buy_weight, sell_weight);
% Wind only, industrial
mg2 = define_microgrid(2, L, 51.93, 4.5, 800, 3, 12, 25, 0, 0.7, 0.2, 2000, rte, ...
    max_grid_buy, max_grid_sell, max_mg_buy, max_mg_sell, max_charge, min_charge, ...
    state_weight, buy_weight, sell_weight);   
% Wind + solar - "public"
mg3 = define_microgrid(3, L, 51.93, 4.5, 400, 3, 12, 25, 3000, 0.7, 0.2, 2000, rte, ...
    max_grid_buy, max_grid_sell, max_mg_buy, max_mg_sell, max_charge, min_charge, ...
    state_weight, buy_weight, sell_weight); 
mgs = [mg1, mg2, mg3];

% Initialize MPC config
[A,B,C,Q,R] = get_system_matrices(n, m, mgs);
[P,K,L] = idare(A,B,Q,R);
mpc_config = get_mpc_config(n, m, N, dt, A, B, C, Q, R, P);

% Load the demand data
demand = load_demand(2016);

%% Initialize simulation parameters

num_time_steps = 72;
num_hours = num_time_steps + N;
sel_month = 6;
use_real_data = true;

% Load energy and demand data
if (use_real_data == true)
    [wt, pv] = get_energy_data(mgs, sel_month, 1, num_hours);
    mg1_demand = 40.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_residential1_grid_import");
    mg1_demand = mg1_demand + 40.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_residential2_grid_import");
    mg2_demand = 3.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_industrial1_grid_import");
    mg3_demand = 10.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_public1_grid_import");
    D = [mg1_demand; mg2_demand; mg3_demand];
else
    wt = zeros(n,num_hours);
    pv = zeros(n,num_hours);
    D = zeros(n,num_hours);
end

% Initialize state and inputs
x = zeros(n, num_time_steps+1);
u = zeros(m, num_time_steps+1);
x(:,1) = [mgs.min]';

% Initialize observer
estimate_disturbance = 0;
if (estimate_disturbance)
    [A_aug, B_aug, C_aug, L] = init_observer(A, B, C, eye(n), zeros(n,n), [.8, .8, .8, .6, .6, .6]'); 
    xhat_aug = zeros(n + nd, num_time_steps);
    xhat_aug(1:n,1) = x(:,1);
end

% Initialize disturbance
disturbance_type = "no_disturbance";
if (disturbance_type == "const_disturbance")
    d1 = get_const_disturbance(40, num_time_steps);
    d2 = get_const_disturbance(40, num_time_steps);
    d3 = get_const_disturbance(40, num_time_steps);
    d = [d1; d2; d3];
elseif (disturbance_type == "ramp_disturbance")
    d1 = get_ramp_disturbance(-20, 20, num_time_steps);
    d2 = get_ramp_disturbance(-20, 20, num_time_steps);
    d3 = get_ramp_disturbance(-20, 20, num_time_steps);
    d = [d1; d2; d3];
else
    d = zeros(n, num_time_steps);
    estimate_disturbance = 0;
end 

% Simulate system with MPC
for k = 1:num_time_steps
    disp(["Solving MPC time step ", num2str(k)]);
       
    power_bal = wt(:,k:k+N-1) + pv(:,k:k+N-1) - D(:,k:k+N-1);

    if (estimate_disturbance)
        dhat = xhat_aug(n+1:end,k);
        xhat = xhat_aug(1:n,k);
        [xref, uref, power_bal_ref] = ots(mgs, mpc_config, [mgs.ref]', dhat);
        power_bal(:, 1) = power_bal(:, 1) - power_bal_ref;
    else
        xref = [mgs.ref]';
        uref = zeros(m,1);
        xhat = x(:,k);
    end

    % Solve optimal control problem
    u(:,k) = solve_mpc(mpc_config, mgs, xhat-xref, xref, uref, power_bal);
    % State update
    x(:,k+1) = A * x(:,k) + B * u(:,k) + d(:,k);
    
    if (estimate_disturbance)
        xhat_aug(:,k+1) = estimate(x(:,k+1), xhat_aug(:,k), u(:,k), A_aug, B_aug, C_aug, L);
    end
end

% Solves the MPC optimal control problem 
function first_u = solve_mpc(cfg, mgs, state, xref, uref, power_bal)

    N = cfg.N;
    m = cfg.m;
    n = cfg.n;
    A = cfg.A;
    B = cfg.B;
    Q = cfg.Q;
    R = cfg.R;

    cvx_begin quiet
        variable u(m,N); % Control sequence
        variable x(n,N); % State sequence

        % Binary variables for hybrid constraints
        % delta = 1 when we sell, 0 when we buy
        variable delta(n,N) binary

        expression ubal_lhs(n, N)
        expression ubal_rhs(n, N)

        % Define cost function
        J = 0;
        % Stage cost
        for k = 1:N
            J = J + x(:,k)'*Q*x(:,k);
            % Have to do this in a loop, cvx won't accept such a large
            % matrix
            for u_idx = 1:m
                J = J + u(u_idx, k)^2 * R(u_idx, u_idx);
            end
        end
        minimize(J)

        % Constraints
        subject to 

            % Initial state
            x(:,1) == A * state + B * u(:,1);

            % State update constraints
            for k = 1:N-1
                x(:,k+1) == A * x(:,k) + B * u(:,k);
            end

            % Terminal constraints 
            x(:,N) == zeros(n,1);

            mg_offset = 0;
            % State constraints     
            for i = 1:n
                mg = mgs(i);

                for k = 1:N
                    u_k = u(:,k);

                    % State constraints
                    mg.min - xref(i) <= x(i,k) <= mg.max - xref(i);

                    % Input constraints
                    mg.u_min <= u_k(mg_offset+1 : mg_offset+mg.num_mg_inputs) <= mg.u_max;

                    % Energy balance constraints 
                    % Left side of equation 
                    ubal_lhs(i,k) = power_bal(i,k);
                    % Right side of equation
                    energy_sold = sum(u_k(mg.grid_sell_idx + mg_offset : mg.mg_sell_idx_end + mg_offset));
                    energy_bought = sum(u_k(mg.grid_buy_idx + mg_offset : mg.mg_buy_idx_end + mg_offset));
                    ubal_rhs(i,k) = energy_sold - energy_bought + mg.rte*u_k(mg.charge_idx + mg_offset);
                    % Left and right side must be equal
                    ubal_lhs(i,k) - ubal_rhs(i,k) == 0;
                
                    % Energy sold from one grid must match energy bought from other grid
                    % If two MGs are not connected, they should not
                    % exchange power
                    % Loop through the microgrid's connections
                    connection_offset = 0;
                    for j = 1:n
                        % Skip the current microgrid
                        if (j == i) 
                            connection_offset = connection_offset + mgs(j).num_mg_inputs;
                            continue
                        end
                        if (j < i)
                            offseti = -2;
                            offsetj = -1;
                        else
                            offseti = -1;
                            offsetj = -2;
                        end
                        % Not connected - power bought/sold is constrained
                        % to 0
                        if isempty(find(mg.connections == j))
                            u_k(mg.mg_sell_idx_start + j + offsetj + mg_offset) == 0;
                            u_k(mgs(j).mg_buy_idx_start + i + offseti + connection_offset) == 0;
                        % Connected - power sold from one grid must equal
                        % power bought by another
                        else
                            u_k(mg.mg_sell_idx_start + j + offsetj + mg_offset) == u_k(mgs(j).mg_buy_idx_start + i + offseti + connection_offset);
                        end
                        connection_offset = connection_offset + mgs(j).num_mg_inputs;
                    end

                    % We cannot buy and sell at the same time
                    % If we are selling, u_bal > 0
                    u_k(mg.grid_sell_idx + mg_offset) <= mg.max_grid_sell * delta(i,k);                    
                    u_k(mg.mg_sell_idx_start + mg_offset) <= mg.max_mg_sell * delta(i,k);    
                    u_k(mg.mg_sell_idx_end + mg_offset) <= mg.max_mg_sell * delta(i,k);
                    ubal_lhs(i,k) >= mg.min_power_bal * (1 - delta(i,k));
                    % If we are buying, u_bal < 0
                    u_k(mg.grid_buy_idx + mg_offset) <= mg.max_grid_buy * (1 - delta(i,k));
                    u_k(mg.mg_buy_idx_start + mg_offset) <= mg.max_mg_buy * (1 - delta(i,k));    
                    u_k(mg.mg_buy_idx_end + mg_offset) <= mg.max_mg_buy * (1 - delta(i,k));
                    ubal_lhs(i,k) <= mg.max_power_bal * delta(i,k);

                end 
                
                mg_offset = mg_offset + mg.num_mg_inputs;

            end

    cvx_end

    sprintf("Optimal val: %0d", J)
    first_u = u(:,1) + uref;

end

%% Plotting


close all

figure(1)
stairs(100.*x(1,:)/mg1.cap)
hold on
stairs(100.*x(2,:)/mg2.cap)
hold on
stairs(100.*x(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
legend(["MG1", "MG2", "MG3"])
title("Battery Percent Charged")
ylim([20 70])
ylabel("Battery Percent Charged (%)")
xlabel("Time step (hour)")
fontsize(16,"points")
grid on

figure(4)
subplot(3,1,1)
stairs(100.*x(1,:)/mg1.cap)
hold on
stairs(100.*x(2,:)/mg2.cap)
hold on
stairs(100.*x(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged")
ylim([20 70])
ylabel("Battery Percent Charged (%)")
fontsize(14,"points")
grid on
lgd1 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd1,10,'points')

subplot(3,1,2)
stairs(u(2,:))
hold on
stairs(u(9,:))
hold on
stairs(u(16,:))
hold on
stairs(u(5,:))
hold on
stairs(u(12,:))
hold on
stairs(u(19,:))
title("Power exchanged with DNO")
ylabel("Power (kW)")
fontsize(14,"points")
grid on
ylim([0 400])
lgd2 = legend(["MG1 to DNO", "MG2 to DNO", "MG3 to DNO", "DNO to MG1", "DNO to MG2", "DNO to MG3"]);

subplot(3,1,3)
stairs(u(3,:))
hold on
stairs(u(4,:))
hold on
stairs(u(10,:))
hold on
stairs(u(11,:))
hold on
stairs(u(17,:))
hold on
stairs(u(18,:))
hold on
title("Power sold to MGs")
ylabel("Power sold (kW)")
fontsize(14,"points")
grid on
ylim([0 400])
lgd3 = legend(["MG1 to MG2", "MG1 to MG3", "MG2 to MG1", "MG2 to MG3", "MG3 to MG1", "MG3 to MG2"]);

fontsize(lgd1,12,'points')
fontsize(lgd2,12,'points')
fontsize(lgd3,12,'points')

figure(6)
stairs(wt(1,:) + pv(1,:) - D(1,:))
hold on
stairs(wt(2,:) + pv(2,:) - D(2,:))
hold on
stairs(wt(3,:) + pv(3,:) - D(3,:), 'Color',[0.4,0.7,0.3])
legend(["MG1", "MG2", "MG3"])
title("Power balance")
ylabel("Power balance (kW)")
xlabel("Time step (hour)")
fontsize(16,"points")

figure(7)
subplot(3,1,1)
stairs(100.*x(1,:)/mg1.cap)
hold on
stairs(100.*x(2,:)/mg2.cap)
hold on
stairs(100.*x(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, No Observer or OTS")
ylim([20 70])
ylabel("Battery Percent Charged (%)")
fontsize(14,"points")
grid on
lgd1 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd1,10,'points')

subplot(3,1,2)
stairs(100.*x_ots(1,:)/mg1.cap)
hold on
stairs(100.*x_ots(2,:)/mg2.cap)
hold on
stairs(100.*x_ots(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, Observer and OTS")
ylim([20 70])
ylabel("Battery Percent Charged (%)")
fontsize(14,"points")
grid on
lgd1 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd1,10,'points')

subplot(3,1,3)
stairs(xhat_aug_ots(4,:))
hold on
stairs(xhat_aug_ots(5,:))
hold on
stairs(xhat_aug_ots(6,:))
hold on
plot(d(1,:), "Color", [0.5, 0.5, 0.5])
legend(["MG1 dhat", "MG2 dhat", "MG3 dhat", "True disturbance value"])
title("Disturbance Estimate")
ylabel("Disturbance Estimate (kWh)")
xlabel("Time step (hour)")
ylim([-25,25])
fontsize(16,"points")

%% Initialization functions

% Compiles the MPC config parameters into a struct
function mpc_config = get_mpc_config(n, m, N, dt, A, B, C, Q, R, P)
    mpc_config.n = n;
    mpc_config.m = m;
    mpc_config.N = N;
    mpc_config.dt = dt;
    mpc_config.A = A;
    mpc_config.B = B;
    mpc_config.C = C;
    mpc_config.Q = Q;
    mpc_config.R = R;
    mpc_config.P = P;
end

% Create a microgrid struct with the given parameters
function mg = define_microgrid(index, adj_matrix, latitude, longitude, Pr, vc, vr, vf, Spv, Pf, eff, cap, rte, ...
    max_grid_buy, max_grid_sell, max_mg_buy, max_mg_sell, max_charge, min_charge, state_weight, buy_weight, sell_weight)
    
    % Connections
    mg.num_connections = sum(adj_matrix(index,:));
    mg.connections = find(adj_matrix(index,:));

    % Microgrid parameters
    mg.latitude = latitude;   % Latitude of the MG
    mg.longitude = longitude; % Longitude of the MG
    mg.Pr = Pr;               % Wind power rating (kW) (100 kW - 1 MW for microgrids)
    mg.vc = vc;               % Cut-in speed (3-4 m/s)
    mg.vr = vr;               % Rated speed (12-13 m/s)
    mg.vf = vf;               % Cut-out speed (25 m/s)
    mg.Spv = Spv;             % Solar cell area (m^2) (500 - 800 m^2 for microgrids)
    mg.Pf = Pf;               % Packing factor (30 - 50%)
    mg.eff = eff;             % Module efficiency (10 - 23% depending on the type of panel)
    mg.cap = cap;             % Storage capacity, in kWh
    mg.rte = rte;             % Round trip efficiency

    % Constraints
    mg.max = 0.8*cap;
    mg.min = 0.2*cap;
    mg.ref = 0.6*cap;
    mg.max_grid_buy = max_grid_buy;
    mg.max_grid_sell = max_grid_sell;
    mg.max_mg_buy = max_mg_buy;
    mg.max_mg_sell = max_mg_sell;
    mg.max_charge = max_charge; 
    mg.min_charge = min_charge;
    mg.max_power_bal = mg.num_connections * mg.max_mg_sell + mg.max_grid_sell + mg.max_charge;
    mg.min_power_bal = -(mg.num_connections * mg.max_mg_buy + mg.max_grid_buy) + mg.min_charge;
    mg.u_max = [mg.max_charge; mg.max_grid_sell; mg.max_mg_sell; mg.max_mg_sell; mg.max_grid_buy; mg.max_mg_buy; mg.max_mg_buy];
    mg.u_min = [mg.min_charge; 0; 0; 0; 0; 0; 0];

    % Input indices 
    mg.charge_idx = 1;
    mg.grid_sell_idx = 2;
    mg.mg_sell_idx_start = 3;
    mg.mg_sell_idx_end = mg.mg_sell_idx_start + mg.num_connections - 1;
    mg.grid_buy_idx = mg.mg_sell_idx_end + 1;
    mg.mg_buy_idx_start = mg.grid_buy_idx + 1;
    mg.mg_buy_idx_end = mg.mg_buy_idx_start + mg.num_connections - 1;
    mg.num_mg_inputs = mg.mg_buy_idx_end;

    % Cost function parameters
    mg.state_weight = state_weight;
    mg.grid_buy_weight = buy_weight;
    mg.grid_sell_weight = sell_weight;

end

% Gets the system matrices 
function [A,B,C,Q,R] = get_system_matrices(n, m, mgs)
    A = eye(n);
    B = zeros(n,m);
    Q = zeros(n,n);
    R = zeros(m,m);
    C = eye(n);
    mg_offset = 0;

    for i = 1:n
        B(i, mgs(i).charge_idx + mg_offset) = mgs(i).rte;
        Q(i,i) = mgs(i).state_weight;
        R(mg_offset + mgs(i).grid_sell_idx, mg_offset + mgs(i).grid_sell_idx) = mgs(i).grid_sell_weight;
        R(mg_offset + mgs(i).grid_buy_idx, mg_offset + mgs(i).grid_buy_idx) = mgs(i).grid_buy_weight;
        mg_offset = mg_offset + mgs(i).num_mg_inputs;
    end
    % PSD for solvers
    R = R + 1e-6*eye(m,m);

end

%% Observer and OTS functions

% Calculates the optimal targets xref and uref, and the offset to the power
% balance, power_bal_ref
function [xref, uref, power_bal_ref] = ots(mgs, mpc_config, yref, dhat)
    n = mpc_config.n;
    m = mpc_config.m;
    R = mpc_config.R;
    xref = yref;
    uref = zeros(m,1);
    mg_offset = 0;
    for i = 1:n
        mg = mgs(i);
        uref(mg_offset + mg.charge_idx) = -dhat(i)/mg.rte;
        power_bal_ref = mg.rte*uref(mg_offset + mg.charge_idx);
        mg_offset = mg_offset + mg.num_mg_inputs;
    end
end

% Gets a constant disturbance 
function disturbance = get_const_disturbance(val, num_hours)
    t = ones(1, num_hours);
    disturbance = t.*val;
end

function disturbance = get_ramp_disturbance(start, stop, num_hours)
    t = ones(1, num_hours);
    disturbance = linspace(start, stop, num_hours);
end

% Initialize Luenberger observers for all MGs
function [A_aug, B_aug, C_aug, L] = init_observer(A, B, C, Bd, Cd, poles) 
    n = size(A, 1);
    m = size(B, 2);
    A_aug = [A, Bd; zeros(n,n), C];
    B_aug = [B; zeros(n,m)];
    C_aug = [C, Cd];
    L = place(A_aug', C_aug', poles)';
end

% Estimates the state using the Luenberger observer 
function xhat_new = estimate(y, xhat, u, A_aug, B_aug, C_aug, L)
    xhat_new = A_aug * xhat + B_aug * u + L * (y - C_aug * xhat);
end

%% Power generation and demand modeling functions

% Loads the beta solar parameters and weibull wind parameters
function [beta_params, weibull_params] = load_params(latitudes, longitudes)
    beta_params = containers.Map(); 
    weibull_params = containers.Map(); 

    for i = 1:length(latitudes)
        coordinates = mat2str([latitudes(i), longitudes(i)]); 
        beta_params(coordinates) = load(latitudes(i) + "_" + longitudes(i) + "/beta_params.mat");
        weibull_params(coordinates) = load(latitudes(i) + "_" + longitudes(i) + "/weibull_params.mat");
    end

end

% Gets the weibull parameters for a given month and hour
function params = get_weibull_params(weibull_params, latitude, longitude, month, hour)
    coordinates = mat2str([latitude, longitude]); 
    data = weibull_params(coordinates);
    params = [data.weibull_params_a{month, hour}, data.weibull_params_b{month, hour}];
end

% Gets the beta parameters for a given month and hour
function params = get_beta_params(beta_params, latitude, longitude, month, hour)
    coordinates = mat2str([latitude, longitude]); 
    data = beta_params(coordinates);
    params = [data.beta_params_a{month, hour}, data.beta_params_b{month, hour}, data.beta_params_scaling{month, hour}];
end

% Loads the power demand for a given year
function demand = load_demand(selected_year)
    demand = load(sprintf("energy_demand_%d.mat", selected_year));
end

% Gets the demand data for a given month, dat, hour, and customer for num_hours
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

% Uses the weibull and beta parameters to calculate the probabilistic power
% generation
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
            sp = mg.Pf*mg.Spv*mg.eff.*i;
            
            wt(m, hour) = trapz(v, wp .* wbl_pdf);
            pv(m, hour) = trapz(i, sp .* beta_pdf);
        
        end
    end
end
