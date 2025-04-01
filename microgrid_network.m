% Define parameters
cvx_solver Mosek
% Control parameters
N = 24; % Control and prediction horizon (hours)
dt = 1;  % Time step (hours)
n = 3;   % Number of state variables (number of microgrids)
m = 21;  % Number of input variables
dt = 1;  % Time step in hours
% Microgrid network topology 
M = 3;
L = ones(n, n) - eye(n, n); % Links between microgrids
% General microgrid parameters 
beta_c = 0.9; % Charging efficiency
beta_d = 0.9;  % Discharging efficiency

% Initialize microgrids
% Solar only, residential
mg1 = define_microgrid(1, L, 51.93, 4.5, 0, 3, 12, 25, 5000, 0.8, 0.2, 0.8, 2000, beta_c, beta_d); % Solar only 
% Wind only, industrial
mg2 = define_microgrid(2, L, 51.93, 4.5, 800, 3, 12, 25, 0, 0.8, 0.2, 0.8, 2000, beta_c, beta_d);  % Wind only 
% Wind + solar - "public"
mg3 = define_microgrid(3, L, 51.93, 4.5, 400, 3, 12, 25, 2500, 0.8, 0.2, 0.8, 2000, beta_c, beta_d); % Hybrid
mgs = [mg1, mg2, mg3];

% Initialize MPC config
mpc_config = get_mpc_config(n, m, N, dt);
demand = load_demand(2016);

%% Solve using our own optimization loop with CVX
num_time_steps = 72;
num_hours = num_time_steps + N;
sel_month = 6;

% Load energy and demand data
[wt, pv] = get_energy_data(mgs, sel_month, 1, num_hours);
mg1_demand = 80.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_residential1_grid_import");
mg2_demand = 4.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_industrial1_grid_import");
mg3_demand = 10.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_public1_grid_import");
D_true = [mg1_demand; mg2_demand; mg3_demand];
D_pred = 20 + D_true;

% Initialize state and inputs
x = zeros(n, num_time_steps+1);
u = zeros(m, num_time_steps+1);
for i = 1:n
    x(i,1) = mgs(i).min;
end

for k = 1:num_time_steps
    % Get optimal input from MPC controller
    sprintf("Solving MPC time step %0d", k)
    u(:,k) = solve_mpc(x(:,k), mpc_config, mgs, wt(:,k:k+N), pv(:,k:k+N), D_pred(:,k:k+N));
    % Apply state update
    x(:,k+1) = state_function(x(:,k), u(:,k), wt(:,k), pv(:,k), D_true(:,k), mgs);
end

function first_u = solve_mpc(state, mpc_config, mgs, wt, pv, D)

    N = mpc_config.N;
    m = mpc_config.m;
    n = mpc_config.n;

    charge_idx = 1;
    grid_sell_idx = 2;
    mg_sell_idx_start = 3;
    mg_sell_idx_end = 4;
    grid_buy_idx = 5;
    mg_buy_idx_start = 6;
    mg_buy_idx_end = 7;
    num_mg_inputs = 7;

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
        for k = 1:N
            mg_offset = 0;
            for i = 1:n
                mg = mgs(i);
                J = J + 4*(x(i,k) - mg.ref)^2;
                J = J + 0.5*u(mg.grid_buy_idx + mg_offset, k)^2 + 0.25*u(mg.grid_sell_idx + mg_offset, k)^2;
                mg_offset = mg_offset + mg.num_mg_inputs;
            end
        end
        minimize(J)

        % Constraints
        subject to 

            mg_offset = 0;
            % State constraints     
            for i = 1:n
                mg = mgs(i);

                % Initial state
                x(i,1) == state(i);

                % State update constraints
                for k = 1:N-1
                    x(i,k+1) == x(i,k) + mg.beta_c*u(mg.charge_idx + mg_offset, k);
                end

                for k = 1:N
                    u_k = u(:,k);

                    % State constraints
                    mg.min <= x(i,k) <= mg.max;

                    % Input constraints
                    % Constraints on charging and discharging
                    mg.min_charge <= u_k(mg.charge_idx + mg_offset) <= mg.max_charge;
                    % Constraints on selling and buying 
                    0 <= u_k(mg.grid_sell_idx + mg_offset) <= mg.max_grid_sell;
                    0 <= u_k(mg.grid_buy_idx + mg_offset) <= mg.max_grid_buy;
                    for j = 1:mg.num_connections
                        0 <= u_k(mg.mg_sell_idx_start + mg_offset + j - 1) <= mg.max_mg_sell;
                        0 <= u_k(mg.mg_buy_idx_start + mg_offset + j - 1) <= mg.max_mg_buy;
                    end

                    % Energy balance constraints 
                    % Left side of equation 
                    ubal_lhs(i,k) = wt(i,k) + pv(i,k) - D(i,k);
                    % Right side of equation
                    energy_sold = sum(u_k(mg.grid_sell_idx + mg_offset : mg.mg_sell_idx_end + mg_offset));
                    energy_bought = sum(u_k(mg.grid_buy_idx + mg_offset : mg.mg_buy_idx_end + mg_offset));
                    ubal_rhs(i,k) = energy_sold - energy_bought + mg.beta_c*u_k(mg.charge_idx + mg_offset);
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
    % Get first control input
    first_u = u(:,1);

end

%% Plotting
figure(1)
plot(100.*x(1,:)/mg1.cap)
hold on
plot(100.*x(2,:)/mg2.cap)
hold on
plot(100.*x(3,:)/mg3.cap)
legend(["MG1", "MG2", "MG3"])
title("Battery Percent Charged (%)")
ylabel("Battery Percent Charged (%)")
xlabel("Time step (hour)")

figure(2)
plot(u(3,:))
hold on
plot(u(4,:))
hold on
plot(u(10,:))
hold on
plot(u(11,:))
hold on
plot(u(17,:))
hold on
plot(u(18,:))
hold on
legend(["MG1 to MG2", "MG1 to MG3", "MG2 to MG1", "MG2 to MG3", "MG3 to MG1", "MG3 to MG2"])
title("Power sold to MGs")
ylabel("Power sold (kW)")
xlabel("Time step (hour)")

figure(3)
plot(u(6,:))
hold on
plot(u(7,:))
hold on
plot(u(13,:))
hold on
plot(u(14,:))
hold on
plot(u(20,:))
hold on
plot(u(21,:))
hold on
legend(["MG1 to MG2", "MG1 to MG3", "MG2 to MG1", "MG2 to MG3", "MG3 to MG1", "MG3 to MG2"])
title("Power purchased from MGs")
ylabel("Power purchased (kW)")
xlabel("Time step (hour)")

figure(4)
plot(u(2,:))
hold on
plot(u(9,:))
hold on
plot(u(16,:))
legend(["MG1", "MG2", "MG3"])
title("Power sold to DNO")
ylabel("Power sold (kW)")
xlabel("Time step (hour)")

figure(5)
plot(u(5,:))
hold on
plot(u(12,:))
hold on
plot(u(19,:))
legend(["MG1", "MG2", "MG3"])
title("Power purchased from DNO")
ylabel("Power purchased (kW)")
xlabel("Time step (hour)")

figure(6)
plot(D_true(1,:))
hold on
plot(D_true(2,:))
hold on
plot(D_true(3,:))
legend(["MG1", "MG2", "MG3"])
title("Power demand")
ylabel("Power demand (kW)")
xlabel("Time step (hour)")

figure(7)
plot(wt(1,:) + pv(1,:) - D_true(1,:))
hold on
plot(wt(2,:) + pv(2,:) - D_true(2,:))
hold on
plot(wt(3,:) + pv(3,:) - D_true(3,:))
legend(["MG1", "MG2", "MG3"])
title("Power balance")
ylabel("Power balance (kW)")
xlabel("Time step (hour)")

figure(8)
plot(u(2,:)+u(3,:)+u(4,:))
hold on
plot(u(5,:)+u(6,:)+u(7,:))
title("MG1 total power bought and sold")
legend(["Sold", "Bought"])


%% Functions

% State function
function next_x = state_function(x, u, wt, pv, D, mgs)
    n = length(x);
    mg_offset = 0;
    for mg_idx = 1:n
        mg = mgs(mg_idx);
        energy_sold = sum(u(mg.grid_sell_idx + mg_offset : mg.mg_sell_idx_end + mg_offset));
        energy_bought = sum(u(mg.grid_buy_idx + mg_offset : mg.mg_buy_idx_end + mg_offset));
        next_x(mg_idx) = x(mg_idx) + energy_bought - energy_sold + wt(mg_idx) + pv(mg_idx) - D(mg_idx);
        mg_offset = mg_offset + mgs(mg_idx).num_mg_inputs;
    end
    reshape(next_x, [n,1]);
end

% Create a microgrid struct with the given parameters
function mg = define_microgrid(index, adj_matrix, latitude, longitude, Pr, vc, vr, vf, Spv, Pf, epv, epc, cap, beta_c, beta_d)
    
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
    mg.epv = epv;             % Module reference efficiency (10 - 23% depending on the type of panel)
    mg.epc = epc;             % Power conditioning efficiency (need to figure this out, just put 1 for now)
    mg.cap = cap;             % Storage capacity, in kWh
    mg.beta_c = beta_c;       % Charging efficiency
    mg.beta_d = beta_d;       % Discharging efficiency

    % Constraints
    % Need to find sources for these values
    mg.max = 0.8*cap;
    mg.min = 0.2*cap;
    mg.ref = 0.6*cap;
    mg.max_grid_buy = 500;
    mg.max_grid_sell = 500;
    mg.max_mg_buy = 300;
    mg.max_mg_sell = 300;
    mg.max_charge = 1000;
    mg.min_charge = -1000;
    mg.max_power_bal = mg.num_connections * mg.max_mg_sell + mg.max_grid_sell + mg.max_charge;
    mg.min_power_bal = -(mg.num_connections * mg.max_mg_buy + mg.max_grid_buy) + mg.min_charge;

    % Input indices 
    mg.charge_idx = 1;
    mg.grid_sell_idx = 2;
    mg.mg_sell_idx_start = 3;
    mg.mg_sell_idx_end = mg.mg_sell_idx_start + mg.num_connections - 1;
    mg.grid_buy_idx = mg.mg_sell_idx_end + 1;
    mg.mg_buy_idx_start = mg.grid_buy_idx + 1;
    mg.mg_buy_idx_end = mg.mg_buy_idx_start + mg.num_connections - 1;
    mg.num_mg_inputs = mg.mg_buy_idx_end;

end

function mpc_config = get_mpc_config(n, m, N, dt)
    mpc_config.n = n;
    mpc_config.m = m;
    mpc_config.N = N;
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

