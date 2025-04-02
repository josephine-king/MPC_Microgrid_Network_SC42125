% Define parameters
%cvx_solver Mosek
% Control parameters
N = 24; % Control and prediction horizon (hours)
dt = 1;  % Time step (hours)
n = 3;   % Number of state variables (number of microgrids)
m = 21;  % Number of input variables
dt = 1;  % Time step in hours
% Microgrid network topology 
L = ones(n, n) - eye(n, n); % Links between microgrids
% General microgrid parameters 
rte = 0.9; % Round trip efficiency 

% Initialize microgrids
% Solar only, residential
mg1 = define_microgrid(1, L, 51.93, 4.5, 0, 3, 12, 25, 5000, 0.8, 0.2, 0.8, 2000, rte); % Solar only 
% Wind only, industrial
mg2 = define_microgrid(2, L, 51.93, 4.5, 800, 3, 12, 25, 0, 0.8, 0.2, 0.8, 2000, rte);  % Wind only 
% Wind + solar - "public"
mg3 = define_microgrid(3, L, 51.93, 4.5, 400, 3, 12, 25, 2500, 0.8, 0.2, 0.8, 2000, rte); % Hybrid
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
mg1_demand = 80.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_residential1_grid_import") + 10;
mg2_demand = 4.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_industrial1_grid_import") + 10;
mg3_demand = 10.*get_demand_data(demand.demand_data, sel_month, 15, 1, num_hours, "DE_KN_public1_grid_import") + 10;
D_pred = [mg1_demand; mg2_demand; mg3_demand];

% amp = 10;
% freq = 2 * pi / 24;
% mg1_disturbance = get_sine_disturbance(amp, freq, 0, num_hours);
% mg2_disturbance = get_sine_disturbance(-amp, freq, 0, num_hours);
% mg3_disturbance = get_sine_disturbance(1.5*amp, freq, 0, num_hours);
% D_true = [mg1_demand + mg1_disturbance; mg2_demand + mg2_disturbance; mg3_demand + mg3_disturbance];
D_true = D_pred - 20;

% Initialize observer
xhat_aug = cell(n,1);
xhat_log = cell(n,1);
[A_base, B_base, C_base, L] = init_observer(1, 2*pi/24, [0.7, 0.6, 0.5]);
A_aug_all = cell(n,1);
B_aug_all = cell(n,1);
C_aug_all = cell(n,1);
L_all = cell(n,1);
dhat = cell(n,1);
for i = 1:n
    A_aug_all{i} = A_base;
    B_aug_all{i} = B_base;
    C_aug_all{i} = C_base;
    L_all{i} = L;
    xhat_aug{i} = [mgs(i).min; 0; 0];
    xhat_log{i} = zeros(3, num_time_steps+1);
    xhat_log{i}(:,1) = xhat_aug{i};
end

% Initialize state and inputs
x = zeros(n, num_time_steps+1);
u = zeros(m, num_time_steps+1);
for i = 1:n
    x(i,1) = mgs(i).min;
end

A = zeros(n,1);
freq = zeros(n,1);
for k = 1:num_time_steps
    disp(["Solving MPC time step ", num2str(k)]);
    
    u(:,k) = solve_mpc(mpc_config, mgs, x(:,k), [mgs.ref], wt(:,k:k+N), pv(:,k:k+N), D_pred(:,k:k+N), 20);
    %u(:,k) = solve_mpc(x(:,k), mpc_config, mgs, wt(:,k:k+N), pv(:,k:k+N), D_pred(:,k:k+N), est_disturbance);
    x(:,k+1) = state_function(x(:,k), u(:,k), wt(:,k), pv(:,k), D_true(:,k), mgs);

    for i = 1:n
        y_k = x(i,k);
        A_aug = A_aug_all{i};
        B_aug = B_aug_all{i};
        C_aug = C_aug_all{i};
        L = L_all{i};

        u_idx = (i-1)*(m/n) + 1;  % Assuming charge index is first in block
        xhat_aug{i} = A_aug * xhat_aug{i} + B_aug * u(u_idx,k) + L * (y_k - C_aug * xhat_aug{i});
        xhat_log{i}(:,k+1) = xhat_aug{i};

    end
end


function terminal_set = estimate_terminal_set()

    rand_ics = rand(1,3);
    ic = [mgs.min] + ([mgs.max] - [mgs.min]).*rand_ics;


end

function first_u = solve_mpc(mpc_config, mgs, state, xref, wt, pv, D, dhat)

    N = mpc_config.N;
    m = mpc_config.m;
    n = mpc_config.n;

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
                J = J + 4*(x(i,k) - xref(i))^2;
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
                x(i,1) == state(i) + mg.rte*u(mg.charge_idx + mg_offset, 1) + dhat;

                % State update constraints
                for k = 1:N-1
                    x(i,k+1) == x(i,k) + mg.rte*u(mg.charge_idx + mg_offset, k) + dhat;
                end

                for k = 1:N
                    u_k = u(:,k);

                    % State constraints
                    mg.min <= x(i,k) <= mg.max;

                    % Input constraints
                    mg.u_min <= u_k(mg_offset+1 : mg_offset+mg.num_mg_inputs) <= mg.u_max;

                    % Energy balance constraints 
                    % Left side of equation 
                    ubal_lhs(i,k) = wt(i,k) + pv(i,k) - D(i,k);
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
    first_u = u(:,1);

end

%% Plotting

close all

figure(1)
plot(100.*x(1,:)/mg1.cap,"--")
hold on
plot(100.*x(2,:)/mg2.cap,"--")
hold on
plot(100.*x(3,:)/mg3.cap,"--")
legend(["MG1", "MG2", "MG3"])
title("Battery Percent Charged (%)")
%ylim([55 65])
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
hold on 
plot(D_pred(1,:), '--')
hold on
plot(D_pred(2,:), '--')
hold on
plot(D_pred(3,:), '--')
legend(["MG1", "MG2", "MG3", "MG1 Pred", "MG2 Pred", "MG3 Pred"])
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
function mg = define_microgrid(index, adj_matrix, latitude, longitude, Pr, vc, vr, vf, Spv, Pf, epv, epc, cap, rte)
    
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
    mg.rte = rte;             % Round trip efficiency

    % Constraints
    % Need to find sources for these values
    mg.max = 0.8*cap;
    mg.min = 0.2*cap;
    mg.ref = 0.6*cap;
    mg.max_grid_buy = 500;
    mg.max_grid_sell = 500;
    mg.max_mg_buy = 300;
    mg.max_mg_sell = 300;
    mg.max_charge = 600; % Tesla megapack
    mg.min_charge = -600;
    mg.max_power_bal = mg.num_connections * mg.max_mg_sell + mg.max_grid_sell + mg.max_charge;
    mg.min_power_bal = -(mg.num_connections * mg.max_mg_buy + mg.max_grid_buy) + mg.min_charge;

    mg.u_max = [mg.max_charge; mg.max_grid_sell; mg.max_mg_sell; mg.max_mg_sell; mg.max_grid_buy; mg.max_mg_buy; mg.max_mg_buy]
    mg.u_min = [mg.min_charge; 0; 0; 0; 0; 0; 0;]

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

function disturbance = get_sine_disturbance(amp, freq, phase, num_hours)
    t = 0:num_hours-1;
    disturbance = amp * sin(freq * t + phase);
end

function disturbance = get_const_disturbance(val, num_hours)
    t = ones(1, num_hours);
    disturbance = t.*val;
end

% Initialize Luenberger observers for all MGs
function [A_base, B_base, C_base, L] = init_observer(dt, omega, poles) 
    A_base = [1, 1, 0; 0, 1, dt*omega; 0, -dt*omega, 1];
    B_base = [1; 0; 0];
    C_base = [1, 0, 0];
    L = place(A_base', C_base', poles)';
end