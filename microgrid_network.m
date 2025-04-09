% Define parameters
%cvx_solver Mosek
% Control parameters
N = 24; % Control and prediction horizon (hours)
dt = 1;  % Time step (hours)
n = 3;   % Number of state variables (number of microgrids)
m = 21;  % Number of input variables
dt = 1;  % Time step in hours
nd = 3;  % Number of disturbances
% Microgrid network topology 
L = ones(n, n) - eye(n, n); % Links between microgrids
% General microgrid parameters 
rte = sqrt(0.87); % Round trip efficiency, square root because we only look at one direction (charge or discharge)

% Initialize microgrids
% Solar only, residential
mg1 = define_microgrid(1, L, 51.93, 4.5, 0, 3, 12, 25, 6000, 0.7, 0.2, 2000, rte); % Solar only 
% Wind only, industrial
mg2 = define_microgrid(2, L, 51.93, 4.5, 800, 3, 12, 25, 0, 0.7, 0.2, 2000, rte);  % Wind only 
% Wind + solar - "public"
mg3 = define_microgrid(3, L, 51.93, 4.5, 400, 3, 12, 25, 3000, 0.7, 0.2, 2000, rte); % Hybrid
mgs = [mg1, mg2, mg3];

% Initialize MPC config
[A,B,C,Q,R] = get_system_matrices(n, m, mgs);
[P,K,L] = idare(A,B,Q,R);
mpc_config = get_mpc_config(n, m, N, dt, A, B, C, Q, R, P);
%Ak = A + B*K;

demand = load_demand(2016);


%% Solve using our own optimization loop with CVX
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
    mg1_demand = get_const_disturbance(10, num_hours);
    mg2_demand = get_const_disturbance(20, num_hours);
    mg3_demand = get_const_disturbance(30, num_hours);
    wt = [30; 20; 10] .* ones(3, num_hours);
    pv = [20; 10; 20] .* ones(3, num_hours);
    D = [mg1_demand; mg2_demand; mg3_demand];
end

% Initialize state and inputs
x = zeros(n, num_time_steps+1);
u = zeros(m, num_time_steps+1);
for i = 1:n
    x(i,1) = mgs(i).max/2;
end

% Initialize observer
[A_aug, B_aug, C_aug, L] = init_observer(A, B, C, eye(n), zeros(n,n), [.95, .95, .95, .7, .7, .7]'); 
xhat_aug = zeros(n + nd,1);
xhat_aug(1:n,:) = x(:,1);

% Initialize disturbance
d = [0; 0; 0];

for k = 1:num_time_steps
    disp(["Solving MPC time step ", num2str(k)]);
        
    dhat = xhat_aug(n+1:end,:);
    xhat = xhat_aug(1:n,:);
    dhat = zeros(n,1);
    [xref, uref] = ots(mgs, mpc_config, [mgs.ref]', d);
    u(:,k) = solve_mpc(mpc_config, mgs, xhat-xref, xref, uref, wt(:,k:k+N), pv(:,k:k+N), D(:,k:k+N));
    
    x(:,k+1) = A * x(:,k) + B * u(:,k) + d;
    
    xhat_aug = estimate(x(:,k+1), xhat_aug, u(:,k), A_aug, B_aug, C_aug, L);
end

function [xref, uref] = ots(mgs, mpc_config, yref, dhat)
    n = mpc_config.n;
    m = mpc_config.m;
    % The equations simplify to the following:
    xref = yref;
    % -B * u = Bd * dhat
    % This just gives us -rte * uch = dhat, so uch is -dhat/rte
    % All the other inputs are ideally 0
    uref = zeros(m, 1);
    mg_offset = 0;
    for i = 1:n
        mg = mgs(i);
        uref(mg_offset + mg.charge_idx) = -dhat(i)/mg.rte;
        mg_offset = mg_offset + mg.num_mg_inputs;
    end
end

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

function first_u = solve_mpc(cfg, mgs, state, xref, uref, wt, pv, D)

    N = cfg.N;
    m = cfg.m;
    n = cfg.n;
    A = cfg.A;
    B = cfg.B;
    C = cfg.C;
    Q = cfg.Q;
    R = cfg.R;
    P = cfg.P;

    cvx_begin quiet
        variable u(m,N); % Control sequence
        variable x(n,N); % State sequence

        % Binary variables for hybrid constraints
        % delta = 1 when we sell, 0 when we buy
        variable delta(n,N) binary

        expression ubal_lhs(n, N)
        expression ubal_rhs(n, N)

        % Define cost function
        % Terminal cost
        %J = 0.5*x(:,N)'*P*x(:,N);
        J = 0;
        % Stage cost
        for k = 1:N
            J = J + x(:,k)'*Q*x(:,k);
            % Have to do this in a loop because cvx will not accept non-PD
            % R matrix (ours is PSD but not PD)
            mg_offset = 0;
            for i = 1:n
                mg = mgs(i);
                J = J + u(mg_offset + mg.grid_buy_idx, k)^2 * R(mg_offset + mg.grid_buy_idx, mg_offset + mg.grid_buy_idx);
                J = J + u(mg_offset + mg.grid_sell_idx, k)^2 * R(mg_offset + mg.grid_sell_idx, mg_offset + mg.grid_sell_idx);
                mg_offset = mg_offset + mg.num_mg_inputs;
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
            %x(:,N) == zeros(n,1);

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
lgd1 = legend(["MG1", "MG2", "MG3"])
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
title("Power balance, October 15-17")
ylabel("Power balance (kW)")
xlabel("Time step (hour)")
fontsize(16,"points")


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
function mg = define_microgrid(index, adj_matrix, latitude, longitude, Pr, vc, vr, vf, Spv, Pf, eff, cap, rte)
    
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
    % Need to find sources for these values
    mg.max = 0.8*cap;
    mg.min = 0.2*cap;
    mg.ref = 0.6*cap;
    mg.max_grid_buy = 500;
    mg.max_grid_sell = 500;
    mg.max_mg_buy = 500;
    mg.max_mg_sell = 500;
    mg.max_charge = 1000; % Tesla megapack
    mg.min_charge = -1000;
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
    mg.state_weight = 16;
    mg.grid_buy_weight = 2;
    mg.grid_sell_weight = 1;

end

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
            sp = mg.Pf*mg.Spv*mg.eff.*i;
            
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
function [A_aug, B_aug, C_aug, L] = init_observer(A, B, C, Bd, Cd, poles) 
    n = size(A, 1);
    m = size(B, 2);
    A_aug = [A, Bd; zeros(n,n), C];
    B_aug = [B; zeros(n,m)];
    C_aug = [C, Cd];
    L = place(A_aug', C_aug', poles)';
end

function xhat_new = estimate(y, xhat, u, A_aug, B_aug, C_aug, L)
    xhat_new = A_aug * xhat + B_aug * u + L * (y - C_aug * xhat);
end