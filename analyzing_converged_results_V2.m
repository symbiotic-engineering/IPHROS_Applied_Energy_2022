% Analyzing converged results, take 2

% Goes with debugging_pareto_IPHRO_RO_for_paper


%% Loading in results from debugging_pareto_... file

data = readmatrix('Opt Results.xlsx','Sheet','gen = 187, new');
x_pareto = data(1:70,1:8);
obj_vals = data(1:70,14:16);
results = [x_pareto,obj_vals];

%% Post Processing

E_er_array = obj_vals(:,1);
V_fwRO_array = obj_vals(:,2);
C_array = obj_vals(:,3);
E_r_array = x_pareto(:,1);
gamma_array = x_pareto(:,2);
gamma_RO_array = x_pareto(:,3);
h_L_array = x_pareto(:,4);
N_e1_array = x_pareto(:,5);
N_e2_array = x_pareto(:,6);
N_pv1_array = x_pareto(:,7);
N_pv2_array = x_pareto(:,8);

% Normalizing
E_er_array_normalized = zeros(70,1);
V_fwRO_array_normalized = zeros(70,1);
C_array_normalized = zeros(70,1);
E_r_array_normalized = zeros(70,1);
gamma_array_normalized = zeros(70,1);
gamma_RO_array_normalized = zeros(70,1);
h_L_array_normalized = zeros(70,1);
N_e1_array_normalized = zeros(70,1);
N_e2_array_normalized = zeros(70,1);
N_pv1_array_normalized = zeros(70,1);
N_pv2_array_normalized = zeros(70,1);
for i = 1:70
    E_er_array_normalized(i) = (E_er_array(i) - min(E_er_array))/(max(E_er_array)-min(E_er_array));
    V_fwRO_array_normalized(i) = (V_fwRO_array(i) - min(V_fwRO_array))/(max(V_fwRO_array)-min(V_fwRO_array));
    C_array_normalized(i) = (C_array(i) - min(C_array))/(max(C_array)-min(C_array));
    E_r_array_normalized(i) = (E_r_array(i) - min(E_r_array))/(max(E_r_array)-min(E_r_array));
    gamma_array_normalized(i) = (gamma_array(i) - min(gamma_array))/(max(gamma_array)-min(gamma_array));
    gamma_RO_array_normalized(i) = (gamma_RO_array(i) - min(gamma_RO_array))/(max(gamma_RO_array)-min(gamma_RO_array));
    h_L_array_normalized(i) = (h_L_array(i) - min(h_L_array))/(max(h_L_array)-min(h_L_array));
    N_e1_array_normalized(i) = (N_e1_array(i) - min(N_e1_array))/(max(N_e1_array)-min(N_e1_array));
    N_e2_array_normalized(i) = (N_e2_array(i) - min(N_e2_array))/(max(N_e2_array)-min(N_e2_array));
    N_pv1_array_normalized(i) = (N_pv1_array(i) - min(N_pv1_array))/(max(N_pv1_array)-min(N_pv1_array));
    N_pv2_array_normalized(i) = (N_pv2_array(i) - min(N_pv2_array))/(max(N_pv2_array)-min(N_pv2_array));
end

% Compiling normalized results
results_norm = [E_r_array_normalized,gamma_array_normalized,...
    gamma_RO_array_normalized,h_L_array_normalized,...
    N_e1_array_normalized,N_e2_array_normalized,N_pv1_array_normalized,...
    N_pv2_array_normalized,E_er_array_normalized,V_fwRO_array_normalized,...
    C_array_normalized];

% Finding the best 10% for each objective
[B_E,I_E] = maxk(E_er_array,7);
[B_V,I_V] = maxk(V_fwRO_array,7);
[B_C,I_C] = mink(C_array,7);

bin_best_E = sort(I_E);
bin_best_V = sort(I_V);
bin_best_C = sort(I_C);
bins = ones(70,1);
for i = 1:7
    bins(bin_best_E(i)) = 2;
    bins(bin_best_V(i)) = 3;
    bins(bin_best_C(i)) = 4;
end

% Rearranging rows for Parallel Coordinate Plot
results_par_coor = [results,bins];
results_par_coor = sortrows(results_par_coor,12);
bins_par_coor = results_par_coor(:,12);
results_par_coor = results_par_coor(:,1:11);


%% Filling out the table
[M_max_E,I_max_E] = max(E_er_array);
[M_max_V,I_max_V] = max(V_fwRO_array);
[M_min_C,I_min_C] = min(C_array);

[M_max_E_norm,I_max_E_norm] = max(E_er_array_normalized);
[M_max_V_norm,I_max_V_norm] = max(V_fwRO_array_normalized);
[M_min_C_norm,I_min_C_norm] = min(C_array_normalized);

utopia_point = [M_max_E_norm,M_max_V_norm,M_min_C_norm];

distances = zeros(70,1);
for i = 1:70
    pareto_point = [E_er_array_normalized(i),V_fwRO_array_normalized(i),C_array_normalized(i)];
    distances(i) = norm(utopia_point - pareto_point);
end

[~,I_best] = min(distances);

%% Running Functions for Constraint Values

% Getting constraint values for all these pareto points
c_collection = zeros(70,11);
%salinity_and_pressure_collection = zeros(70,2);

for i = 1:70
    x_individual = x_pareto(i,:);
    [c,ceq] = constraints(x_individual);
    c_collection(i,:) = c;
    %salinity_and_pressure_collection(i,:) = ceq;
end

%% Constraint Analysis
% Making the constraints just variable values
Sht_constraint = c_collection(:,1)+40;
Qf_min_constraint = (c_collection(:,2)-21620)*(-1);
Qf_max_constraint = c_collection(:,3)+98272;
Qc_min_constraint = (c_collection(:,4)-21620)*(-1);
Qc_max_constraint = c_collection(:,5)+98272;
Qp_min_constraint = (c_collection(:,6)-0)*(-1);
Qp_max_constraint = c_collection(:,7)+8369;
rr_min_constraint = (c_collection(:,8)-0)*(-1);
rr_max_constraint = c_collection(:,9)+0.13;
plane_constraint = c_collection(:,11)+1334.934;

constraint_vals = [Sht_constraint,Qf_min_constraint,Qf_max_constraint,...
    Qc_min_constraint,Qc_max_constraint,Qp_min_constraint,...
    Qp_max_constraint,rr_min_constraint,rr_max_constraint,plane_constraint];


%% Plotting Pareto Front
% % 3D Pareto Front, Scatter
% figure(1)
% scatter3(E_er_array,V_fwRO_array,C_array,[],'k')
% hold on
% scatter3(E_er_array(2),V_fwRO_array(2),C_array(2),[],[0.8 0.4 0],'filled')
% scatter3(E_er_array(1),V_fwRO_array(1),C_array(1),[],[0.9 0.6 0],'filled')
% scatter3(E_er_array(5),V_fwRO_array(5),C_array(5),[],[0.8 0.6 0.7],'filled')
% scatter3(E_er_array(48),V_fwRO_array(48),C_array(48),[],[0.95 0.9 0.25],'filled')
% xlabel('$\dot{E}_{er}$', 'Interpreter','latex','FontSize',14)
% ylabel('$\dot{V}_{fw,RO}$', 'Interpreter','latex','FontSize',14)
% zlabel('$C$', 'Interpreter','latex','FontSize',14)
% 
% % 2D Pareto Front, Scatter, E and V
% figure(2)
% scatter(E_er_array,V_fwRO_array,[],[0 0.45 0.7],'LineWidth',2)
% hold on
% scatter(E_er_array(2),V_fwRO_array(2),[],[0.8 0.4 0],'filled')
% scatter(E_er_array(1),V_fwRO_array(1),[],[0.9 0.6 0],'filled')
% scatter(E_er_array(5),V_fwRO_array(5),[],[0.8 0.6 0.7],'filled')
% scatter(E_er_array(48),V_fwRO_array(48),[],[0.95 0.9 0.25],'filled')
% xlabel('$\dot{E}_{er}$', 'Interpreter','latex','FontSize',18)
% ylabel('$\dot{V}_{fw,RO}$', 'Interpreter','latex','FontSize',18)
% 
% % 2D Pareto Front, Scatter, E and C
% figure(3)
% scatter(E_er_array,C_array,[],[0 0.45 0.7],'LineWidth',2)
% hold on
% scatter(E_er_array(2),C_array(2),[],[0.8 0.4 0],'filled')
% scatter(E_er_array(1),C_array(1),[],[0.9 0.6 0],'filled')
% scatter(E_er_array(5),C_array(5),[],[0.8 0.6 0.7],'filled')
% scatter(E_er_array(48),C_array(48),[],[0.95 0.9 0.25],'filled')
% xlabel('$\dot{E}_{er}$', 'Interpreter','latex','FontSize',18)
% ylabel('$C$', 'Interpreter','latex','FontSize',18)
% 
% % 2D Pareto Front, Scatter, V and C
% figure(4)
% scatter(V_fwRO_array,C_array,[],[0 0.45 0.7],'LineWidth',2)
% hold on
% scatter(V_fwRO_array(2),C_array(2),[],[0.8 0.4 0],'filled')
% scatter(V_fwRO_array(1),C_array(1),[],[0.9 0.6 0],'filled')
% scatter(V_fwRO_array(5),C_array(5),[],[0.8 0.6 0.7],'filled')
% scatter(V_fwRO_array(48),C_array(48),[],[0.95 0.9 0.25],'filled')
% xlabel('$\dot{V}_{fw,RO}$', 'Interpreter','latex','FontSize',18)
% ylabel('$C$', 'Interpreter','latex','FontSize',18)
% 
% % Interpolated surface plot for pareto front
% F = scatteredInterpolant(E_er_array,V_fwRO_array,C_array,'natural','none');
% xgr = linspace(min(E_er_array),max(E_er_array));
% ygr = linspace(min(V_fwRO_array),max(V_fwRO_array));
% [XX,YY] = meshgrid(xgr,ygr);
% ZZ = F(XX,YY);
% 
% figure(5)
% surf(XX,YY,ZZ,'LineStyle','none')
% hold on
% scatter3(E_er_array,V_fwRO_array,C_array,40,'k')
% hold on
% scatter3(E_er_array(2),V_fwRO_array(2),C_array(2),40,[0.8 0.4 0],'filled')
% scatter3(E_er_array(1),V_fwRO_array(1),C_array(1),40,[0.9 0.6 0],'filled')
% scatter3(E_er_array(5),V_fwRO_array(5),C_array(5),40,[0.8 0.6 0.7],'filled')
% scatter3(E_er_array(48),V_fwRO_array(48),C_array(48),40,[0.95 0.9 0.25],'filled')
% xlabel('$\dot{E}_{er}$', 'Interpreter','latex','FontSize',14)
% ylabel('$\dot{V}_{fw,RO}$', 'Interpreter','latex','FontSize',14)
% zlabel('$C$', 'Interpreter','latex','FontSize',14)
% hleg = legend('$Pareto \: Surface$','$Pareto \: Points$','$max(\dot{E}_{er})$','$max(\dot{V}_{fw,RO})$','$min(C)$','$rel. \: max(J)$');
% set(hleg, 'Interpreter', 'latex');
% hleg.FontSize = 14;
% 
% view(30,30)
% 
% Parallel Coordinate Plot
figure(6)
p = parallelplot(results_par_coor);
p.CoordinateTickLabels = ["DVar1";"DVar2";"DVar3";"DVar4";"DVar5";"DVar6";"DVar7";"DVar8";"Obj1";"Obj2";"Obj3"];
p.Jitter = 0;
p.GroupData = bins_par_coor;
p.LegendVisible = 0;


%% Functions
% Objective Function
function J = objective(x)
    % General Parameters
    g = 9.81;             % gravitational acceleration [m/s^2]
    rho_sw = 1023.6;      % generic density of saltwater at 1 atm and 25 degrees celsius [kg/m^3]
    rho_fw = 996.9;       % generic density of freshwater at 1 atm and 25 degrees celsius [kg/m^3]
    eta_hp = 0.894;       % pump-side efficiency
    eta_ht = 0.894;       % pump-side efficiency [-]
    T = 25.0;             % temperature [degrees celsius]
    S_sw = 35.0;          % generic seawater salinity [g salt/kg seawater]
    M_salt = 58.44;       % molar mass of NaCl [g/mol]
    P_p = 14.696;         % permeate pressure, atmospheric for now [psi]
    R = 8.3145;           % universal gas constant [J/(mol*K)]
    MW_salt_avg = 30.548; % average molar mass of NaCl [g/mol]
    Am = 440;             % membrane surface area [ft^2]
    M_cost_e1 = 1.7;      % Proxy for membrane cost of 1st element in RO unit [$/unit]
    M_cost_e2 = 1;        % Proxy for membrane cost of other elements in RO unit [$/unit] 
    
    % Unpack x
    E_r = x(1); gamma = x(2); gamma_RO = x(3); h = x(4);
    N_e1 = x(5); N_e2 = x(6); N_pv1 = x(7); N_pv2 = x(8);
    
    % Go through the modules
    [E_rp,E_rd] = initial_energy_division(E_r,gamma);
    V_wp = pump(E_rp,h,eta_hp,rho_sw,g);
    [V_wRO, V_swht, Q_f] = mountaintop_reservoir(V_wp,gamma_RO,N_pv1);
    [V_fwRO,V_oRO,S_oRO,rho_oRO,eta_ROio,eta_RO,rr_max,rr_min,Q_f_max,Q_f_min,Q_c_max,Q_c_min,Q_p_max,Q_p_min,LHS_SVM] = ro_system(h,Q_f,N_e1,N_e2,N_pv1,N_pv2,rho_sw,rho_fw,T,S_sw,M_salt,P_p,R,MW_salt_avg,Am);

    E_swht = turbine(V_swht,h,rho_sw,g,eta_ht);
    E_htRO = energy_recapture(V_oRO,h,eta_ROio,rho_oRO,g,eta_ht);
    %[rho_ht,S_ht,V_ht] = mixing(V_oRO,V_swht,rho_oRO,S_oRO,gamma_RO,eta_RO,rho_sw,S_sw);
    E_er = energy_summation(E_rd,E_swht,E_htRO);
    C = RO_cost(N_e1,N_e2,N_pv1,N_pv2,M_cost_e1,M_cost_e2);

    E_er_min = -E_er; % maximizing
    V_fwRO_min = -V_fwRO; % maximizing
    C_min = C; % minimizing
    
    % Compiling
    J = [E_er_min,V_fwRO_min,C_min];
end

% Constraint Function
function [c,ceq] = constraints(x)
    % General Parameters
    g = 9.81;             % gravitational acceleration [m/s^2]
    rho_sw = 1023.6;      % generic density of saltwater at 1 atm and 25 degrees celsius [kg/m^3]
    rho_fw = 996.9;       % generic density of freshwater at 1 atm and 25 degrees celsius [kg/m^3]
    eta_hp = 0.894;       % pump-side efficiency
    eta_ht = 0.894;       % pump-side efficiency [-]
    T = 25.0;             % temperature [degrees celsius]
    S_sw = 35.0;          % generic seawater salinity [g salt/kg seawater]
    M_salt = 58.44;       % molar mass of NaCl [g/mol]
    P_p = 14.696;         % permeate pressure, atmospheric for now [psi]
    R = 8.3145;           % universal gas constant [J/(mol*K)]
    MW_salt_avg = 30.548; % average molar mass of NaCl [g/mol]
    Am = 440;             % membrane surface area [ft^2]
    M_cost_e1 = 1.7;      % Proxy for membrane cost of 1st element in RO unit [$/unit]
    M_cost_e2 = 1;        % Proxy for membrane cost of other elements in RO unit [$/unit] 
    
    % Unpack x
    E_r = x(1); gamma = x(2); gamma_RO = x(3); h = x(4);
    N_e1 = x(5); N_e2 = x(6); N_pv1 = x(7); N_pv2 = x(8);
    
    % Go through the modules
    [E_rp,E_rd] = initial_energy_division(E_r,gamma);
    V_wp = pump(E_rp,h,eta_hp,rho_sw,g);
    [V_wRO, V_swht, Q_f] = mountaintop_reservoir(V_wp,gamma_RO,N_pv1);
    [V_fwRO,V_oRO,S_oRO,rho_oRO,eta_ROio,eta_RO,rr_max,rr_min,Q_f_max,Q_f_min,Q_c_max,Q_c_min,Q_p_max,Q_p_min,LHS_SVM] = ro_system(h,Q_f,N_e1,N_e2,N_pv1,N_pv2,rho_sw,rho_fw,T,S_sw,M_salt,P_p,R,MW_salt_avg,Am);

    E_swht = turbine(V_swht,h,rho_sw,g,eta_ht);
    E_htRO = energy_recapture(V_oRO,h,eta_ROio,rho_oRO,g,eta_ht);
    [rho_ht,S_ht,V_ht] = mixing(V_oRO,V_swht,rho_oRO,S_oRO,gamma_RO,eta_RO,rho_sw,S_sw);

    % Salinity constraint
    S_ht_con = S_ht-40;

    % Feed flow rate to pressure vessel constraint
    Q_f_min_con = 21620-Q_f_min;
    Q_f_max_con = Q_f_max-98272;

    % Concentrate flow rate constraint
    Q_c_min_con = 21620-Q_c_min;
    Q_c_max_con = Q_c_max-98272;

    % Permeate flow rate constraint
    Q_p_min_con = -Q_p_min;
    Q_p_max_con = Q_p_max-8369;

    % Recovery rate for element constraints
    rr_min_con = -rr_min;
    rr_max_con = rr_max-0.13;

    % If N_e2 = 0, then N_pv2 = 0
    e2_pv2_con = abs(min([N_e2,1]) - min([N_pv2,1]));

    % Equation for SVM plane
    SVM_con = LHS_SVM;

    % Compiling
    c = [S_ht_con;Q_f_min_con;Q_f_max_con;Q_c_min_con;Q_c_max_con;Q_p_min_con;Q_p_max_con;rr_min_con;rr_max_con;e2_pv2_con;SVM_con];
    ceq = [];
end

% Modules
function [E_rp,E_rd] = initial_energy_division(E_r,gamma)
    E_rp = gamma*E_r; % kWh
    E_rd = (1-gamma)*E_r; % kWh
end

function V_wp = pump(E_rp,h,eta_hp,rho_sw,g)
    E_rp = E_rp*(3.6*10^6); % joules (from kWh)
    V_wp = E_rp*eta_hp/(rho_sw*g*h); % m^3
end

function [V_wRO,V_swht,Q_f] = mountaintop_reservoir(V_wp,gamma_RO,N_pv1)
    V_wRO = gamma_RO*V_wp; % m^3
    V_swht = (1-gamma_RO)*V_wp; % m^3
    Q_f = (V_wRO*264.172)/N_pv1; % gpd (from m^3/day)
end

function [V_fwRO,V_oRO,S_oRO,rho_oRO,eta_ROio,eta_RO,rr_max,rr_min,Q_f_max,Q_f_min,Q_c_max,Q_c_min,Q_p_max,Q_p_min,LHS_SVM] = ro_system(h,Q_f_input,N_e1,N_e2,N_pv1,N_pv2,rho_sw,rho_fw,T,S_sw,M_salt,P_p,R,MW_salt_avg,Am)

    % Initialize oopsies relative to Q_f_input
    if Q_f_input <= 98272 && Q_f_input >= 21620
        oopsies = 0;
        oopsies1 = 0;
        oopsies2 = 0;
    else
        oopsies = 1;
        oopsies1 = 1;
        oopsies2 = 1;
    end

    % Collection Lists
    Q_p_collection = []; % gpd
    P_c_collection = []; % psi
    c_c_collection = []; % mol/L
    Q_c_collection = []; % gpd
    rr_collection = [];  % -
    P_f_collection = []; % psi
    S_f_collection = []; % g/kg
    Q_f_collection = []; % gpd
    
    % Initialize 1st stage element counter
    i1 = 1;
    
    % First stage
    while i1 <= N_e1 && oopsies ~= 1
        % Beginning with feed conditions into an element
        if i1 == 1
            % First element
            P_f = 1.4228424*h*(rho_sw/rho_fw); % psi
            c_f = (S_sw*rho_sw)/(1000*M_salt); % mol/L
            Q_f = Q_f_input;
        else
            % Feed conditions from previous element
            P_f = P_c_collection(i1-1); % psi
            Q_f = Q_c_collection(i1-1); % gpd
            c_f = c_c_collection(i1-1); % mol/L
        end

        pi_f = 1.12*(273+T)*2*c_f;
        S_f = salinity_from_osmotic(pi_f,rho_fw,T,R,MW_salt_avg);

        P_f_collection = [P_f_collection;P_f];
        S_f_collection = [S_f_collection;S_f];
        Q_f_collection = [Q_f_collection;Q_f];

        % Solving for Q_p and other numeric values %%%%%%%%%%%%%%%%%%%%%%%
        [Q_p_num,P_c_num,c_c_num,Q_c_num,rr_num,oopsies1] = symbolic_Q_p_stuff(c_f,T,Q_f,Am,P_f,P_p);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if oopsies1 == 1
            i1 = N_e1 + 1;
            oopsies2 = 1;
        else
            % Storing element outputs for further calculations
            Q_p_collection = [Q_p_collection;Q_p_num]; % gpd
            P_c_collection = [P_c_collection;P_c_num]; % psi
            c_c_collection = [c_c_collection;c_c_num]; % mol/L
            Q_c_collection = [Q_c_collection;Q_c_num]; % gpd
            rr_collection = [rr_collection;rr_num];    % -
            
            % Increasing counter
            i1 = i1 + 1;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up storage variables
    Q_p_collection2 = []; % gpd
    P_c_collection2 = []; % psi
    c_c_collection2 = []; % mol/L
    Q_c_collection2 = []; % gpd
    rr_collection2 = [];  % - 
    P_f_collection2 = []; % psi
    S_f_collection2 = []; % g/kg
    Q_f_collection2 = []; % gpd

    % Second stage
    if N_e2 == 0 || N_pv2 == 0 || oopsies == 1 || oopsies1 == 1
        i2 = N_e2+1;
    else        
        % initialize second stage element counter
        i2 = 1;
    end

    while i2 <= N_e2
        if i2 == 1
            % Feed conditions from previous stage (feed flowrate reconfigured)
            P_f = P_c_collection(end); % psi
            c_f = c_c_collection(end); % mol/L
            
            V_intermediate = Q_c_collection(end)*N_pv1; % gpd
            Q_f = V_intermediate/N_pv2; % gpd
        else
            % Feed conditions from previous element
            P_f = P_c_collection2(end); % psi
            Q_f = Q_c_collection2(end); % gpd
            c_f = c_c_collection2(end); % mol/L
        end

        pi_f = 1.12*(273+T)*2*c_f;
        S_f = salinity_from_osmotic(pi_f,rho_fw,T,R,MW_salt_avg);

        P_f_collection = [P_f_collection;P_f];
        S_f_collection = [S_f_collection;S_f];
        Q_f_collection = [Q_f_collection;Q_f];
        
        % Solving for Q_p and other numeric values %%%%%%%%%%%%%%%%%%%%%%%
        [Q_p_num,P_c_num,c_c_num,Q_c_num,rr_num,oopsies2] = symbolic_Q_p_stuff(c_f,T,Q_f,Am,P_f,P_p);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if oopsies2 == 1
            i2 = N_e2 + 1;
        else
            % Storing element outputs for further calculations
            Q_p_collection2 = [Q_p_collection2;Q_p_num]; % gpd
            P_c_collection2 = [P_c_collection2;P_c_num]; % psi
            c_c_collection2 = [c_c_collection2;c_c_num]; % mol/L
            Q_c_collection2 = [Q_c_collection2;Q_c_num]; % gpd
            rr_collection2 = [rr_collection2;rr_num];    % -
            
            % Increasing counter
            i2 = i2 + 1;
        end
    end
    if oopsies == 1 || oopsies1 == 1 || oopsies2 == 1
        % Dependent vars
        V_fwRO = 0; V_oRO = 0; S_oRO = 1e5; rho_oRO = 1e5; eta_ROio = 0; eta_RO = 0;

        % Constraint vars
        rr_max = 1e6; rr_min = -1e6; Q_f_max = 1e6; Q_f_min = -1e6; 
        Q_c_max = 1e6; Q_c_min = -1e6; Q_p_max = 1e6; Q_p_min = -1e6;
        LHS_SVM = 1e6;
        
    else
        % Processing for RO System Output Values
        V_fwRO = N_pv1*(sum(Q_p_collection)/264.172) + N_pv2*(sum(Q_p_collection2)/264.172); % m^3/day
        
        brine_collective = [Q_c_collection;Q_c_collection2];
        pressure_vessel_collective = [N_pv1;N_pv2];
        V_oRO = (brine_collective(end)/264.172)*pressure_vessel_collective(end);  % m^3/day (from gpd)
        
        c_c_collective = [c_c_collection;c_c_collection2];
        pi_c = 1.12*(273+T)*2*c_c_collective(end); % psi
        S_oRO = salinity_from_osmotic(pi_c,rho_fw,T,R,MW_salt_avg); % g/kg
    
        P_c_collective = [P_c_collection;P_c_collection2];
        rho_oRO = density_from_salinity(S_oRO,P_c_collective(end),rho_fw,T,P_p); % kg/m^3
        
        P_f_1st_element = 1.4228424*h*(rho_sw/rho_fw);                        % psi
        eta_ROio = 1-((P_f_1st_element-P_c_collective(end))/P_f_1st_element); % -
    
        rr_collective = [rr_collection;rr_collection2];
        eta_RO_product_term = 1;
        i_eta_RO = 1;
        N_e = N_e1 + N_e2;
        while i_eta_RO <= N_e
            eta_RO_product_term = eta_RO_product_term*(1-rr_collective(i_eta_RO));
            i_eta_RO = i_eta_RO + 1;
        end
        eta_RO = 1-eta_RO_product_term; % -
    
        %%%%%%%%%%%%%%
        P_f_array = [P_f_collection;P_f_collection2];
        S_f_array = [S_f_collection;S_f_collection2];
        Q_f_array = [Q_f_collection;Q_f_collection2];
        LHS_SVM = -1e6;
        for i = 1:length(P_f_array)
            var = -1334.934 + (18.222*P_f_array(i)) - (282.748*S_f_array(i)) - (200*(Q_f_array(i)/(24*264.172)));
            if var > LHS_SVM
                LHS_SVM = var;
            end
        end
        %%%%%%%%%%%%%%
    
        % Constraint Variables
        rr_max = max(rr_collective); 
        rr_min = min(rr_collective); 
        Q_f_max = Q_f_input; 
        Q_f_min = Q_f; 
        Q_c_max = max([Q_c_collection;Q_c_collection2]); 
        Q_c_min = min([Q_c_collection;Q_c_collection2]); 
        Q_p_max = max([Q_p_collection;Q_p_collection2]); 
        Q_p_min = min([Q_p_collection;Q_p_collection2]);
    end
end

function E_swht = turbine(V_swht,h,rho_sw,g,eta_ht)    
    E_swht_joules = rho_sw*V_swht*g*h*eta_ht; % joules
    E_swht = E_swht_joules/(3.6e6); % kWh
end

function E_htRO = energy_recapture(V_oRO,h,eta_ROio,rho_oRO,g,eta_ht)    
    E_htRO_joules = rho_oRO*V_oRO*g*h*eta_ht*eta_ROio; % joules
    E_htRO = E_htRO_joules/(3.6e6); % kWh 
end

function [rho_ht,S_ht,V_ht] = mixing(V_oRO,V_swht,rho_oRO,S_oRO,gamma_RO,eta_RO,rho_sw,S_sw)
    rho_ht = ((rho_sw*(1-gamma_RO)) + (gamma_RO*rho_oRO*(1-eta_RO)))/(1-(eta_RO*gamma_RO));
    S_ht = ((S_sw*rho_sw*(1-gamma_RO)) + (gamma_RO*S_oRO*rho_oRO*(1-eta_RO)))/((rho_sw*(1-gamma_RO)) + (gamma_RO*rho_oRO*(1-eta_RO)));
    V_ht = V_swht + V_oRO;
end

function E_er = energy_summation(E_rd,E_swht,E_htRO)
    E_er = E_rd+E_swht+E_htRO; % kWh
end

function C = RO_cost(N_e1,N_e2,N_pv1,N_pv2,M_cost_e1,M_cost_e2)
    C1 = N_pv1*(M_cost_e1 + ((N_e1-1)*M_cost_e2));
    C2 = N_pv2*(M_cost_e1 + ((N_e2-1)*M_cost_e2));
        
    C = C1 + C2;
end

% Miscellaneous Functions
function S_brine = salinity_from_osmotic(pi,rho_fw,T,R,MW_salt_avg)
    % Get Temperature and freshwater density
    t = T; % celsius
    T = (T)+273; % kelvin (from celsius)
    
    % Osmotic Coefficient Correlation (good for S>10)
    a1 = 8.9453233003*10^-1;
    a2 = 4.1560737424*10^-4;
    a3 = -4.6262121398*10^-6;
    a4 = 2.2211195897*10^-11;
    a5 = -1.1445456438*10^-4;
    a6 = -1.4783462366*10^-6;
    a7 = -1.3526263499*10^-11;
    a8 = 7.0132355546*10^-6;
    a9 = 5.6960486681*10^-8;
    a10 = -2.8624032584*10^-10;
    
    syms S real % g/kg
    
    phi = a1+(a2*t)+(a3*(t^2))+(a4*(t^4))+(a5*S)+(a6*S*t)+(a7*S*(t^3))+...
        (a8*(S^2))+(a9*t*(S^2))+(a10*(S^2)*(t^2));
    
    % Converting osmotic pressure units
    pi = pi/145.038; % MPa (from psi)
    
    % Getting numerical value for S
    eqn = pi == ((phi*R*T*rho_fw)/(10^6))*(1000/MW_salt_avg)*...
        (S/(1000-S));
    S_num = double(solve(eqn,S));
    
    S_brine = S_num;
end

function rho_brine = density_from_salinity(S,P,rho_fw,T,P_p)
    t = T; % celsius
    P_o = P_p/145.038; % MPa (from psi)
    
    % Density at Atmospheric
    b1 = 8.020*10^2;
    b2 = -2.001;
    b3 = 1.677*10^-2;
    b4 = -3.060*10^-5;
    b5 = -1.613*10^-5;
    
    S_kg_kg = S/1000; % kg/kg
    
    rho_sw_atm = rho_fw+((b1*S_kg_kg)+(b2*S_kg_kg*t)+(b3*S_kg_kg*(t^2))+...
        (b4*S_kg_kg*(t^3))+(b5*(S_kg_kg^2)*(t^2))); % kg/m^3
    
    % Other Correction Term
    c1 = 5.0792*10^-4;
    c2 = -3.4168*10^-6;
    c3 = 5.6931*10^-8;
    c4 = -3.7263*10^-10;
    c5 = 1.4465*10^-12;
    c6 = -1.7058*10^-15;
    c7 = -1.3389*10^-6;
    c8 = 4.8603*10^-9;
    c9 = -6.8039*10^-13;
    
    d1 = -1.1077*10^-6;
    d2 = 5.5584*10^-9;
    d3 = -4.2539*10^-11;
    d4 = 8.3702*10^-9;
    
    P = P/145.038; % MPa (from psi)
    
    Fp = exp(((P-P_o)*(c1+(c2*t)+(c3*(t^2))+(c4*(t^3))+(c5*(t^4))+...
        (c6*(t^5))+(S*(d1+(d2*t)+(d3*(t^2))))))+((0.5*((P^2)-(P_o^2)))*...
        (c7+(c8*t)+(c9*(t^3))+(d4*S))));
    
    rho_brine = rho_sw_atm*Fp; % kg/m^3
end

% RO helper function (handles symbolic Q_p stuff)
function [Q_p_num,P_c_num,c_c_num,Q_c_num,rr_num,oopsies] = symbolic_Q_p_stuff(c_f,T,Q_f,Am,P_f,P_p)
    % Set up symbolic Q_p
    syms Q_p real
        
    % Getting Permeate Osmotic Pressure
    R = R_value(Q_p);            % -
    c_p = c_f*(1-R);             % mol/L
    pi_p = 1.12*(273+T)*(2*c_p); % psi
    
    % Getting Osmotic Pressure along Membrane
    rr = Q_p/Q_f;                    % -
    pf = exp(0.7*rr);                % -
    Q_c = Q_f-Q_p;                   % gpd
    c_c = ((Q_f*c_f)-(Q_p*c_p))/Q_c; % mol/L
    c_fc = (c_f+c_c)/2;              % mol/L
    c_m = (pf*(c_fc-c_p))+c_p;       % mol/L
    pi_m = 1.12*(273+T)*(2*c_m);     % psi
    
    % Getting Pressure Drop along Element
    Q_fc = 0.5*(Q_f+Q_c);           % gpd
    dP_fc = 0.01*((Q_fc/1440)^1.7); % psi
    
    % Getting Aw
    Aw = Aw_value(Q_p,c_p); % gfd/psi
    
%     disp('Here is another one')
%     disp(P_f)
%     Q_f_for_display = Q_f/6340.128;
%     disp(Q_f_for_display)

    % Getting Numerical Value for Q_p
    eqn = Q_p == Aw*Am*(P_f-(0.5*dP_fc)-P_p-pi_m+pi_p);
    %tic
    Q_p_num = double(vpasolve(eqn,Q_p,[60 8370])); % gpd
    %toc

%     disp(Q_p_num)

    if isempty(Q_p_num) == 1
        oopsies = 1;
        Q_p_num = 0;
        P_c_num = 0;
        c_c_num = 0;
        Q_c_num = 0;
        rr_num = 0;
    else
        oopsies = 0;
        % Subbing in Q_p to element outputs
        P_c_num = P_f-double(subs(dP_fc,Q_p,Q_p_num)); % psi
        c_c_num = double(subs(c_c,Q_p,Q_p_num));       % mol/L
        Q_c_num = double(subs(Q_c,Q_p,Q_p_num));       % gpd
        rr_num = double(subs(rr,Q_p,Q_p_num));         % -
    end
end

% Function for calculating R (symbolic)
function R = R_value(Q_p)
    R = 1.0034 - (12.2124*(Q_p^(-0.8122)));
end

% Function for calculating Aw (symbolic)
function Aw = Aw_value(Q_p,c_p)
    term1 = -85.71804787;
    term2 = 3.3699665311*log(c_p);
    term3 = 53.163612047*log(Q_p);
    term4 = 3.2938449078*log(c_p)*log(c_p);
    term5 = 0.726912142*log(c_p)*log(c_p)*log(c_p);
    term6 = 0.0516395755*log(c_p)*log(c_p)*log(c_p)*log(c_p);
    term7 = -12.22459358*log(Q_p)*log(Q_p);
    term8 = 1.1774872012*log(Q_p)*log(Q_p)*log(Q_p);
    term9 = -0.041266007*log(Q_p)*log(Q_p)*log(Q_p)*log(Q_p);
    
    term = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9;
    Aw = exp(term);
end
