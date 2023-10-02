% Matthew Haefner
% SEA Lab
% Optimization Results Analysis
% 8/10/23

%% Loading in Optimization Results

data = readmatrix('..\data\Haefner_Figure_12_13_14_15_Table_5_Data.xlsx','Sheet','Sheet1');
x_pareto = data(1:70,1:8);
obj_vals = data(1:70,14:16);
results = [x_pareto,obj_vals];


%% Post Processing

E_er_array = obj_vals(:,1);
V_fwRO_array = obj_vals(:,2);
eta_RO_array = obj_vals(:,3);
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
eta_RO_array_normalized = zeros(70,1);
E_r_array_normalized = zeros(70,1);
gamma_array_normalized = zeros(70,1);
gamma_RO_array_normalized = zeros(70,1);
h_L_array_normalized = zeros(70,1);
N_e1_array_normalized = zeros(70,1);
N_e2_array_normalized = zeros(70,1);
N_pv1_array_normalized = zeros(70,1);
N_pv2_array_normalized = zeros(70,1);

% Design variable bounds
min_E_r = 1; max_E_r = 100e6;
min_gamma = 0.01; max_gamma = 0.99;
min_gamma_RO = 0.01; max_gamma_RO = 0.99;
min_h_L = 240; max_h_L = 821;
min_N_e1 = 1; max_N_e1 = 8;
min_N_e2 = 0; max_N_e2 = 8;
min_N_pv1 = 1; max_N_pv1 = 1.86e6;
min_N_pv2 = 0; max_N_pv2 = 1.86e6;
min_E_er = min(E_er_array); max_E_er = max(E_er_array);
min_V_fwRO = min(V_fwRO_array); max_V_fwRO = max(V_fwRO_array);
min_eta_RO = min(eta_RO_array); max_eta_RO = max(eta_RO_array);

for i = 1:70
    E_er_array_normalized(i) = (E_er_array(i) - min_E_er)/(max_E_er-min_E_er);
    V_fwRO_array_normalized(i) = (V_fwRO_array(i) - min_V_fwRO)/(max_V_fwRO-min_V_fwRO);
    eta_RO_array_normalized(i) = (eta_RO_array(i) - min_eta_RO)/(max_eta_RO-min_eta_RO);
    E_r_array_normalized(i) = (E_r_array(i) - min_E_r)/(max_E_r-min_E_r);
    gamma_array_normalized(i) = (gamma_array(i) - min_gamma)/(max_gamma-min_gamma);
    gamma_RO_array_normalized(i) = (gamma_RO_array(i) - min_gamma_RO)/(max_gamma_RO-min_gamma_RO);
    h_L_array_normalized(i) = (h_L_array(i) - min_h_L)/(max_h_L-min_h_L);
    N_e1_array_normalized(i) = (N_e1_array(i) - min_N_e1)/(max_N_e1-min_N_e1);
    N_e2_array_normalized(i) = (N_e2_array(i) - min_N_e2)/(max_N_e2-min_N_e2);
    N_pv1_array_normalized(i) = (N_pv1_array(i) - min_N_pv1)/(max_N_pv1-min_N_pv1);
    N_pv2_array_normalized(i) = (N_pv2_array(i) - min_N_pv2)/(max_N_pv2-min_N_pv2);
end

% Compiling normalized results
results_norm = [E_r_array_normalized,gamma_array_normalized,...
    gamma_RO_array_normalized,h_L_array_normalized,...
    N_e1_array_normalized,N_e2_array_normalized,N_pv1_array_normalized,...
    N_pv2_array_normalized,E_er_array_normalized,V_fwRO_array_normalized,...
    eta_RO_array_normalized];


%% Filling out the table
[M_max_E,I_max_E] = max(E_er_array);
[M_max_V,I_max_V] = max(V_fwRO_array);
[M_max_eta_RO,I_max_eta_RO] = max(eta_RO_array);

[M_max_E_norm,I_max_E_norm] = max(E_er_array_normalized);
[M_max_V_norm,I_max_V_norm] = max(V_fwRO_array_normalized);
[M_max_eta_RO_norm,I_max_eta_RO_norm] = max(eta_RO_array_normalized);

utopia_point = [M_max_E_norm,M_max_V_norm,M_max_eta_RO_norm];

distances = zeros(70,1);
for i = 1:70
    pareto_point = [E_er_array_normalized(i),V_fwRO_array_normalized(i),eta_RO_array_normalized(i)];
    distances(i) = norm(utopia_point - pareto_point);
end

[~,I_best] = min(distances);


%% Processing data for parallel coordinates plot

% Finding the best 10% for each objective
[~,I_E] = maxk(E_er_array,7);
[~,I_V] = maxk(V_fwRO_array,7);
[~,I_eta] = maxk(eta_RO_array,7);

bin_best_E = sort(I_E);
bin_best_V = sort(I_V);
bin_best_eta_RO = sort(I_eta);
bin_best_pareto = I_best;
bins = ones(70,1);
for i = 1:7
    bins(bin_best_E(i)) = 2;
    bins(bin_best_V(i)) = 3;
    bins(bin_best_eta_RO(i)) = 4;
end
bins(bin_best_pareto) = 5;

% Rearranging rows for Parallel Coordinate Plot
results_par_coor = [results,bins];
results_par_coor = sortrows(results_par_coor,12);
bins_par_coor = results_par_coor(:,12);
results_par_coor = results_par_coor(:,1:11);
% results_par_coor = [results_norm,bins];
% results_par_coor = sortrows(results_par_coor,12);
% bins_par_coor = results_par_coor(:,12);
% results_par_coor = results_par_coor(:,1:11);

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
Qf_min_constraint = (c_collection(:,2)-3.41)*(-1);
Qf_max_constraint = c_collection(:,3)+15.5;
Qc_min_constraint = (c_collection(:,4)-3.41)*(-1);
Qc_max_constraint = c_collection(:,5)+15.5;
Qp_min_constraint = (c_collection(:,6)-0)*(-1);
Qp_max_constraint = c_collection(:,7)+1.32;
rr_min_constraint = (c_collection(:,8)-0)*(-1);
rr_max_constraint = c_collection(:,9)+0.13;
ndp_min_constraint = (c_collection(:,11)-(1e-6))*(-1);
ndp_min_constraint = ndp_min_constraint*6894.76; % to convert psi to pa

constraint_vals = [Sht_constraint,Qf_min_constraint,Qf_max_constraint,...
    Qc_min_constraint,Qc_max_constraint,Qp_min_constraint,...
    Qp_max_constraint,rr_min_constraint,rr_max_constraint,ndp_min_constraint];

%% Plotting Pareto Front
% 3D Pareto Front, Scatter
figure(1)
scatter3(E_er_array,V_fwRO_array,eta_RO_array,[],'k')
hold on
scatter3(E_er_array(I_max_E),V_fwRO_array(I_max_E),eta_RO_array(I_max_E),[],[0.8 0.4 0],'filled')
scatter3(E_er_array(I_max_V),V_fwRO_array(I_max_V),eta_RO_array(I_max_V),[],[0.8 0.6 0.7],'filled')%[0.9 0.6 0],'filled')
scatter3(E_er_array(I_max_eta_RO),V_fwRO_array(I_max_eta_RO),eta_RO_array(I_max_eta_RO),[],[0 0.6 0.5],'filled')%[0.8 0.6 0.7],'filled')
scatter3(E_er_array(I_best),V_fwRO_array(I_best),eta_RO_array(I_best),[],[0 0 0],'filled')%[0.95 0.9 0.25],'filled')
xlabel({'$$Energy \; Sent \; to$$','$$Consumer \; Per \; Day$$','$$\dot{E}_{er} \; [kWh/day]$$'}, 'Interpreter','latex','FontSize',14)
ylabel({'$$Freshwater \; Volumetric$$','$$Flowrate$$','$$\dot{V}_{fw,RO} \; [m^3/day]$$'}, 'Interpreter','latex','FontSize',14)
zlabel({'$$RO \; System \; Recovery \; Ratio$$','$$\eta_{RO} \; [-]$$'}, 'Interpreter','latex','FontSize',14)
view(60,31)

% 2D Pareto Front, Scatter, E and V
figure(2)
scatter(E_er_array,V_fwRO_array,[],[0 0.45 0.7],'LineWidth',2)
hold on
scatter(E_er_array(I_max_E),V_fwRO_array(I_max_E),[],[0.8 0.4 0],'filled')
scatter(E_er_array(I_max_V),V_fwRO_array(I_max_V),[],[0.8 0.6 0.7],'filled')%[0.9 0.6 0],'filled')
scatter(E_er_array(I_max_eta_RO),V_fwRO_array(I_max_eta_RO),[],[0 0.6 0.5],'filled')%[0.8 0.6 0.7],'filled')
scatter(E_er_array(I_best),V_fwRO_array(I_best),[],[0 0 0],'filled')%[0.95 0.9 0.25],'filled')
xlabel({'$$Energy \; Sent \; to \; Consumer \; Per \; Day$$','$$\dot{E}_{er} \; [kWh/day]$$'}, 'Interpreter','latex','FontSize',14) % 14
ylabel({'$$Freshwater \; Volumetric \; Flowrate$$','$$\dot{V}_{fw,RO} \; [m^3/day]$$'}, 'Interpreter','latex','FontSize',14) % 14

% 2D Pareto Front, Scatter, E and eta_RO
figure(3)
scatter(E_er_array,eta_RO_array,[],[0 0.45 0.7],'LineWidth',2)
hold on
scatter(E_er_array(I_max_E),eta_RO_array(I_max_E),[],[0.8 0.4 0],'filled')
scatter(E_er_array(I_max_V),eta_RO_array(I_max_V),[],[0.8 0.6 0.7],'filled')%[0.9 0.6 0],'filled')
scatter(E_er_array(I_max_eta_RO),eta_RO_array(I_max_eta_RO),[],[0 0.6 0.5],'filled')%[0.8 0.6 0.7],'filled')
scatter(E_er_array(I_best),eta_RO_array(I_best),[],[0 0 0],'filled')%[0.95 0.9 0.25],'filled')
xlabel({'$$Energy \; Sent \; to \; Consumer \; Per \; Day$$','$$\dot{E}_{er} \; [kWh/day]$$'}, 'Interpreter','latex','FontSize',14) % 14
ylabel({'$$RO \; System \; Recovery \; Ratio$$','$$\eta_{RO} \; [-]$$'}, 'Interpreter','latex','FontSize',14) % 14

% 2D Pareto Front, Scatter, V and eta_RO
figure(4)
scatter(V_fwRO_array,eta_RO_array,[],[0 0.45 0.7],'LineWidth',2)
hold on
scatter(V_fwRO_array(I_max_E),eta_RO_array(I_max_E),[],[0.8 0.4 0],'filled')
scatter(V_fwRO_array(I_max_V),eta_RO_array(I_max_V),[],[0.8 0.6 0.7],'filled')%[0.9 0.6 0],'filled')
scatter(V_fwRO_array(I_max_eta_RO),eta_RO_array(I_max_eta_RO),[],[0 0.6 0.5],'filled')%[0.8 0.6 0.7],'filled')
scatter(V_fwRO_array(I_best),eta_RO_array(I_best),[],[0 0 0],'filled')%[0.95 0.9 0.25],'filled')
xlabel({'$$Freshwater \; Volumetric \; Flowrate$$','$$\dot{V}_{fw,RO} \; [m^3/day]$$'}, 'Interpreter','latex','FontSize',14) % 14
ylabel({'$$RO \; System \; Recovery \; Ratio$$','$$\eta_{RO} \; [-]$$'}, 'Interpreter','latex','FontSize',14) % 14

% Interpolated surface plot for pareto front
F = scatteredInterpolant(E_er_array,V_fwRO_array,eta_RO_array,'natural','none');
xgr = linspace(min(E_er_array),max(E_er_array));
ygr = linspace(min(V_fwRO_array),max(V_fwRO_array));
[XX,YY] = meshgrid(xgr,ygr);
ZZ = F(XX,YY);

figure(5)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(E_er_array,V_fwRO_array,eta_RO_array,40,'k')
hold on
scatter3(E_er_array(I_max_E),V_fwRO_array(I_max_E),eta_RO_array(I_max_E),40,[0.8 0.4 0],'filled')
scatter3(E_er_array(I_max_V),V_fwRO_array(I_max_V),eta_RO_array(I_max_V),40,[0.8 0.6 0.7],'filled')%[0.9 0.6 0],'filled')
scatter3(E_er_array(I_max_eta_RO),V_fwRO_array(I_max_eta_RO),eta_RO_array(I_max_eta_RO),40,[0 0.6 0.5],'filled')%[0.8 0.6 0.7],'filled')
scatter3(E_er_array(I_best),V_fwRO_array(I_best),eta_RO_array(I_best),40,[0 0 0],'filled')%[0.95 0.9 0.25],'filled')
xlabel({'$$Energy \; Sent \; to$$','$$Consumer \; Per \; Day$$','$$\dot{E}_{er} \; [kWh/day]$$'}, 'Interpreter','latex','FontSize',14)
ylabel({'$$Freshwater \; Volumetric$$','$$Flowrate$$','$$\dot{V}_{fw,RO} \; [m^3/day]$$'}, 'Interpreter','latex','FontSize',14)
zlabel({'$$RO \; System \; Recovery \; Ratio$$','$$\eta_{RO} \; [-]$$'}, 'Interpreter','latex','FontSize',14)
hleg = legend('$Pareto \: Surface$','$Pareto \: Points$','$max(\dot{E}_{er})$','$max(\dot{V}_{fw,RO})$','$max(\eta_{RO})$','$rel. \: max(J)$');
set(hleg, 'Interpreter', 'latex');
hleg.FontSize = 14;

view(208,12)

% Parallel Coordinate Plot
figure(6)
p = parallelplot(results_par_coor);
%p.Color = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
p.Color = [0.6784 0.8471 0.902;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0 0 0];
p.DataNormalization = 'range';
p.CoordinateTickLabels = {'DVar1';'DVar2';'DVar3';'DVar4';'DVar5';'DVar6';'DVar7';'DVar8';'Obj1';'Obj2';'Obj3'};
p.Jitter = 0;
p.GroupData = bins_par_coor;
p.LegendVisible = 1;
%exportgraphics(gca,'zzzzz_legend.eps','ContentType','vector')

%% Plotting Constraint Stuff
x_con_array = linspace(1,70,70);

% Salinity
S_f_max_con_plotting = constraint_vals(:,1);
figure(7)
yline(40,'-','Max','LineWidth',2)
hold on
scatter(x_con_array,S_f_max_con_plotting,30,[0 0.45 0.7])
scatter(x_con_array(I_E),S_f_max_con_plotting(I_E),30,[0.8 0.4 0],'filled')
scatter(x_con_array(I_V),S_f_max_con_plotting(I_V),30,[0.8 0.6 0.7],'filled')
scatter(x_con_array(I_eta),S_f_max_con_plotting(I_eta),30,[0 0.6 0.5],'filled')
scatter(x_con_array(I_best),S_f_max_con_plotting(I_best),30,[0 0 0],'filled')
xlabel('$Pareto \; Solution$','Interpreter','latex','FontSize',14)
ylabel({'$$Discharge \; Salinity$$','$$S_{ht} \; [g/kg]$$'},'Interpreter','latex','FontSize',14)
ylim([35 40.5])
hleg = legend('$S_{ht, \: max}$','$Pareto \: Points$','$Top \: 10\% \: max(\dot{E}_{er})$', ...
    '$Top \: 10\% \: max(\dot{V}_{fw,RO})$','$Top \: 10\% \: max(\eta_{RO})$','$rel. \: max(J)$');
set(hleg, 'Interpreter', 'latex');
hleg.FontSize = 14;
hleg.Location = 'eastoutside';

% Q_f
Q_f_min_con_plotting = constraint_vals(:,2);
Q_f_max_con_plotting = constraint_vals(:,3);
figure(8)
yline(15.5,'-','Max','LineWidth',2)
hold on
yline(3.41,'-','Min','LineWidth',2)
scatter(x_con_array,Q_f_max_con_plotting,30,[0 0.45 0.7])
scatter(x_con_array(I_E),Q_f_max_con_plotting(I_E),30,[0.8 0.4 0],"filled")
scatter(x_con_array(I_V),Q_f_max_con_plotting(I_V),30,[0.8 0.6 0.7],"filled")
scatter(x_con_array(I_eta),Q_f_max_con_plotting(I_eta),30,[0 0.6 0.5],"filled")
scatter(x_con_array(I_best),Q_f_max_con_plotting(I_best),30,[0 0 0],'filled')
scatter(x_con_array,Q_f_min_con_plotting,30,[0 0.45 0.7],'v')
scatter(x_con_array(I_E),Q_f_min_con_plotting(I_E),30,[0.8 0.4 0],"filled",'v')
scatter(x_con_array(I_V),Q_f_min_con_plotting(I_V),30,[0.8 0.6 0.7],"filled",'v')
scatter(x_con_array(I_eta),Q_f_min_con_plotting(I_eta),30,[0 0.6 0.5],"filled",'v')
scatter(x_con_array(I_best),Q_f_min_con_plotting(I_best),30,[0 0 0],'filled','v')
xlabel('$Pareto \; Solution$','Interpreter','latex','FontSize',14)
ylabel({'$$Feed \; Volumetric \; Flowrate$$','$$Q_f \; [m^3/h]$$'},'Interpreter','latex','FontSize',14)
hleg = legend('$Q_{f, \: max}$','$Q_{f, \: min}$','$Pareto \: Points, \: max$', ...
    '$Top \: 10\% \: max(\dot{E}_{er}), \: max$','$Top \: 10\% \: max(\dot{V}_{fw,RO}), \: max$', ...
    '$Top \: 10\% \: max(\eta_{RO}), \: max$','$rel. \: max(J), \: max$','$Pareto \: Points, \: min$', ...
    '$Top \: 10\% \: max(\dot{E}_{er}), \: min$','$Top \: 10\% \: max(\dot{V}_{fw,RO}), \: min$', ...
    '$Top \: 10\% \: max(\eta_{RO}), \: min$','$rel. \: max(J), \: min$');
set(hleg, 'Interpreter', 'latex');
hleg.FontSize = 14;
hleg.Location = 'eastoutside';

% Q_c
Q_c_min_con_plotting = constraint_vals(:,4);
Q_c_max_con_plotting = constraint_vals(:,5);
figure(9)
yline(15.5,'-','Max','LineWidth',2)
hold on
yline(3.41,'-','Min','LineWidth',2)
scatter(x_con_array,Q_c_max_con_plotting,30,[0 0.45 0.7])
scatter(x_con_array(I_E),Q_c_max_con_plotting(I_E),30,[0.8 0.4 0],"filled")
scatter(x_con_array(I_V),Q_c_max_con_plotting(I_V),30,[0.8 0.6 0.7],"filled")
scatter(x_con_array(I_eta),Q_c_max_con_plotting(I_eta),30,[0 0.6 0.5],"filled")
scatter(x_con_array(I_best),Q_c_max_con_plotting(I_best),30,[0 0 0],'filled')
scatter(x_con_array,Q_c_min_con_plotting,30,[0 0.45 0.7],'v')
scatter(x_con_array(I_E),Q_c_min_con_plotting(I_E),30,[0.8 0.4 0],"filled",'v')
scatter(x_con_array(I_V),Q_c_min_con_plotting(I_V),30,[0.8 0.6 0.7],"filled",'v')
scatter(x_con_array(I_eta),Q_c_min_con_plotting(I_eta),30,[0 0.6 0.5],"filled",'v')
scatter(x_con_array(I_best),Q_c_min_con_plotting(I_best),30,[0 0 0],'filled','v')
xlabel('$Pareto \; Solution$','Interpreter','latex','FontSize',14)
ylabel({'$$Concentrate \; Volumetric \; Flowrate$$','$$Q_c \; [m^3/h]$$'},'Interpreter','latex','FontSize',14)
hleg = legend('$Q_{c, \: max}$','$Q_{c, \: min}$','$Pareto \: Points, \: max$', ...
    '$Top \: 10\% \: max(\dot{E}_{er}), \: max$','$Top \: 10\% \: max(\dot{V}_{fw,RO}), \: max$', ...
    '$Top \: 10\% \: max(\eta_{RO}), \: max$','$rel. \: max(J), \: max$','$Pareto \: Points, \: min$', ...
    '$Top \: 10\% \: max(\dot{E}_{er}), \: min$','$Top \: 10\% \: max(\dot{V}_{fw,RO}), \: min$', ...
    '$Top \: 10\% \: max(\eta_{RO}), \: min$','$rel. \: max(J), \: min$');
set(hleg, 'Interpreter', 'latex');
hleg.FontSize = 14;
hleg.Location = 'eastoutside';

% Q_p
Q_p_min_con_plotting = constraint_vals(:,6);
Q_p_max_con_plotting = constraint_vals(:,7);
figure(10)
yline(1.32,'-','Max','LineWidth',2)
hold on
yline(0,'-','Min','LineWidth',2)
scatter(x_con_array,Q_p_max_con_plotting,30,[0 0.45 0.7])
scatter(x_con_array(I_E),Q_p_max_con_plotting(I_E),30,[0.8 0.4 0],"filled")
scatter(x_con_array(I_V),Q_p_max_con_plotting(I_V),30,[0.8 0.6 0.7],"filled")
scatter(x_con_array(I_eta),Q_p_max_con_plotting(I_eta),30,[0 0.6 0.5],"filled")
scatter(x_con_array(I_best),Q_p_max_con_plotting(I_best),30,[0 0 0],'filled')
scatter(x_con_array,Q_p_min_con_plotting,30,[0 0.45 0.7],'v')
scatter(x_con_array(I_E),Q_p_min_con_plotting(I_E),30,[0.8 0.4 0],"filled",'v')
scatter(x_con_array(I_V),Q_p_min_con_plotting(I_V),30,[0.8 0.6 0.7],"filled",'v')
scatter(x_con_array(I_eta),Q_p_min_con_plotting(I_eta),30,[0 0.6 0.5],"filled",'v')
scatter(x_con_array(I_best),Q_p_min_con_plotting(I_best),30,[0 0 0],'filled','v')
xlabel('$Pareto \; Solution$','Interpreter','latex','FontSize',14)
ylabel({'$$Permeate \; Volumetric \; Flowrate$$','$$Q_p \; [m^3/h]$$'},'Interpreter','latex','FontSize',14)
hleg = legend('$Q_{p, \: max}$','$Q_{p, \: min}$','$Pareto \: Points, \: max$', ...
    '$Top \: 10\% \: max(\dot{E}_{er}), \: max$','$Top \: 10\% \: max(\dot{V}_{fw,RO}), \: max$', ...
    '$Top \: 10\% \: max(\eta_{RO}), \: max$','$rel. \: max(J), \: max$','$Pareto \: Points, \: min$', ...
    '$Top \: 10\% \: max(\dot{E}_{er}), \: min$','$Top \: 10\% \: max(\dot{V}_{fw,RO}), \: min$', ...
    '$Top \: 10\% \: max(\eta_{RO}), \: min$','$rel. \: max(J), \: min$');
set(hleg, 'Interpreter', 'latex');
hleg.FontSize = 14;
hleg.Location = 'eastoutside';

% rr
rr_min_con_plotting = constraint_vals(:,8);
rr_max_con_plotting = constraint_vals(:,9);
figure(11)
yline(0.13,'-','Max','LineWidth',2)
hold on
yline(0,'-','Min','LineWidth',2)
scatter(x_con_array,rr_max_con_plotting,30,[0 0.45 0.7])
scatter(x_con_array(I_E),rr_max_con_plotting(I_E),30,[0.8 0.4 0],"filled")
scatter(x_con_array(I_V),rr_max_con_plotting(I_V),30,[0.8 0.6 0.7],"filled")
scatter(x_con_array(I_eta),rr_max_con_plotting(I_eta),30,[0 0.6 0.5],"filled")
scatter(x_con_array(I_best),rr_max_con_plotting(I_best),30,[0 0 0],'filled')
scatter(x_con_array,rr_min_con_plotting,30,[0 0.45 0.7],'v')
scatter(x_con_array(I_E),rr_min_con_plotting(I_E),30,[0.8 0.4 0],"filled",'v')
scatter(x_con_array(I_V),rr_min_con_plotting(I_V),30,[0.8 0.6 0.7],"filled",'v')
scatter(x_con_array(I_eta),rr_min_con_plotting(I_eta),30,[0 0.6 0.5],"filled",'v')
scatter(x_con_array(I_best),rr_min_con_plotting(I_best),30,[0 0 0],'filled','v')
xlabel('$Pareto \; Solution$','Interpreter','latex','FontSize',14)
ylabel({'$$Recovery \; Ratio$$','$$rr \; [-]$$'},'Interpreter','latex','FontSize',14)
hleg = legend('$rr_{max}$','$rr_{min}$','$Pareto \: Points, \: max$', ...
    '$Top \: 10\% \: max(\dot{E}_{er}), \: max$','$Top \: 10\% \: max(\dot{V}_{fw,RO}), \: max$', ...
    '$Top \: 10\% \: max(\eta_{RO}), \: max$','$rel. \: max(J), \: max$','$Pareto \: Points, \: min$', ...
    '$Top \: 10\% \: max(\dot{E}_{er}), \: min$','$Top \: 10\% \: max(\dot{V}_{fw,RO}), \: min$', ...
    '$Top \: 10\% \: max(\eta_{RO}), \: min$','$rel. \: max(J), \: min$');
set(hleg, 'Interpreter', 'latex');
hleg.FontSize = 14;
hleg.Location = 'eastoutside';

% NDP
NDP_min_con_plotting = constraint_vals(:,10);
figure(12)
yline(0,'-','Min','LineWidth',2)
hold on
scatter(x_con_array,NDP_min_con_plotting,30,[0 0.45 0.7],'v')
scatter(x_con_array(I_E),NDP_min_con_plotting(I_E),30,[0.8 0.4 0],"filled",'v')
scatter(x_con_array(I_V),NDP_min_con_plotting(I_V),30,[0.8 0.6 0.7],"filled",'v')
scatter(x_con_array(I_eta),NDP_min_con_plotting(I_eta),30,[0 0.6 0.5],"filled",'v')
scatter(x_con_array(I_best),NDP_min_con_plotting(I_best),30,[0 0 0],'filled','v')
xlabel('$Pareto \; Solution$','Interpreter','latex','FontSize',14)
ylabel({'$$Net \; Driving \; Pressure$$','$$NDP \; [Pa]$$'},'Interpreter','latex','FontSize',14)
hleg = legend('$NDP_{min}$','$Pareto \: Points, \: min$', ...
    '$Top \: 10\% \: max(\dot{E}_{er}), \: min$','$Top \: 10\% \: max(\dot{V}_{fw,RO}), \: min$', ...
    '$Top \: 10\% \: max(\eta_{RO}), \: min$','$rel. \: max(J), \: min$');
set(hleg, 'Interpreter', 'latex');
hleg.FontSize = 14;
hleg.Location = 'eastoutside';

%% Functions
% Matthew Haefner
% Applied Energy
% IPHROS MOGA

% Note: all V and E values are m^3/day and kWh/day

% Objective Function
function J = objective(x)
    % General Parameters
    g = 9.81;             % gravitational acceleration [m/s^2]
    rho_sw = 1023.6;      % generic density of saltwater at 1 atm and 25 degrees celsius [kg/m^3]
    rho_fw = 996.9;       % generic density of freshwater at 1 atm and 25 degrees celsius [kg/m^3]
    eta_hp = 0.894;       % pump-side efficiency
    eta_ht = 0.894;       % turbine-side efficiency [-]
    T = 25.0;             % temperature [degrees celsius]
    S_sw = 35.0;          % generic seawater salinity [g salt/kg seawater]
    M_salt = 58.44;       % molar mass of NaCl [g/mol]
    P_p = 14.696;         % permeate pressure, atmospheric for now [psi]
    R = 8.3145;           % universal gas constant [J/(mol*K)]
    MW_salt_avg = 30.548; % average molar mass of NaCl [g/mol]
    M_cost_e1 = 1.7;      % Proxy for membrane cost of 1st element in RO unit [$/unit]
    M_cost_e2 = 1;        % Proxy for membrane cost of other elements in RO unit [$/unit] 

    % Load in RO NN Weights
    load('..\data\ro_weights.mat','weights1','weights2','weights3','weights4','weights5','weights6','weights7','weights8','weights9','weights10')
    RO_weights.weights1 = weights1; RO_weights.weights2 = weights2;
    RO_weights.weights3 = weights3; RO_weights.weights4 = weights4;
    RO_weights.weights5 = weights5; RO_weights.weights6 = weights6;
    RO_weights.weights7 = weights7; RO_weights.weights8 = weights8;
    RO_weights.weights9 = weights9; RO_weights.weights10 = weights10;

    % Unpack x
    E_r = x(1); gamma = x(2); gamma_RO = x(3); h = x(4);
    N_e1 = x(5); N_e2 = x(6); N_pv1 = x(7); N_pv2 = x(8);
    
    % Go through the relevant modules
    [E_rp,E_rd] = initial_energy_division(E_r,gamma);
    V_wp = pump(E_rp,h,eta_hp,rho_sw,g);
    [V_wRO, V_swht, Q_f] = mountaintop_reservoir(V_wp,gamma_RO,N_pv1);
    [V_fwRO,V_oRO,S_oRO,rho_oRO,eta_ROio,eta_RO,rr_max,rr_min,Q_f_max,Q_f_min,Q_c_max,Q_c_min,Q_p_max,Q_p_min,NDP_min] = ro_system(h,Q_f,N_e1,N_e2,N_pv1,N_pv2,rho_sw,rho_fw,T,S_sw,M_salt,P_p,R,MW_salt_avg,RO_weights);

    E_swht = turbine(V_swht,h,rho_sw,g,eta_ht);
    E_htRO = energy_recapture(V_oRO,h,eta_ROio,rho_oRO,g,eta_ht);
    E_er = energy_summation(E_rd,E_swht,E_htRO);
    %C = RO_cost(N_e1,N_e2,N_pv1,N_pv2,M_cost_e1,M_cost_e2);

    E_er_min = -E_er; % maximizing (minimizing the negative)
    V_fwRO_min = -V_fwRO; % maximizing (minimizing the negative)
    era_RO_min = -eta_RO; % maximizing (minimizing the negative)
    
    % Compiling
    J = [E_er_min,V_fwRO_min,era_RO_min];
end

% Constraint Function
function [c,ceq] = constraints(x)
    % General Parameters
    g = 9.81;             % gravitational acceleration [m/s^2]
    rho_sw = 1023.6;      % generic density of saltwater at 1 atm and 25 degrees celsius [kg/m^3]
    rho_fw = 996.9;       % generic density of freshwater at 1 atm and 25 degrees celsius [kg/m^3]
    eta_hp = 0.894;       % pump-side efficiency
    eta_ht = 0.894;       % turbine-side efficiency [-]
    T = 25.0;             % temperature [degrees celsius]
    S_sw = 35.0;          % generic seawater salinity [g salt/kg seawater]
    M_salt = 58.44;       % molar mass of NaCl [g/mol]
    P_p = 14.696;         % permeate pressure, atmospheric for now [psi]
    R = 8.3145;           % universal gas constant [J/(mol*K)]
    MW_salt_avg = 30.548; % average molar mass of NaCl [g/mol]
    M_cost_e1 = 1.7;      % Proxy for membrane cost of 1st element in RO unit [$/unit]
    M_cost_e2 = 1;        % Proxy for membrane cost of other elements in RO unit [$/unit] 
    
    % Bring in RO NN Weights
    load('..\data\ro_weights.mat','weights1','weights2','weights3','weights4','weights5','weights6','weights7','weights8','weights9','weights10')
    RO_weights.weights1 = weights1; RO_weights.weights2 = weights2;
    RO_weights.weights3 = weights3; RO_weights.weights4 = weights4;
    RO_weights.weights5 = weights5; RO_weights.weights6 = weights6;
    RO_weights.weights7 = weights7; RO_weights.weights8 = weights8;
    RO_weights.weights9 = weights9; RO_weights.weights10 = weights10;

    % Unpack x
    E_r = x(1); gamma = x(2); gamma_RO = x(3); h = x(4);
    N_e1 = x(5); N_e2 = x(6); N_pv1 = x(7); N_pv2 = x(8);
    
    % Go through the relevant modules
    [E_rp,E_rd] = initial_energy_division(E_r,gamma);
    V_wp = pump(E_rp,h,eta_hp,rho_sw,g);
    [V_wRO, V_swht, Q_f] = mountaintop_reservoir(V_wp,gamma_RO,N_pv1);
    [V_fwRO,V_oRO,S_oRO,rho_oRO,eta_ROio,eta_RO,rr_max,rr_min,Q_f_max,Q_f_min,Q_c_max,Q_c_min,Q_p_max,Q_p_min,NDP_min] = ro_system(h,Q_f,N_e1,N_e2,N_pv1,N_pv2,rho_sw,rho_fw,T,S_sw,M_salt,P_p,R,MW_salt_avg,RO_weights);

    [rho_ht,S_ht,V_ht] = mixing(V_oRO,V_swht,rho_oRO,S_oRO,gamma_RO,eta_RO,rho_sw,S_sw);

    % Salinity constraint
    S_ht_con = S_ht-40;

    % Feed flow rate to pressure vessel constraint
    Q_f_min_con = 3.41-Q_f_min;
    Q_f_max_con = Q_f_max-15.5;

    % Concentrate flow rate constraint (gpd constraint)
    Q_c_min_con = 3.41-Q_c_min;
    Q_c_max_con = Q_c_max-15.5;

    % Permeate flow rate constraint (gpd constraint)
    Q_p_min_con = -Q_p_min;
    Q_p_max_con = Q_p_max-1.32;

    % Recovery rate for element constraints
    rr_min_con = -rr_min;
    rr_max_con = rr_max-0.13;

    % If N_e2 = 0, then N_pv2 = 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Different version of constraint %%%%%%%%%%%
    %e2_pv2_con = abs(min([N_e2,1]) - min([N_pv2,1]));
    e2_pv2_con = 2-((N_pv2^N_e2)+(N_e2^N_pv2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% New Constraint %%%%%%%%%%%%%%%
    % Equation for SVM plane
    %SVM_con = LHS_SVM;
    NDP_con = (1e-6)-NDP_min;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compiling
    %c = [S_ht_con;Q_f_min_con;Q_f_max_con;Q_c_min_con;Q_c_max_con;Q_p_min_con;Q_p_max_con;rr_min_con;rr_max_con;e2_pv2_con;SVM_con];
    c = [S_ht_con;Q_f_min_con;Q_f_max_con;Q_c_min_con;Q_c_max_con;Q_p_min_con;Q_p_max_con;rr_min_con;rr_max_con;e2_pv2_con;NDP_con];
    ceq = [];
end

% Modules
function [E_rp,E_rd] = initial_energy_division(E_r,gamma)
    E_rp = gamma*E_r; % kWh/day
    E_rd = (1-gamma)*E_r; % kWh/day
end

function V_wp = pump(E_rp,h,eta_hp,rho_sw,g)
    E_rp = E_rp*(3.6*10^6); % joules/day (from kWh/day)
    V_wp = E_rp*eta_hp/(rho_sw*g*h); % m^3/day
end

function [V_wRO,V_swht,Q_f] = mountaintop_reservoir(V_wp,gamma_RO,N_pv1)
    V_wRO = gamma_RO*V_wp; % m^3/day
    V_swht = (1-gamma_RO)*V_wp; % m^3/day
    Q_f = (V_wRO/24)/N_pv1; % m^3/hr (from m^3/day)
end

function [V_fwRO,V_oRO,S_oRO,rho_oRO,eta_ROio,eta_RO,rr_max,rr_min,Q_f_max,Q_f_min,Q_c_max,Q_c_min,Q_p_max,Q_p_min,NDP_min] = ro_system(h,Q_f_input,N_e1,N_e2,N_pv1,N_pv2,rho_sw,rho_fw,T,S_sw,M_salt,P_p,R,MW_salt_avg,RO_weights)

    % Initialize violation variable
    violation = 0; % Will equal one if input to NN is outside the WAVE simulation bounds

    % Collection Lists
    Q_p_collection = zeros(1,N_e1); % m^3/hr
    P_c_collection = zeros(1,N_e1); % psi
    c_c_collection = zeros(1,N_e1); % mol/L
    Q_c_collection = zeros(1,N_e1); % m^3/hr
    rr_collection = zeros(1,N_e1);  % -
    NDP_collection = zeros(1,N_e1); % psi
    Q_f_collection = zeros(1,N_e1); % m^3/hr
    
    % Initialize 1st stage element counter
    i1 = 1;
    
    % First stage
    while i1 <= N_e1 && violation ~= 1
        % Beginning with feed conditions into an element
        if i1 == 1
            % First element
            P_f = 0.433*(h*3.281)*(rho_sw/rho_fw); % psi
            c_f = (S_sw*rho_sw)/(1000*M_salt); % mol/L
            Q_f = Q_f_input; % m^3/hr
        else
            % Feed conditions from previous element
            P_f = P_c_collection(i1-1); % psi
            Q_f = Q_c_collection(i1-1); % m^3/hr
            c_f = c_c_collection(i1-1); % mol/L
        end

        pi_f = 1.12*(273+T)*2*c_f; % psi
        S_f = salinity_from_osmotic(pi_f,rho_fw,T,R,MW_salt_avg); % g/kg

        % Solving for Q_p and other related variables %%%%%%%%%%%%%%%%%%%
        [Q_p,P_c,c_c,Q_c,rr,NDP,violation] = Q_p_NN(P_f,S_f,Q_f,P_p,c_f,T,RO_weights);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Storing element outputs for further calculations
        Q_p_collection(i1) = Q_p; % m^3/hr
        P_c_collection(i1) = P_c; % psi
        c_c_collection(i1) = c_c; % mol/L
        Q_c_collection(i1) = Q_c; % gpd
        rr_collection(i1) = rr;   % -
        NDP_collection(i1) = NDP; % psi
        Q_f_collection(i1) = Q_f; % m^3/hr
        
        % Increasing counter
        i1 = i1 + 1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up storage variables for second stage
    Q_p_collection2 = zeros(1,N_e2); % m^3/hr
    P_c_collection2 = zeros(1,N_e2); % psi
    c_c_collection2 = zeros(1,N_e2); % mol/L
    Q_c_collection2 = zeros(1,N_e2); % m^3/hr
    rr_collection2 = zeros(1,N_e2);  % - 
    NDP_collection2 = zeros(1,N_e2); % psi
    Q_f_collection2 = zeros(1,N_e2); % m^3/hr

    % Second stage
    if N_e2 == 0 || N_pv2 == 0
        i2 = N_e2+1;
    else        
        % initialize second stage element counter
        i2 = 1;
    end

    while i2 <= N_e2 && violation ~= 1
        if i2 == 1
            % Feed conditions from previous stage (feed flowrate reconfigured)
            P_f = P_c_collection(end); % psi
            c_f = c_c_collection(end); % mol/L
            
            V_intermediate = Q_c_collection(end)*N_pv1; % m^3/hr
            Q_f = V_intermediate/N_pv2; % m^3/hr
        else
            % Feed conditions from previous element
            P_f = P_c_collection2(i2-1); % psi
            Q_f = Q_c_collection2(i2-1); % m^3/hr
            c_f = c_c_collection2(i2-1); % mol/L
        end

        pi_f = 1.12*(273+T)*2*c_f; % psi
        S_f = salinity_from_osmotic(pi_f,rho_fw,T,R,MW_salt_avg); % g/kg
        
        % Solving for Q_p and other numeric values %%%%%%%%%%%%%%%%%%%%%%%
        [Q_p,P_c,c_c,Q_c,rr,NDP,violation] = Q_p_NN(P_f,S_f,Q_f,P_p,c_f,T,RO_weights);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Storing element outputs for further calculations
        Q_p_collection2(i2) = Q_p; % m^3/hr
        P_c_collection2(i2) = P_c; % psi
        c_c_collection2(i2) = c_c; % mol/L
        Q_c_collection2(i2) = Q_c; % m^3/hr
        rr_collection2(i2) = rr;   % -
        NDP_collection2(i2) = NDP; % psi
        Q_f_collection2(i2) = Q_f; % m^3/hr
        
        % Increasing counter
        i2 = i2 + 1;
    end

    % If there was a violation, return zeros for dependent variables we want to be high, 1e6 for dependent variables we want to be low, 
    % and extreme constraint values
    if violation == 1
        % Dependent vars
        V_fwRO = 0; V_oRO = 0; S_oRO = 1e6; rho_oRO = 1e6; eta_ROio = 0; eta_RO = 0;

        % Constraint vars
        rr_max = 1e6; rr_min = -1e6; Q_f_max = 1e6; Q_f_min = -1e6; 
        Q_c_max = 1e6; Q_c_min = -1e6; Q_p_max = 1e6; Q_p_min = -1e6;
        NDP_min = -1e6;
    else
        % Processing for RO System Output Values
        V_fwRO = (N_pv1*sum(Q_p_collection) + N_pv2*sum(Q_p_collection2))*24; % m^3/day (from m^3/hr)
        
        brine_collective = [Q_c_collection,Q_c_collection2];
        if N_pv2 == 0 || N_e2 == 0
            V_oRO = (brine_collective(end)*N_pv1)*24;  % m^3/day (from m^3/hr)
        else
            V_oRO = (brine_collective(end)*N_pv2)*24;  % m^3/day (from m^3/hr)
        end
        
        c_c_collective = [c_c_collection,c_c_collection2];
        pi_c = 1.12*(273+T)*2*c_c_collective(end); % psi
        S_oRO = salinity_from_osmotic(pi_c,rho_fw,T,R,MW_salt_avg); % g/kg
    
        P_c_collective = [P_c_collection,P_c_collection2];
        rho_oRO = density_from_salinity(S_oRO,P_c_collective(end),rho_fw,T,P_p); % kg/m^3
        
        P_f_1st_element = 0.433*(h*3.281)*(rho_sw/rho_fw);                    % psi
        eta_ROio = 1-((P_f_1st_element-P_c_collective(end))/P_f_1st_element); % -
    
        rr_collective = [rr_collection,rr_collection2];
        eta_RO_product_term = 1;
        i_eta_RO = 1;
        N_e = N_e1 + N_e2;
        while i_eta_RO <= N_e
            eta_RO_product_term = eta_RO_product_term*(1-rr_collective(i_eta_RO));
            i_eta_RO = i_eta_RO + 1;
        end
        eta_RO = 1-eta_RO_product_term; % -
    
        % Constraint Variables
        rr_max = max(rr_collective); 
        rr_min = min(rr_collective); 
        Q_f_max = max([Q_f_collection,Q_f_collection2]); 
        Q_f_min = min([Q_f_collection,Q_f_collection2]); 
        Q_c_max = max([Q_c_collection,Q_c_collection2]); 
        Q_c_min = min([Q_c_collection,Q_c_collection2]); 
        Q_p_max = max([Q_p_collection,Q_p_collection2]); 
        Q_p_min = min([Q_p_collection,Q_p_collection2]);
        NDP_min = min([NDP_collection,NDP_collection2]);
    end
end

function E_swht = turbine(V_swht,h,rho_sw,g,eta_ht)    
    E_swht_joules = rho_sw*V_swht*g*h*eta_ht; % joules/day
    E_swht = E_swht_joules/(3.6e6); % kWh/day
end

function E_htRO = energy_recapture(V_oRO,h,eta_ROio,rho_oRO,g,eta_ht)    
    E_htRO_joules = rho_oRO*V_oRO*g*h*eta_ht*eta_ROio; % joules/day
    E_htRO = E_htRO_joules/(3.6e6); % kWh/day
end

function [rho_ht,S_ht,V_ht] = mixing(V_oRO,V_swht,rho_oRO,S_oRO,gamma_RO,eta_RO,rho_sw,S_sw)
    rho_ht = ((rho_sw*(1-gamma_RO)) + (gamma_RO*rho_oRO*(1-eta_RO)))/(1-(eta_RO*gamma_RO)); % kg/m^3
    S_ht = ((S_sw*rho_sw*(1-gamma_RO)) + (gamma_RO*S_oRO*rho_oRO*(1-eta_RO)))/((rho_sw*(1-gamma_RO)) + (gamma_RO*rho_oRO*(1-eta_RO))); % g/kg
    V_ht = V_swht + V_oRO; % m^3/day
end

function E_er = energy_summation(E_rd,E_swht,E_htRO)
    E_er = E_rd+E_swht+E_htRO; % kWh/day
end

function C = RO_cost(N_e1,N_e2,N_pv1,N_pv2,M_cost_e1,M_cost_e2)
    C1 = N_pv1*(M_cost_e1 + ((N_e1-1)*M_cost_e2)); % $ (proxy)
    C2 = N_pv2*(M_cost_e1 + ((N_e2-1)*M_cost_e2)); % $ (proxy)
        
    C = C1 + C2; % $ (proxy)
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
    
    S_brine = S_num; % g/kg
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

% RO helper function 
function [Q_p,P_c,c_c,Q_c,rr,NDP,violation] = Q_p_NN(P_f,S_f,Q_f,P_p,c_f,T,RO_weights)
    % Use NN to get Q_p from P_f, S_f, and Q_f, first checking the inputs
    if P_f < 300 || P_f > 1200 || S_f < 35 || S_f > 55 || Q_f < 3.41 || Q_f > 15.5
        violation = 1;
        Q_p = -1; P_c = -1; c_c = -1; Q_c = -1; rr = -1; NDP = -1;
    else
        violation = 0;
        x = [P_f,S_f,Q_f];
        Q_p = forward_prop(x,RO_weights); % m^3/hr
            
        % Getting Permeate Osmotic Pressure
        R_eqn = 1.0034 - (0.00997*(Q_p^(-0.8122))); % -
        c_p = c_f*(1-R_eqn);                        % mol/L
        pi_p = 1.12*(273+T)*(2*c_p);            % psi
        
        % Getting Osmotic Pressure along Membrane
        rr = Q_p/Q_f;                    % -
        pf = exp(0.7*rr);                % -
        Q_c = Q_f-Q_p;                   % m^3/hr
        c_c = ((Q_f*c_f)-(Q_p*c_p))/Q_c; % mol/L
        c_fc = (c_f+c_c)/2;              % mol/L
        c_m = (pf*(c_fc-c_p))+c_p;       % mol/L
        pi_m = 1.12*(273+T)*(2*c_m);     % psi
        
        % Getting Pressure Drop along Element
        Q_fc = 0.5*(Q_f+Q_c);           % m^3/hr
        dP_fc = 0.01*((Q_fc*(264.172/60))^1.7); % psi
    
        % Getting Concentrate Pressure
        P_c = P_f-dP_fc; % psi
    
        % Net Driving Pressure
        P_fc = P_f-(0.5*dP_fc);         % psi
        NDP = (P_fc-P_p) - (pi_m-pi_p); % psi
    end
end

%% Forward propagation functions
% Forward propagation function
function y = forward_prop(x,RO_weights)
    % Unpack Weights
    weights1 = RO_weights.weights1; weights2 = RO_weights.weights2;
    weights3 = RO_weights.weights3; weights4 = RO_weights.weights4;
    weights5 = RO_weights.weights5; weights6 = RO_weights.weights6;
    weights7 = RO_weights.weights7; weights8 = RO_weights.weights8;
    weights9 = RO_weights.weights9; weights10 = RO_weights.weights10;

    % First need to normalize x 
    x = normalization(x,weights1,weights2);
    
    % Set up variables
    m = size(x,1); % number of samples
    y = zeros(m,1); % one solution for each sample
    
    for i = 1:m
        a1_i = x(i,:); % one sample at a time, 1x3
        a = a1_i;
        a = single_forward_prop(a,weights3,weights4);
        a = single_forward_prop(a,weights5,weights6);
        a = single_forward_prop(a,weights7,weights8);
        y(i) = dot(a,weights9)+weights10;
    end
end

% Forward propagation helper function
function a_next = single_forward_prop(a_i,Theta,bias)
    % Set up variables
    num_from_nodes = size(Theta,2);
    a_next = zeros(1,num_from_nodes);
    
    for j = 1:num_from_nodes
        T_j = Theta(:,j); % num_to_nodesx1
        z_next_j = dot(a_i,T_j)+bias(j); % 1x1
        a_next_j = max(z_next_j,0); % ReLU function, 1x1
        a_next(j) = a_next_j;
    end
end

% Function for normalizing the data
function A_normalized = normalization(A,averages,variances)    
    % Set up variables
    m = size(A,1);
    n = size(A,2);
    A_normalized = zeros(m,n);
    
    % Normalize each element of A
    for j = 1:n
        a = averages(j);
        v = variances(j);
        s = sqrt(v);
        for i = 1:m
            A_normalized(i,j) = (A(i,j)-a)/s;
        end
    end
end