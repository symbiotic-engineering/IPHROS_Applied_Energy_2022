% Matthew Haefner
% SEA Lab
% Remaking of python hyperplane figures in matlab
% 8/14/23

%% Equations

% Getting the original svm plane equation (these match with lines 120-125
% of Figure_6_hyperplane.py)
% V12 = [200;4;12.5671];
% V13 = [500;11;30.0041];
% C = cross(V12,V13) = [-18.22208934,282.74837627,200];

% Original Plane Equation (the point used is P1 (from line 116 of
% Figure_6_hyperplane.py))
% -18.22208934(P_f-500) + 282.74837627(S_f-40) + 200(Q_f+10.031609776253827) = 0

% New (shifted) plane equation
% -18.22208934(P_f-1033) + 282.74837627(S_f-54) + 200(Q_f+11.1) = 0

% Note: above eqns are with pressure in units in psi

%% Getting Data Points

data = readmatrix("Haefner_Figure_6_Data.xlsx");
data_trimmed = data(20:end,:);

P_f = data_trimmed(:,2)*6894.76; % Pa (from psi)
S_f = data_trimmed(:,3); % g/kg
Q_f = data_trimmed(:,6); % m^3/hr
good_bad = data_trimmed(:,13); % -

P_f_good = P_f(good_bad ~= 0);
S_f_good = S_f(good_bad ~= 0);
Q_f_good = Q_f(good_bad ~= 0);
P_f_bad = P_f(good_bad == 0);
S_f_bad = S_f(good_bad == 0);
Q_f_bad = Q_f(good_bad == 0);

%% Symbolic expressions for planes
% old means not shifted, new is shifted
% see above for origin of plane equations

syms P_f_sym S_f_sym

Q_f_old = ((192.0280447 + 18.2217*(P_f_sym/6894.76) - 282.73*S_f_sym)/200); % m^3/hr
Q_f_new = ((-1334.934 + 18.222*(P_f_sym/6894.76) - 282.748*S_f_sym)/200); % m^3/hr

%% Plotting

figure(1)
scatter3(P_f_good,S_f_good,Q_f_good,'LineWidth',2)
hold on
scatter3(P_f_bad,S_f_bad,Q_f_bad,'LineWidth',2)
fmesh(Q_f_old,'EdgeColor','k','LineWidth',3)
fmesh(Q_f_new,'EdgeColor','m','LineWidth',3)
xlabel({'$$Feed \; Pressure$$', '$$P_f \; [Pa]$$'},'Interpreter','latex','FontSize',14)
ylabel({'$$Feed \; Salinity$$', '$$S_f \; [g/kg]$$'},'Interpreter','latex','FontSize',14)
zlabel({'$$Feed \; Volumetric \; Flowrate$$', '$$Q_f \; [m^3/hr]$$'},'Interpreter','latex','FontSize',14)
legend('Legal Input Configuration','Illegal Output Configuration','SVM Classification Plane','Shifted SVM Classification Plane','FontSize',14)
view(45,33.5)
