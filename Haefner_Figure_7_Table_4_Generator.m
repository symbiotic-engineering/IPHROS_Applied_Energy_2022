% Matthew Haefner
% SEA Lab
% Equation Generation for R
% 8/14/23

%% Loading in data, getting Q_p and R arrays

% Reading in Excel table
data = readmatrix('Haefner_WAVE_Simulations_Trimmed.xlsx');

% Getting just the Q_p and R columns
Q_p_array = data(:,11); % m^3/hr
R_array = data(:,27); % -

% Putting arrays into a table
tabularized_arrays = table(Q_p_array,R_array);

%% Building the models

% Biexponential 5p
model_fun_5p = @(b,x)b(1) + b(2)*exp(-b(3)*x(:,1)) + b(4)*exp(-b(5)*x(:,1));
beta0_5p = [0.99 -0.18 0.0019 -0.72 0.017];
mdl_5p = fitnlm(tabularized_arrays,model_fun_5p,beta0_5p);

% Power
model_fun_power = @(b,x)b(1) + b(2)*x(:,1).^b(3);
beta0_power = [-50 500 -1];
mdl_power = fitnlm(tabularized_arrays,model_fun_power,beta0_power);

% Rational
model_fun_rat = @(p,x)((p(1)*x(:,1).^2) + (p(2)*x(:,1)) + p(3))./(x(:,1) + p(4));
beta0_rat = [0 1 -10 20];
mdl_rat = fitnlm(tabularized_arrays,model_fun_rat,beta0_rat);

%% Comparing the fits

% Setting up arrays
x_array = linspace(min(Q_p_array),max(Q_p_array),1001);
y_5p = zeros(1,1001);
y_power = zeros(1,1001);
y_rat = zeros(1,1001);

% Getting the coefficients for the models
coeff_5p = table2array(mdl_5p.Coefficients(:,1));
coeff_power = table2array(mdl_power.Coefficients(:,1));
coeff_rat = table2array(mdl_rat.Coefficients(:,1));

% Data for plotting the fits
for i = 1:1001
    y_5p(i) = coeff_5p(1) + coeff_5p(2)*exp(-coeff_5p(3)*x_array(i)) + ...
        coeff_5p(4)*exp(-coeff_5p(5)*x_array(i));
    y_power(i) = coeff_power(1) + coeff_power(2)*x_array(i)^coeff_power(3);
    y_rat(i) = ((coeff_rat(1)*x_array(i)^2) + (coeff_rat(2)*...
        x_array(i)) + coeff_rat(3))/(x_array(i) + coeff_rat(4));
end

% Plotting (for visual inspection)
figure(1)
scatter(Q_p_array,R_array)
hold on
plot(x_array,y_5p)
xlabel('Q_p')
ylabel('R')
title('Biexponential 5p')

figure(2)
scatter(Q_p_array,R_array)
hold on
plot(x_array,y_power)
xlabel('Q_p')
ylabel('R')
title('Power')

figure(3)
scatter(Q_p_array,R_array)
hold on
plot(x_array,y_rat)
xlabel('Q_p')
ylabel('R')
title('Rational')

% R-squared values
R2_5p = mdl_5p.Rsquared.Ordinary;
R2_power = mdl_power.Rsquared.Ordinary;
R2_rat = mdl_rat.Rsquared.Ordinary;

%% Removing points based on Cook's D for power fit

% Getting Cook's distances, removing points
cook_d = table2array(mdl_power.Diagnostics(:,2));
yes_no = cook_d <= (4/2507);
new_R_array = R_array(yes_no,:);
new_Q_p_array = Q_p_array(yes_no,:);

% Building new power model
new_tabularized_arrays = table(new_Q_p_array,new_R_array);
new_mdl_power = fitnlm(new_tabularized_arrays,model_fun_power,beta0_power);

% New R-squared value
new_R2_power = new_mdl_power.Rsquared.Ordinary;

%% Final R Figure for Paper

% Setting up arrays
new_x_array = linspace(min(new_Q_p_array),max(new_Q_p_array),1001);
new_y_power = zeros(1,1001);

% Getting the coefficients for the model
new_coeff_power = table2array(new_mdl_power.Coefficients(:,1));

% Data for plotting the fits
for i = 1:1001
    new_y_power(i) = new_coeff_power(1) + new_coeff_power(2)*...
        new_x_array(i)^new_coeff_power(3);
end

% Plotting
figure(4)
scatter(new_Q_p_array,new_R_array)
hold on
plot(new_x_array,new_y_power,'k','LineWidth',2)
xlabel({'$$Permeate \; Volumetric \; Flowrate$$','$$Q_p \; [m^3/hr]$$'},'Interpreter','latex','FontSize',14)
ylabel({'$$Fractional \; Salt \; Rejection \; Rate$$','$$R \; [-]$$'},'Interpreter','latex','FontSize',14)
legend('Data','Power Fit','location','southeast','FontSize',12)
