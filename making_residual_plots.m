% Matthew Haefner (mwh85)
% SEA Lab
% Making residual plots for Aw Eqn

%% Load in data
data = readmatrix("data_for_residual_plot.xlsx");

%% Arrange the data
row = data(:,1);
log_Aw_residuals = data(:,2);
%predicted_Aw = data(:,3);
predicted_Aw = data(:,6);
studentized_residuals = data(:,4);
actual_Aw = data(:,5);

%% Get linear fit of actual vs predicted
p = polyfit(actual_Aw,predicted_Aw,1);
x1 = linspace(0,max(predicted_Aw),1001);
y1 = polyval(p,x1);

%% Plotting
% Studentized Residuals
figure(1)
scatter(predicted_Aw,studentized_residuals)
yline(0,'-.','LineWidth',2)
pbaspect([3.5 1 1])
yticks([-5 -2.5 0 2.5 5])
xlabel('A_{w, predicted}','FontSize',14)
ylabel('Studentized Residuals','FontSize',14)
% xlabel('$A_{w, \; predicted}$','Interpreter','latex','FontSize',14)
% ylabel('$Studentized \; Residuals$','Interpreter','latex','FontSize',14)
legend('Data','y=0')

% Actual vs predicted
% figure(2)
% scatter(predicted_Aw,actual_Aw)
% hold on
% plot(x1,y1,'LineWidth',2)
% xlabel('Predicted A_w','FontSize',14)
% ylabel('Actual A_w','FontSize',14)
% title('Actual A_w vs. Predicted A_w','FontSize',16)