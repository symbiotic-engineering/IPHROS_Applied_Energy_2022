% Matthew Haefner
% SEA Lab
% CAPEX, OPEX, AVE, AVW (Individual Implementations) Figure Generator
% 8/14/23

%% Read in data, compile for plotting
data = readcell("Haefner_Figure_16_Data.xlsx");

CAPEX_PSH_values = cell2mat(data(19,3:6));
OPEX_PSH_values = cell2mat(data(20,3:6));
CAPEX_RO_values = cell2mat(data(21,3:6));
OPEX_RO_values = cell2mat(data(22,3:6));
AVE_values = cell2mat(data(16,3:6));
AVW_values = cell2mat(data(17,3:6));

y = [CAPEX_PSH_values;OPEX_PSH_values;CAPEX_RO_values;OPEX_RO_values;AVE_values;AVW_values];
catnames = {'CAPEX_{PSH}','OPEX_{PSH}','CAPEX_{RO}','OPEX_{RO}','AVE','AVW'};

%% Plotting
figure(1)
b = bar(y);
set(gca,'xticklabel',catnames,'YScale','log','ylim',[1e6,1e10])
ylabel('\$, \$/year','Interpreter','latex','FontSize',14)
hleg = legend('$max(E_{er})$','$max(\dot{V}_{fw,RO})$','$max(\eta_{RO})$','$rel. \: max(J)$','Location','eastoutside');
set(hleg, 'Interpreter', 'latex','FontSize',14);
b(1).FaceColor = [0.8 0.4 0];
b(2).FaceColor = [0.8 0.6 0.7];
b(3).FaceColor = [0 0.6 0.5];
b(4).FaceColor = [0 0 0];
