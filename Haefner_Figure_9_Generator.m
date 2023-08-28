% Matthew Haefner
% SEA Lab
% Plot of percent error NN Q_p to WAVE Q_p 
% 8/14/23

data = readmatrix('Haefner_Figure_9_Data.xlsx');
percent_error = data(:,5);

figure(1)
histogram(percent_error,19,'Normalization','probability')
xlabel('$Neural \; Network \; Percent \; Error$','Interpreter','latex','FontSize',14)
ylabel('$Frequency$','Interpreter','latex','FontSize',14)