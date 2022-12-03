% Matthew Haefner
% Applied Energy
% IPHROS MOGA

% Plot of percent error NN Q_p to WAVE Q_p 

data = readmatrix('Error_Qp_NN_WAVE_NEW.xlsx');
percent_error = data(:,5);

figure(1)
histogram(percent_error,19,'Normalization','probability')
xlabel('$Neural \; Network \; Percent \; Error$','Interpreter','latex','FontSize',14)
ylabel('$Frequency$','Interpreter','latex','FontSize',14)