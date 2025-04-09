% Load data from a .mat file
filename = 'cf1_june.mat';  % Replace with your actual filename
data = load(filename);

% Access variables from the loaded file
% Example: if the file contains a variable called 'x'
x = data.x;
u = data.u;

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
lgd2 = legend(["MG1 to DNO", "MG2 to DNO", "MG3 to DNO", "DNO to MG1", "DNO to MG2", "DNO to MG3"])