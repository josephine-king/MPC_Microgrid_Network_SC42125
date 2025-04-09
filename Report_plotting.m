% Load data from a .mat file
filename = 'cf1_oct.mat';  % Replace with your actual filename
data = load(filename);

data2 = load('cf2_oct.mat');

data3 = load('cf3_oct.mat');
data4 = load('non_coop_oct.mat');

% Access variables from the loaded file
% Example: if the file contains a variable called 'x'
x = data.x;
u = data.u;

x2 = data2.x;
u2 = data2.u;

x3 = data3.x;
u3 = data3.u;

x4 = data4.x;
u4 = data4.u;

battery_size = 2000;

close all

figure(1)
stairs(100.*x(1,:)/battery_size)
hold on
stairs(100.*x(2,:)/battery_size)
hold on
stairs(100.*x(3,:)/battery_size, 'Color',[0.4,0.7,0.3])
legend(["MG1", "MG2", "MG3"])
title("Battery Percent Charged")
ylim([20 70])
ylabel("Battery Percent Charged (%)")
xlabel("Time step (hour)")
fontsize(16,"points")
grid on

figure(4)
subplot(3,1,1)
stairs(100.*x(1,:)/battery_size)
hold on
stairs(100.*x(2,:)/battery_size)
hold on
stairs(100.*x(3,:)/battery_size, 'Color',[0.4,0.7,0.3])
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
lgd2 = legend(["MG1 to DNO", "MG2 to DNO", "MG3 to DNO", "DNO to MG1", "DNO to MG2", "DNO to MG3"])

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
lgd3 = legend(["MG1 to MG2", "MG1 to MG3", "MG2 to MG1", "MG2 to MG3", "MG3 to MG1", "MG3 to MG2"])

fontsize(lgd1,12,'points')
fontsize(lgd2,12,'points')
fontsize(lgd3,12,'points')

figure(5)
mg1bar = sum(u(2,:)) + sum(u(5,:));
mg2bar = sum(u(9, :)) + sum(u(12,:));
mg3bar = sum(u(16, :)) + sum(u(19,:));

mg1bar2 = sum(u2(2,:)) + sum(u2(5,:));
mg2bar2 = sum(u2(9, :)) + sum(u2(12,:));
mg3bar2 = sum(u2(16, :)) + sum(u2(19,:));

mg1bar3 = sum(u3(2,:)) + sum(u3(5,:));
mg2bar3 = sum(u3(9, :)) + sum(u3(12,:));
mg3bar3 = sum(u3(16, :)) + sum(u3(19,:));

mg1bar4 = sum(u4(2,:)) + sum(u4(5,:));
mg2bar4 = sum(u4(9, :)) + sum(u4(12,:));
mg3bar4 = sum(u4(16, :)) + sum(u4(19,:));

bargraph = [mg1bar mg1bar2 mg1bar3 mg1bar4; mg2bar mg2bar2 mg2bar3 mg2bar4; mg3bar mg3bar2 mg3bar3 mg3bar4];
bar(bargraph)

% figure(6)
% stairs(wt(1,:) + pv(1,:) - D(1,:))
% hold on
% stairs(wt(2,:) + pv(2,:) - D(2,:))
% hold on
% stairs(wt(3,:) + pv(3,:) - D(3,:), 'Color',[0.4,0.7,0.3])
% legend(["MG1", "MG2", "MG3"])
% title("Power balance, October 15-17")
% ylabel("Power balance (kW)")
% xlabel("Time step (hour)")
% fontsize(16,"points")
% 
% figure(7)
% plot(mg1.rte*u(1,:))
% %+u(2,:)+u(3,:)+u(4,:)-u(5,:)-u(6,:)-u(7,:)
% hold on 
% plot(wt(1,:) + pv(1,:) - D(1,:))
