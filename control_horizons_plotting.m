load("n1_june.mat")
x_n1 = x;
u_n1 = u;

load("n12_june.mat")
x_n12 = x;
u_n12 = u;

load("n24_june.mat")
x_n24 = x;
u_n24 = u;


figure(4)
subplot(3,1,1)
stairs(100.*x_n1(1,:)/mg1.cap)
hold on
stairs(100.*x_n1(2,:)/mg2.cap)
hold on
stairs(100.*x_n1(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, N=1")
ylim([20 70])
ylabel("Battery Percent Charged (%)")
fontsize(14,"points")
grid on
lgd1 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd1,10,'points')
ylim([59,61])

subplot(3,1,2)
stairs(100.*x_n12(1,:)/mg1.cap)
hold on
stairs(100.*x_n12(2,:)/mg2.cap)
hold on
stairs(100.*x_n12(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, N=12")
ylim([20 70])
ylabel("Battery Percent Charged (%)")
fontsize(14,"points")
grid on
lgd2 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd2,10,'points')
ylim([59,61])

subplot(3,1,3)
stairs(100.*x_n24(1,:)/mg1.cap)
hold on
stairs(100.*x_n24(2,:)/mg2.cap)
hold on
stairs(100.*x_n24(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, N=24")
ylim([20 70])
ylabel("Battery Percent Charged (%)")
fontsize(14,"points")
grid on
lgd3 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd3,10,'points')
ylim([59,61])


fontsize(lgd1,12,'points')
fontsize(lgd2,12,'points')
fontsize(lgd3,12,'points')

figure(5)
subplot(3,1,1)
stairs(u_n1(2,:))
hold on
stairs(u_n1(9,:))
hold on
stairs(u_n1(16,:))
hold on
stairs(u_n1(5,:))
hold on
stairs(u_n1(12,:))
hold on
stairs(u_n1(19,:))
title("Power exchanged with DNO")
ylabel("Power (kW)")
fontsize(14,"points")
grid on
ylim([0 400])
lgd1 = legend(["MG1 to DNO", "MG2 to DNO", "MG3 to DNO", "DNO to MG1", "DNO to MG2", "DNO to MG3"]);

subplot(3,1,2)
stairs(u_n12(2,:))
hold on
stairs(u_n12(9,:))
hold on
stairs(u_n12(16,:))
hold on
stairs(u_n12(5,:))
hold on
stairs(u_n12(12,:))
hold on
stairs(u_n12(19,:))
title("Power exchanged with DNO")
ylabel("Power (kW)")
fontsize(14,"points")
grid on
ylim([0 400])
lgd2 = legend(["MG1 to DNO", "MG2 to DNO", "MG3 to DNO", "DNO to MG1", "DNO to MG2", "DNO to MG3"]);


subplot(3,1,3)
stairs(u_n24(2,:))
hold on
stairs(u_n24(9,:))
hold on
stairs(u_n24(16,:))
hold on
stairs(u_n24(5,:))
hold on
stairs(u_n24(12,:))
hold on
stairs(u_n24(19,:))
title("Power exchanged with DNO")
ylabel("Power (kW)")
fontsize(14,"points")
grid on
ylim([0 400])
lgd3 = legend(["MG1 to DNO", "MG2 to DNO", "MG3 to DNO", "DNO to MG1", "DNO to MG2", "DNO to MG3"]);


fontsize(lgd1,12,'points')
fontsize(lgd2,12,'points')
fontsize(lgd3,12,'points')