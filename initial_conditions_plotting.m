
%load("ic_20_june.mat")
load("n24_june.mat")
save("ic_20_june.mat", "x", "u")
x_20 = x;
u_20 = u;

load("ic_40_june.mat")
x_40 = x;
u_40 = u;

load("ic_60_june.mat")
x_60 = x;
u_60 = u;

load("ic_80_june.mat")
x_80 = x;
u_80 = u;


figure(4)
subplot(4,1,1)
stairs(100.*x_20(1,:)/mg1.cap)
hold on
stairs(100.*x_20(2,:)/mg2.cap)
hold on
stairs(100.*x_20(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, Initial Charge = 20%")
ylim([20 80])
ylabel("Percent Charged (%)")
fontsize(14,"points")
grid on
lgd1 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd1,10,'points')

subplot(4,1,2)
stairs(100.*x_40(1,:)/mg1.cap)
hold on
stairs(100.*x_40(2,:)/mg2.cap)
hold on
stairs(100.*x_40(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, Initial Charge = 40%")
ylim([20 80])
ylabel("Percent Charged (%)")
fontsize(14,"points")
grid on
lgd2 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd2,10,'points')

subplot(4,1,3)
stairs(100.*x_60(1,:)/mg1.cap)
hold on
stairs(100.*x_60(2,:)/mg2.cap)
hold on
stairs(100.*x_60(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, Initial Charge = 60%")
ylim([20 80])
ylabel("Percent Charged (%)")
fontsize(14,"points")
grid on
lgd3 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd3,10,'points')

subplot(4,1,4)
stairs(100.*x_80(1,:)/mg1.cap)
hold on
stairs(100.*x_80(2,:)/mg2.cap)
hold on
stairs(100.*x_80(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, Initial Charge = 80%")
ylim([20 80])
ylabel("Percent Charged (%)")
fontsize(14,"points")
grid on
lgd4 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd4,10,'points')

fontsize(lgd1,12,'points')
fontsize(lgd2,12,'points')
fontsize(lgd3,12,'points')
fontsize(lgd4,12,'points')
