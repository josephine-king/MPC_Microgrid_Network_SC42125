% Load data to plot different control horizons

load("n3_june.mat")
x_n3 = x;
u_n3 = u;

load("n6_june.mat")
x_n6 = x;
u_n6 = u;

load("n12_june.mat")
x_n12 = x;
u_n12 = u;

load("n24_june.mat")
x_n24 = x;
u_n24 = u;


figure(1)
subplot(4,1,1)
stairs(100.*x_n3(1,:)/mg1.cap)
hold on
stairs(100.*x_n3(2,:)/mg2.cap)
hold on
stairs(100.*x_n3(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, N=3")
ylabel("Percent Charged (%)")
fontsize(14,"points")
grid on
lgd1 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd1,10,'points')
ylim([59,61])

subplot(4,1,2)
stairs(100.*x_n6(1,:)/mg1.cap)
hold on
stairs(100.*x_n6(2,:)/mg2.cap)
hold on
stairs(100.*x_n6(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, N=6")
ylabel("Percent Charged (%)")
fontsize(14,"points")
grid on
lgd2 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd2,10,'points')
ylim([59,61])

subplot(4,1,3)
stairs(100.*x_n12(1,:)/mg1.cap)
hold on
stairs(100.*x_n12(2,:)/mg2.cap)
hold on
stairs(100.*x_n12(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, N=12")
ylabel("Percent Charged (%)")
fontsize(14,"points")
grid on
lgd3 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd3,10,'points')
ylim([59,61])

subplot(4,1,4)
stairs(100.*x_n24(1,:)/mg1.cap)
hold on
stairs(100.*x_n24(2,:)/mg2.cap)
hold on
stairs(100.*x_n24(3,:)/mg3.cap, 'Color',[0.4,0.7,0.3])
title("Battery Percent Charged, N=24")
ylabel("Percent Charged (%)")
fontsize(14,"points")
grid on
lgd4 = legend(["MG1", "MG2", "MG3"]);
fontsize(lgd4,10,'points')
ylim([59,61])
xlabel("Time step (hours)")



fontsize(lgd1,12,'points')
fontsize(lgd2,12,'points')
fontsize(lgd3,12,'points')
fontsize(lgd4,12,'points')

figure(2)
subplot(4,1,1)
stairs(u_n3(2,:))
hold on
stairs(u_n3(5,:))
hold on
stairs(u_n3(9,:))
hold on
stairs(u_n3(12,:))
hold on
stairs(u_n3(16,:))
hold on
stairs(u_n3(19,:))

subplot(4,1,2)
stairs(u_n6(2,:))
hold on
stairs(u_n6(5,:))
hold on
stairs(u_n6(9,:))
hold on
stairs(u_n6(12,:))
hold on
stairs(u_n6(16,:))
hold on
stairs(u_n6(19,:))

subplot(4,1,3)
stairs(u_n12(2,:))
hold on
stairs(u_n12(5,:))
hold on
stairs(u_n12(9,:))
hold on
stairs(u_n12(12,:))
hold on
stairs(u_n12(16,:))
hold on
stairs(u_n12(19,:))

subplot(4,1,4)
stairs(u_n24(2,:))
hold on
stairs(u_n24(5,:))
hold on
stairs(u_n24(9,:))
hold on
stairs(u_n24(12,:))
hold on
stairs(u_n24(16,:))
hold on
stairs(u_n24(19,:))

sum3 = sum([u_n3(2,:), u_n3(5,:), u_n3(9,:), u_n3(12,:), u_n3(16,:), u_n3(19,:)]);
sum6 = sum([u_n6(2,:), u_n6(5,:), u_n6(9,:), u_n6(12,:), u_n6(16,:), u_n6(19,:)]);
sum12 = sum([u_n12(2,:), u_n12(5,:), u_n12(9,:), u_n12(12,:), u_n12(16,:), u_n12(19,:)]);
sum24 = sum([u_n24(2,:), u_n24(5,:), u_n24(9,:), u_n24(12,:), u_n24(16,:), u_n24(19,:)]);

figure(3)
bar([sum3, sum6, sum12, sum24])