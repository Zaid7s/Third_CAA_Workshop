clc
clf
clear all
format long
xMax    = 10;
epsilon = 10^-5;
omega   = 0.6*pi;
M       = 0.4;
iMax    = 201;
for i = 1 : iMax
    time(i)     = (i - 1)/60;
    rho(i)      = epsilon*cos(omega*(xMax/(1 + M) + time(i)));
    u(i)        = -epsilon*cos(omega*(xMax/(1 + M) + time(i)));
    p(i)        = epsilon*cos(omega*(xMax/(1 + M) + time(i)));
end

figure(1)
    plot(time, rho,'LineWidth', 2.0)
    hold on
    plot(time,u,'LineWidth', 2.0)
    hold on
    xlabel('time')
    ylabel('Peturbation')
    xlim([0 time(end)])
    ylim([-2e-5 2e-5 ])
    grid on
    grid minor