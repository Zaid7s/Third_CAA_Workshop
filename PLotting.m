clc
clear all
clf
figure(1)
datfiles  = dir('*.dat');

data = load(datfiles.name);
plot(data(:,1),data(:,4),'LineWidth',2.0)
    xlabel('Domain')
    grid on
    grid minor
    xlim([-10 10])
    ylabel('p')
    pause(0.001)
