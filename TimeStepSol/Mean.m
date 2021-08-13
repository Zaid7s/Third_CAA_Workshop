clc
clear all
clf
%%
figure(6)
datfiles  = dir('*Mean.dat');
for k = 1 : length(datfiles)
    data = load(datfiles(k).name);
    plot(data(:,1),data(:,4),'LineWidth',2.0)
    hold on
%     plot(data(:,1),data(:,3),'LineWidth',2.0)
%     hold on
%     plot(data(:,1),data(:,4),'LineWidth',2.0)
%     hold on
    plot(data(:,1),data(:,8),'LineWidth',2.0)
    hold off
        xlabel('Domain')
        grid on
        grid minor
        xlim([-1 1])
        xlim([-10 10])
%         ylim([0.4 0.8])
        ylabel('Mean pressure')
        pause(0.001)
end
%%
% figure(7)
% datfiles  = dir('*Mean.dat');
% for k = 1 : length(datfiles)
%     data = load(datfiles(k).name);
%     plot(data(:,1),(data(:,4) - data(:,8)),'LineWidth',2.0)
%     hold on
%     hold off
%         xlabel('Domain')
%         grid on
%         grid minor
%         xlim([-10 10])
% %         ylim([-0.0001 0.0001])
%         ylabel('p')
%         pause(0.001)
% end
