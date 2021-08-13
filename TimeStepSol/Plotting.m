clc
clear all
clf
%% Instantaneous and Envelope
figure(1)
datfiles  = dir('*nTPetu*');
for k = 1: 3: length(datfiles)
    if (k == 1)
        clf
    end
    data = load(datfiles(k).name);
    subplot(2, 1, 1)
        plot(data(:,1),smooth(data(:,2)),'LineWidth',2.0)
        hold off
        xlabel('Domain')
        grid on
        grid minor
        ax = gca;
        xlim([-10 10])
        set(gca,'XTick',-10:(1):10)
        ylim([-0.0001 0.0001])
        set(gca,'YTick',-0.0001:(0.000025):0.0001)
        ax.YAxis.Exponent = 0;
        ylabel('Pressure Peturbation')
    subplot(2, 1, 2)
        plot(data(:,1),smooth(data(:, 4)),'Color',[0, 0.4470, 0.7410],'LineWidth',0.75)
        hold on 
        xlabel('Domain')
        grid on
        grid minor
        xlim([-10 10])
        set(gca,'XTick',-10:(1):10)
        ax = gca;
        ylim([-0.0001 0.0001])
        set(gca,'YTick',-0.0001:(0.000025):0.0001)
        ax.YAxis.Exponent = 0;
        ylabel('Pressure Peturbation')
        pause(0.00001)
end

%% Calculate Mean
Mean_Limit                  = smooth(data(:, 14));
Mean_Limit_Exact            = smooth(data(:, 11) + data(:, 13))/2;
%% PLot Max Petrubation
Max_Disturbance         = smooth(data(:, 10)) - smooth(Mean_Limit);
Max_Disturbance_Exact   = smooth(data(:, 11)) - smooth(Mean_Limit_Exact);
figure(2)
subplot(1, 2, 2)
    plot(data(:,1),smooth(Max_Disturbance),'LineWidth',2.0)
    xlabel('Domain')
    grid on
    grid minor
    ax = gca;
%     ylim([0 0.00011])
%     set(gca,'YTick',0:(0.00001):0.00011)
%     ax.YAxis.Exponent = 0;
    xlim([-1.5 1.5])
    ylabel('Max Pressure Peturbation')
subplot(1, 2, 1)
    plot(data(:,1),smooth(Max_Disturbance),'LineWidth',2.0)
    xlabel('Domain')
    grid on
    grid minor
    ax = gca;
    ylim([5e-6 12e-6])
    xlim([-10 10])
%     ax.YAxis.Exponent = 2;
    ylabel('Max Pressure Peturbation')
%% PLot Mean
% figure(3)
%     plot(data(:,1),smooth(Mean_Limit),'LineWidth',2.0)
%     xlabel('Domain')
%     grid on
%     grid minor
%     ax = gca;
%     ylim([-0.0001 0.0001])
%     set(gca,'YTick',-0.0001:(0.000025):0.0001)
%     ax.YAxis.Exponent = 0;
%     xlim([-1 1])
%     ylabel('Max Pressure Peturbation')
% %     Mean_Limit(1)
