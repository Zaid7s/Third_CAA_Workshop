clc
clear all
% clf
% %%
% figure(2)
% datfiles  = dir('*Mean_Petu*');
% for k = 1 :  length(datfiles)
%     data = load(datfiles(k).name);
%     plot(data(:,1),data(:,4),'LineWidth',2.0)
%         xlabel('Domain')
%         grid on
%         grid minor
%         xlim([-10 10])
%         ylim([0.4 0.8])
%         ylabel('Mean Pressure')
%         pause(0.001)
% %         if (k == 1)
% %             print(['CAA_Mean', num2str(k)], '-djpeg', '-r300')
% %         end
% end
%%
figure(8)
datfiles  = dir('*nTPetu*');
for k = 1 :  length(datfiles)
    data = load(datfiles(k).name);
    plot(data(:,1),smooth(data(:,4)),'LineWidth',2.0)
        xlabel('Domain')
        grid on
        grid minor
        ax = gca;
        xlim([-10 10])
        set(gca,'XTick',-10:(1):10)
        ylim([-0.0001 0.0001])
        set(gca,'YTick',-0.0001:(0.00002):0.0001)
        ax.YAxis.Exponent = 0;
%         set(gca,'XTick',Begin:Increment:End);
        ylabel('Pressure Peturbation')
        title(['Time: ', num2str(data(1, 8))])
        pause(0.001)
%         if (k == 1)
%             print(['CAA_Mean', num2str(k)], '-djpeg', '-r300')
%         end
end
%%
figure(9)
datfiles  = dir('*nTPetu*');
for k = 1 :  length(datfiles)
    data = load(datfiles(k).name);
    [Max_Upper_Limit, Min_Lower_Limit] = envelope(smooth(data(:, 4)));
        plot(data(:,1),smooth(data(:, 4)),'LineWidth',2.0)
        hold on 
        plot(data(:,1),smooth(Max_Upper_Limit),'-.','LineWidth',2.0)
        hold on 
        plot(data(:,1),smooth(Min_Lower_Limit),'-.','LineWidth',2.0)
        hold on 
        xlabel('Domain')
        grid on
        grid minor
        xlim([-10 10])
        set(gca,'XTick',-10:(1):10)
        ax = gca;
        ylim([-0.0001 0.0001])
        set(gca,'YTick',-0.0001:(0.00002):0.0001)
        ax.YAxis.Exponent = 0;
        ylabel('Max Pressure Peturbation')
        title(['Time: ', num2str(data(1, 8))])
        hold off
        pause(0.001)
end
%%
figure(10)
datfiles  = dir('*nTPetu*');
for k = 1 :  length(datfiles)
    data = load(datfiles(k).name);
    [Max_Upper_Limit, Min_Lower_Limit] = envelope(data(:, 4));
    Mean_Limit = smooth(Max_Upper_Limit + Min_Lower_Limit)/2;
    Max_Disturbance = smooth(Max_Upper_Limit) - smooth(Mean_Limit);
        plot(data(:,1),smooth(Max_Disturbance),'LineWidth',2.0)
        xlabel('Domain')
        grid on
        grid minor
        xlim([-1.5 1.5])
        set(gca,'XTick',-1.5:(0.1):1.5)
        ax = gca;
        ylim([0 0.00011])
        set(gca,'YTick',0:(0.00001):0.00011)
        ax.YAxis.Exponent = 0;
        ylabel('Max Pressure Peturbation')
        title(['Time: ', num2str(data(1, 8))])
        hold off
        pause(0.001)
end
