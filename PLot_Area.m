clc
clear all
filename1            = 'Area_RDRPSten_NPTS401.dat';
A1                   = importdata(filename1);
figure(1)
    plot(A1(:, 1), abs(A1(:, 4)),'LineWidth', 2.0)
    grid on
    grid minor
    xlabel('Grid Point')
    ylabel('Jacobian')
    ax = gca;
    xlim([1 size(A1, 1)])
    set(gca,'XTick',1:50:size(A1, 1))
    ax.XAxis.Exponent = 0;
    ylim([0 100])
    set(gca,'YTick',0:10:100)
    ax.YAxis.Exponent = 0;
%     print('i_vs_J_CAA', '-djpeg', '-r150')
figure(2)
    plot(A1(:, 1), abs(A1(:, 5)),'LineWidth', 2.0)
    grid on
    grid minor
    xlabel('Grid Point')
    ylabel('\Delta{x}')
    ax = gca;
    xlim([1 size(A1, 1)])
    set(gca,'XTick',1:50:size(A1, 1))
    ax.XAxis.Exponent = 0;
    ylim([0 0.12])
    set(gca,'YTick',0:0.01:0.12)
    ax.YAxis.Exponent = 0;
%     print('i_vs_dX_CAA', '-djpeg', '-r150')
figure(3)
    plot(A1(:, 2), abs(A1(:, 4)),'LineWidth', 2.0)
    grid on
    grid minor
    xlabel('Domain')
    ylabel('Jacobian')
    ylim([0 100])
    set(gca,'YTick',0:10:100)
    ax.YAxis.Exponent = 0;
    xlim([-10 10])
%     print('x_vs_J_CAA', '-djpeg', '-r150')
figure(4)
    plot(A1(:, 2), abs(A1(:, 5)),'LineWidth', 2.0)
    grid on
    grid minor
    xlabel('Domain')
    ylabel('\Delta{x}')
    ylim([0 0.12])
    set(gca,'YTick',0:0.01:0.12)
    ax.YAxis.Exponent = 0;
    xlim([-10 10])
%     print('x_vs_dX_CAA', '-djpeg', '-r150')
