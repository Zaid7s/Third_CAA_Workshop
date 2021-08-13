clc
clear all
% clf
%%
xMin    = -10;
xMax    = 10;
x0      = 2;
epsilon = 10^-5;
omega   = 0.6*pi;
M       = 0.4;
maxTime = 500;
iMax1    = 201;
iMax2    = 1201;
iMax3    = 201;


delta_x1 = ( xMax - xMin )/(iMax1 - 1);
delta_x2 = ( xMax - xMin )/(iMax2 - 1);
delta_x3 = ( xMax - xMin )/(iMax3 - 1);

iMax_Nx0 = ( -x0 - xMin )/delta_x1;
iMax_Bx0 = ( x0 - -x0  )/delta_x2;
iMax_Px0 = ( xMax - x0  )/delta_x3;

iMax = iMax_Nx0 + iMax_Bx0 + iMax_Px0 + 1

for i = 1 : iMax
    if (i <= iMax_Nx0 + 1)
        x(i) = xMin + (i - 1)*delta_x1;
        dX(i) = delta_x1;
    elseif (i >= iMax_Nx0 + 1) && (i <= iMax_Nx0 + iMax_Bx0 + 1)
        x(i) = x(i - 1) + delta_x2;  
        dX(i) = delta_x2;
    else
        x(i) = x(i - 1) + delta_x3;
        dX(i) = delta_x3;
    end
    ii(i) = i;
    Jac(i) = 1/dX(i);
    Zi(i) = x(i);
end
% figure(1)
% plot(x, dX, 'LineWidth', 2.0)
% figure(2)
% plot(ii, Jac, 'LineWidth', 2.0)
%%
for nT = 1 : maxTime
    time(nT)  = nT/10;
    for i = 1 : iMax
        rho(i)          =  epsilon*sin(omega*(x(i)/(1 - M) + time(nT)));
        u(i)            = -epsilon*sin(omega*(x(i)/(1 - M) + time(nT)));
        p(i)            =  epsilon*sin(omega*(x(i)/(1 - M) + time(nT)));
        rho_peturb(i)   =  epsilon*cos(omega*(xMax/(1 - M)));
        u_peturb(i)     = -epsilon*cos(omega*(xMax/(1 - M)));
        p_peturb(i)     =  epsilon*cos(omega*(xMax/(1 - M)));
        y(i)            = cos(x(i) + time(nT));
    end
    yT(nT) = y(end);
    FFTTrans = fft(yT);
    plot(x, y,'LineWidth', 2.0)
    hold off
    grid minor
    pause(0.1)
end
%%
plot(time, abs(FFTTrans),'LineWidth', 2.0)
grid minor
