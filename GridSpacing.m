clc
clear all
clf
format long
x0      = 1;

xMin    = -10;
xMax    = 10;
xxx     = 10;
iMax_Left1   = 101;
iMax_Mid1    = 3101;
iMax_Righ1   = 101;

iMin        = 1;

dX_Left   = (xMax - xMin)/(iMax_Left1 - 1);
dX_Mid    = (xMax - xMin)/(iMax_Mid1 - 1);
dX_Right  = (xMax - xMin)/(iMax_Righ1 - 1);

iMin_Left   = 1;
iMax_Left   = (xMax - x0)/(dX_Left) + 1;

iMin_Mid   = iMax_Left;
iMax_Mid   = iMin_Mid + (x0 - -x0)/(dX_Mid);

iMin_Right   = iMax_Mid;
iMax_Right   = iMax_Mid + (-x0 - xMin)/(dX_Left);

iiMax = iMax_Right;
iiMin = iMin;


for i = iiMin : iiMax
    if (i >= iMin_Left) && (i <= iMax_Left)
        dX1(i) = dX_Right;
    elseif (i >= iMin_Mid) && (i <= iMax_Mid)
        dX1(i) = dX_Mid;
    elseif (i >= iMin_Right) && (i <= iMax_Right)
        dX1(i) = dX_Left;
    end
end

for i = iiMin : iiMax
    if (i == iMin)
        x(i) = xMin;
    else
        x(i) = x(i - 1) + dX1(i);
    end
    ii(i) = i;
end
for i = iiMin : iiMax
    if (i == iMin)
        x1(i) = xMin;
    else
        x1(i) = x1(i - 1) + dX1(i);
    end
    ii(i) = i;
end
for i = iiMin : iiMax
    if (x1(i) < -x0)
       dX(i) = dX_Left; 
    elseif (x1(i) >= -x0) && (x1(i) < 0)
       dX(i) = double(((tanh(-xxx*(x1(i) + (x0 - 0.5))))));
    elseif (x1(i) >= 0) && (x1(i) <= x0)
       dX(i) = double(((tanh(-xxx*(-x1(i) + (x0 - 0.5))))));
    elseif (x1(i) > x0) && (x1(i) <= xMax)
       dX(i) = dX_Right;
    else
       dX(i) = 0;
    end
    ii(i) = i;
    Jac(i) = 1/dX(i);
end

MMin = min(dX)

for i = iiMin : iiMax
    if (x1(i) < -x0)
       dX(i) = dX_Left; 
    elseif (x1(i) >= -x0) && (x1(i) < 0)
       dX(i) = double(((tanh(-xxx*(x1(i) + (x0 - 0.5))))) - MMin);
    elseif (x1(i) >= 0) && (x1(i) <= x0)
       dX(i) = double(((tanh(-xxx*(-x1(i) + (x0 - 0.5))))) - MMin);
    elseif (x1(i) > x0) && (x1(i) <= xMax)
       dX(i) = dX_Right;
    else
       dX(i) = 0;
    end
    ii(i) = i;
    Jac(i) = 1/dX(i);
end

MMax = max(dX)

%% Get A Smooth Non-Uniform Gird
for nn = 1 : 1000001
    nFac    = (nn - 1)/10000000;
    MaxDX   = MMax;
    dXX     = (nn - 1)/1000000;
    Fac     = double((dXX - nFac)/MaxDX);
    for i = iiMin : iiMax
        if (x(i) < -x0)
           dX(i) = dXX; 
        elseif (x(i) >= -x0) && (x(i) < 0)
           dX(i) = double(((tanh(-xxx*(x(i) + (x0 - 0.5)))-MMin))*Fac + nFac);
        elseif (x(i) >= 0) && (x(i) <= x0)
           dX(i) = double(((tanh(-xxx*(-x(i) + (x0 - 0.5)))-MMin))*Fac + nFac);
        elseif (x(i) > x0) && (x(i) <= xMax)
           dX(i) = dXX;
        else
           dX(i) = 0;
        end
        ii(i) = i;
        Jac(i) = 1/dX(i);
    end
    for i = iiMin : iiMax
        if (i == iMin)
            x1(i) = xMin;
        else
            x1(i) = x1(i - 1) + dX(i);
        end
    ii(i) = i;
    end
    if (abs(x1(end) - 10) <= 10^-3)
        nFac
        Fac
        nn
        dXX
        round(20/nFac)
        round(20/dXX)
        break
    end
end

%% 
figure(5)
plot(x1, dX, 'LineWidth', 2.0)
% hold on 
% plot(x, dX1, 'LineWidth', 2.0)
grid on
grid minor
% xlim([-2 2])
ylim([0 0.12])
hold off
% 


