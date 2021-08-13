clc
clear all
clf

m = 16;
omega = 4*pi/m;
t = [0:1:m-1];
xCos = cos(omega*t);
xSin = sin(omega*t);

XCos = fft(xCos);
XSin = fft(xSin);

X_Mag_Cos = abs(XCos);
X_Mag_Sin = abs(XSin);
%%
figure(1)
plot(t, xCos, 'LineWidth', 2.0)
hold on
plot(t, xSin, 'LineWidth', 2.0)
grid on
hold off
grid minor
%%
figure(3)
plot(X_Mag_Cos, 'LineWidth', 2.0)
hold on
plot(X_Mag_Sin, 'LineWidth', 2.0)
grid on
hold off
grid minor
