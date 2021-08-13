clc
clear all
filename1            = 'LUSGS_RDRPSten_D_10_Time.dat';
A1                   = importdata(filename1);
filename2            = 'LUSGS_RDRPSten_D_10_FFFT.dat';
A2                   = importdata(filename2);
% -------------------------------------------------------------------------------------------------- %
% (Number of Time Steps) || (Exit Density) || (Exit Velocity) || (Exit Pressure)|| (Time) || (Cycle) %
% -------------------------------------------------------------------------------------------------- %
figure(1)
plot(A1(1:255, 5), A1(1:255, 4), 'LineWidth', 2.0)
grid on
grid minor

%% Numerical Fourier
Numerical_Pressure_Fft = fft(A1(:, 4));
Numerical_Pressure_ABS = abs(Numerical_Pressure_Fft);

figure(2)
plot(Numerical_Pressure_ABS, 'LineWidth', 2.0)
grid on
grid minor
%% Exact Fourier
omega       = 0.6*pi;
finalTIme   = (2*pi/omega)*A1(end, 6);
L           = length(A1);
time        = [0: 1/L : finalTIme - 1/L];
Pressure    = (10^-5)*cos(omega*time + omega*(10/0.6));

Pressure_Fft = fft(Pressure);
Pressure_ABS = real(Pressure_Fft);
Pressure_Phase = angle(Pressure_Fft);

figure(3)
plot(time, Pressure, 'LineWidth', 2.0)
grid on
grid minor

figure(4)
plot(Pressure_ABS, 'LineWidth', 2.0)
grid on
grid minor

%% Numerical Fourier From Fortran
figure(5)
plot((A2(:, 3)), 'LineWidth', 2.0)
grid on
grid minor
% xlim([1 150])
