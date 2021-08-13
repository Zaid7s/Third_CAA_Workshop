clc
clear all
syms epsi omega x M t bar_C bar_rho

rho_pet = epsi*cos(omega*(x/(0.6) + t));
u_pet   = -epsi*cos(omega*(x/(0.6) + t));
p_pet   = epsi*cos(omega*(x/(0.6) + t));

rho = diff(rho_pet, t)
u   = diff(u_pet, t)
p   = diff(p_pet, t)

A = p - bar_C*bar_rho*u

