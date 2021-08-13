clc
clear all
syms Q1 Q2 Q3 gamma rho u Etot p c J Area Source
%%
% E_tot        = p\(gamma - 1) + rho*u^2/2;
% Q1          = rho
% Q2          = rho*u
% Q3          = E_tot
% p           = (E_tot - (rho*u^2)/2)*(gamma - 1);
Q_n(1, 1)   = Q1;
Q_n(2, 1)   = Q2;
Q_n(3, 1)   = Q3;
F_1         = Q2;
F_2         = Q2*Q2/Q1 + (gamma - 1)*(Q3-(Q2*Q2/(Q1*2)));
F_3         = Q2/Q1*(Q3 + (gamma - 1)*(Q3-(Q2*Q2/(Q1*2))));
F(:, 1)     = [F_1;F_2;F_3];
dEdQ        = simplify(sym(jacobian([F],[Q_n])))
A1          = simplify(subs(subs(subs(dEdQ, Q1, rho/J), Q2, (rho*u)/J), Q3, (Etot)/J));

EigenValue_A    = eig(A1);
%%
S_1             = Source*Q2/Area;
S_2             = Source*(Q2*Q2/Q1)/Area;
S_3             = Source*F_3/Area;
S(:, 1)         = [S_1; S_2; S_3];
dSdQ            = simplify(sym(jacobian([S],[Q_n])))
Z1               = simplify(subs(subs(subs(dSdQ, Q1, (rho*u)), Q2, (rho*u*u)), Q3, (u*(Etot + p))));

% EigenValue_S    = eig(Z);
% % % dSdQ1           = subs(dSdQ, Q1, 0.98027);
% % % dSdQ2           = subs(dSdQ1, Q2, 0.19605);
% % % dSdQ3           = simplify(subs(dSdQ2, Q3, 1.7821));
% % % dSdQ4           = simplify(subs(dSdQ3, gamma, 1.4));
% % % dSdQ4           = vpa(simplify(subs(dSdQ4, E_tot, [p\(gamma - 1) + rho*u^2/2])), 10)
