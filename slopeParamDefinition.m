% Parameter analysis (Error calculation)

%% Determination of us_star parameter
clear, clc

R = 8.3144626 * 1e-3;  % (KJ/mol K)
T1 = 273.15;           % (K) | 0 °C
T2 = 298.15;           % (K) | 25 °C
T3 = 323.15;           % (K) | 50 °C
Vmax = 307;

% We define: beta = 1 / us_star for simplicity
beta = 0:0.01:1;

dG = 0:0.1:20;
dG = -dG;
syms f1(x,y)
f1(x,y) = piecewise(x>=-1/y, -Vmax* y * x, x<-1/y, Vmax);
f2_T1 = Vmax*(1 - exp(dG./(R * T1)));
f2_T2 = Vmax*(1 - exp(dG./(R * T2)));
f2_T3 = Vmax*(1 - exp(dG./(R * T3)));
X = repmat(dG,size(beta,2),1);
Y = repmat(beta',1,size(dG,2));
F1 = f1(X,Y);
F3_T1 = abs(F1 - f2_T1);
D_T1  = sum(F3_T1,2);
F3_T2 = abs(F1 - f2_T2);
D_T2  = sum(F3_T2,2);
F3_T3 = abs(F1 - f2_T3);
D_T3  = sum(F3_T3,2);
plot(beta,D_T1)
hold on
plot(beta,D_T2)
plot(beta,D_T3)
hold off
legend('0°C','25°C','50°C','Location','Best')

% Pseudo-optimum
[M,P] = min(D_T2);
b_opt = beta(P);
us_star_opt = 1/b_opt;
us_star = 1./beta;
save('us_starResults.mat','D_T1','D_T2','D_T3','beta','us_star')

figure
plot(us_star,D_T2)
xlim([0.1 10])

% Error calculation

g = 4.3478; % Parameter (opt)

real_curve = Vmax*(1 - exp(dG./(R * T2))); % real curve
dG1 = dG(dG>=-g);
dG2 = dG(dG<=-g);
appr1 = -Vmax/g .* dG1; % approximation (line)
appr2 = Vmax * ones(size(dG2)); % approximation (constant)
appr = [appr1 appr2];

mean_error = mean(abs(real_curve-appr)/real_curve * 100);



plot(dG,real_curve)
hold on
plot(dG,appr)
hold off


