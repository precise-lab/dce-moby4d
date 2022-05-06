% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     April 27, 2022
%

A = 0.6;    % mM
B = 0.18;   % min^(-1)
C = 0.45;   % mM
D = 0.5;    % min^(-1)
E = 0.013;  % min^(-1)

Ktrans  = 0.25;     % s^(-1)
kep     = 0.625;     % s^(-1)

t0 = 0;

syms ctoi(t);
ode = diff(ctoi) + kep*ctoi == ...
    Ktrans*(A*t*exp(-B*t) + C*(1 - exp(-D*t))*exp(-E*t));
cond = ctoi(t0) == 0;
ctoiSol(t) = dsolve(ode, cond)

t = ((1:1224) - 1)*0.5/60;
cp = A*t.*exp(-B*t) + C*(1 - exp(-D*t)).*exp(-E*t);
ctoi = (47601125*exp(-(5*t)/8))/30163168 ...
    + exp(-(5*t)/8).*((25*exp((153*t)/250))/136 ...
    - (225*exp((14*t)/125))/224 ...
    + (30*exp((89*t)/200).*(89*t - 200))/7921);

t = t*60; % min to s

figure; plot(t, cp, 'r', 'LineWidth', 3); hold on;
plot(t, ctoi, 'b--', 'LineWidth', 3); hold on;
ylim([0, 2.2]);
ylabel('Concentration (mM)');
xlabel('Time (s)');
legend({'C_p(t)', 'C_{TOI}(t)'});
set(gca, 'FontSize', 24);

save('contrast_agent_curve.mat');