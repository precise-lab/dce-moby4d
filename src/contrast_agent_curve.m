% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     April 27, 2022
%

%Implementation of the TOFTS model from:
% https://onlinelibrary.wiley.com/doi/10.1002/mrm.22861

%endTime = 17*36; %s
%dt      = 0.1;    %s

function [ca, c_perf] = contrast_agent_curve(endTime, dt);
A = 0.6;    % mM
B = 0.18/60.;   % s^(-1)
C = 0.45;   % mM
D = 0.5/60.;    % s^(-1)
E = 0.013/60;  % s^(-1)



Ktrans  = 0.25;     % s^(-1)
kep     = 0.625;     % s^(-1)


t = 0:dt:endTime;
% Input arterial function
ca = A*t.*exp(-B*t) + C*(1 - exp(-D*t)).*exp(-E*t);
% Perfusion
c_perf = Ktrans*conv( exp(-kep*t), ca);

figure; plot(t, cp, 'r', 'LineWidth', 3); hold on;
plot(t, c_perf, 'b--', 'LineWidth', 3); hold on;
ylim([0, 2.2]);
ylabel('Concentration (mM)');
xlabel('Time (s)');
legend({'C_a(t)', 'C_{perf}(t)'});
set(gca, 'FontSize', 24);

save('contrast_agent_curve.mat');

end
