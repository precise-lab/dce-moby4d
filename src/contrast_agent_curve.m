% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     April 27, 2022
%

%Implementation of the TOFTS model from:
% https://onlinelibrary.wiley.com/doi/10.1002/mrm.22861

%endTime = 17*36; %s
%dt      = 0.1;    %s

function [ca, c_perf, t] = contrast_agent_curve(startTime, endTime, dt)
% Injection time is at t=0.
A = 0.6/60;    % mM s^(-1)
B = 0.18/60.;   % s^(-1)
C = 0.45;   % mM
D = 0.5/60.;    % s^(-1)
E = 0.013/60;  % s^(-1)



Ktrans  = 0.25;     % s^(-1)
kep     = 0.625;     % s^(-1)


t = 0:dt:endTime;
% Input arterial function
ca_tp = A*t.*exp(-B*t) + C*(1 - exp(-D*t)).*exp(-E*t);
% Perfusion
c_perf_tp = Ktrans*conv( exp(-kep*t), ca_tp)*dt;
c_perf_tp = c_perf_tp(1:numel(t));

t = [startTime:dt:endTime];
ca = zeros(size(t));
ca(t>=0.) = ca_tp;
c_perf = zeros(size(t));
c_perf(t>=0.) = c_perf_tp;

end
