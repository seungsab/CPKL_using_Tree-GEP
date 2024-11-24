function [y] = Cantilever_hollow_tube(xx)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2, ..., x9]
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t   = xx(:,1); % t: Thickness of tube (mm) //  Normal(5.0,0.5) => [1.5 8.5] OK
d   = xx(:,2); % d: Outside diameter (mm) //  Normal(42.0,4.2) => [10 75] OK
L1  = xx(:,3); % L1:  Distance of F1 (mm) //  Normal(120.0,12.0) => [45 195] OK
L2  = xx(:,4); % L2:  Distance of F2 (mm) //  Normal(50.0,5.0) => [15 85] OK
F1  = xx(:,5); % F1: External Force #1 with the angle of 15 deg (kN) // Log-normal (3.0,0.1)
F2  = xx(:,6); % F2: External Force #2 with the angle of 25 deg (kN) // Log-normal (3.0,0.1)
P   = xx(:,7); % P: External Force #3 (axial force) (kN) // Log-normal (12.0,0.1)
T   = xx(:,8); % T: Tortional Force (kN)  // Log-normal (100.0,0.1)
Sy  = xx(:,9); % Sy: Material Yield Strength (MPa) // Log-normal (220.0,0.1)

% CONVERT DEG to RAD
beta1 = deg2rad(15); beta2 = deg2rad(25);

% COMPUTE SECTIONAL PROPERTIES
A = pi/4*(d.^2+(d-2*t).^2); % mm2
c = d/2; % mm
I = pi/64*(d.^4+(d-2*t).^4); %% mm4

% COMPUTE BENDING MOMENT
M = F1.*L1.*cos(beta1)+F2.*L2.*cos(beta2); % kN*mm

% COMPUTE NORMAL STRESS
sigma_x = (P+F1.*sin(beta1)+F2.*sin(beta2))./A+M.*c./I; %kN/mm2
% COMPUTE TORSION STRESS
tau_zx = T.*d./(4*I); % kN/mm2

% COMPUTE MAXIMUM von-MISES STRESS (kN/mm2)
sigma_max = sqrt(sigma_x.^2+2*tau_zx.^2);
sigma_max = sigma_max * 1e3; % kN/mm2 => N/mm2 (MPa)

% Performance function
y = Sy./(sigma_max); % MPa

end