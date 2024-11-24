function Y = uq_liquidHydrogenTank(X)
% UQ_LIQUIDHYDROGENTANK evaluates the limit state function of a liquid hydrogen tank model.
%
%   Y = UQ_LIQUIDHYDROGENTANK(X) evaluates the limit state function of a 
%   liquid hydrogen tank model for N-by-5 input matrix X, where N is 
%   the number of evaluation points; and returns a column vector of
%   length N. The five input variables are:
%     X(:,1): t_plate, thickness of plate (-)
%     X(:,2): t_h, thickness of honeycomb (-)
%     X(:,3): N_x, load on tank, x-component (-)
%     X(:,4): N_y, Load on tank, y-component (-)
%     X(:,5): N_xy, Load on tank, xy-component (-)
%
%   Note:
%
%   - This function is used as the test function in Example 5.2
%     of Bichon et al. (2011).
%
%   Reference:
%
%   - Bichon, B.J., J. M. McFarland, and S. Mahadevan. (2011).
%     Efficient surrogate models for system reliability analysis of systems
%     with multiple failure modes. Reliability Engineering and System 
%     Safety, vol. 96, pp. 1386-1395.
%
% lb = [0.04,0.03,-375,4500,-720]
% ub = [0.15,0.20,400,5100,-580]
% https://uqworld.org/t/liquid-hydrogen-tank-problem/58
%% Check input
%
narginchk(1,1)
assert(size(X,2)==5,'Only 5 input variables are allowed!')

%% Get input
%
t_plate = X(:,1); % Plate thickness [-]
t_h = X(:,2); % Honeycomb thickness [-]
N_x = X(:,3); % Load on tank (x-component) [-]
N_y = X(:,4); % Load on tank (y-component) [-]
N_xy = X(:,5); % Load on tank (xy-component) [-]

%% Intermediate computations
%
x_1 = 4*(t_plate - 0.075); 
x_2 = 20*(t_h - 0.1); 
x_3 = -6000*(1./N_xy + 0.003); 

%% Compute the three different performance functions
%
% Von Mises stress failure
P_vm = 84000*t_plate./sqrt(N_x.^2 + N_y.^2 - N_x.*N_y + 3*N_xy.^2) - 1; 

% Isotropic strength failure
P_is = 84000*t_plate./abs(N_y) - 1; 

% Honeycomb buckling failure
P_hb = 0.847 + 0.96*x_1 + 0.986*x_2 - 0.216*x_3 + 0.077*x_1.^2 +...
       0.11*x_2.^2 + 0.007*x_3.^2 + 0.378*x_1.*x_2 - 0.106*x_1.*x_3 +...
       0.11*x_2.*x_3; 

%% Evaluate the overall performance function
Y = min([P_vm P_is P_hb], [], 2);  
        
end