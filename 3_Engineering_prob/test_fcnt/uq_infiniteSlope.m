function Y = uq_infiniteSlope(X)
% UQ_INFINITESLOPE computes the Infinite Slope model limit state function.
%
%   Y = UQ_INFINITESLOPE(X) evaluates the 6-dimensional limit state
%   function of the Infinite Slope model for N-by-6 input matrix X, where N
%   is the number of evaluation points; and returns a column vector
%   of length N. The six input variables are:
%      X(:,1): H, depth of soil above bedrock (m)
%      X(:,2): U, relative height of water table (-)
%      X(:,3): phi, effective stress friction angle (-)
%      X(:,4): theta, slope inclination (-)
%      X(:,5): Gs, specific gravity of soil (-)
%      X(:,6): e, void ratio of soil (-)
%
%   Reference:
%
%   - Phoon, K.-K. (2008). Numerical recipes for reliability
%     analysis - a primer. In Reliability-based Design in Geotechnical
%     Engineering: Computations and Applications, K.-K Phoon, Ed.
%     London: CRC Press, pp. 34-35. DOI:10.1201/9781482265811
%
% lb = [2 0 0.4 0.25 2.5 0.3];
% ub = [8 1 0.8 0.5 2.7 0.6];
% https://uqworld.org/t/infinite-slope-model/102

%% Check input
%
narginchk(1,1)
assert(size(X,2) == 6, 'Exactly 6 input variables are needed!')

%% Define constants
%
k = 0.2;        % Degree of saturation of moist soil (-)
gamma_w = 9.81; % unit weight of water (kN/m3)

%% Assign input to local variables
%
H = X(:,1);     % Depth of soil above bedrock (m)
u = X(:,2);     % Relative height of water table(-)
phi = X(:,3);   % Effective stress friction angle (-)
theta = X(:,4); % Slope inclination (-)
Gs = X(:,5);    % Specific gravity of soil (-)
e = X(:,6);     % Void ratio of soil (-)

%% Calculate intermediate quantities
%
h = H.*u;  % Water table height (m)
gamma = gamma_w * (Gs + k*e)./(1 + e);   % Moist unit weight of soil (kN/m3)
gamma_sat = gamma_w * (Gs + e)./(1 + e); % Saturated unit weight of soil (kN/m3)

%% Evaluate the limit state function
%
T = (gamma.*(H-h) + h.*(gamma_sat - gamma_w)).*cos(theta).*tan(phi);
R = (gamma.*(H-h) + h.*gamma_sat).*sin(theta);

Y = T./R - 1;

end