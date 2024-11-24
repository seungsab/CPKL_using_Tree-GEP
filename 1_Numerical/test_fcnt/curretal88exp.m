function [y] = curretal88exp(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CURRIN ET AL. (1988) EXPONENTIAL FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
% Currin, C., Mitchell, T., Morris, M., & Ylvisaker, D. (1988). A Bayesian approach to the design and analysis of computer experiments. Technical Report 6498. Oak Ridge National Laboratory.
% 
% Currin, C., Mitchell, T., Morris, M., & Ylvisaker, D. (1991). Bayesian prediction of deterministic functions, with applications to the design and analysis of computer experiments. Journal of the American Statistical Association, 86(416), 953-963.
% 
% Xiong, S., Qian, P. Z., & Wu, C. J. (2013). Sequential design and analysis of high-accuracy and low-accuracy computer codes. Technometrics, 55(1), 37-46.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2]
%
% xi ¡ô [0, 1], for all i = 1, 2. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = x(:,1);
x2 = x(:,2);

fact1 = 1 - exp(-1./(2*x2));
fact2 = 2300*x1.^3 + 1900*x1.^2 + 2092*x1 + 60;
fact3 = 100*x1.^3 + 500*x1.^2 + 4*x1 + 20;

y = fact1.*fact2./fact3;

end
