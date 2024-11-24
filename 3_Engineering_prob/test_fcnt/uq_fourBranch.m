function Y = uq_fourBranch(X,P)
% UQ_FOURBRANCH computes the Four-branch limit state function.
%
%   Y = UQ_FOURBRANCH(X) evaluates the Four-branch limit state function for
%   N-by-2 input matrix X, where N is the number of evaluation points; and
%   returns a column vector of length N.  The default value of the scalar 
%   parameter P is 6.
%
%   Y = UQ_FOURBRANCH(X,P) evaluates the Four-branch limit state function
%   with a user-specified scalar parameter.
%
%   Reference:
%
%   - Schueremans, L. and D. van Gemert. (2005).
%     Benefit of splines and neural networks in simulation based
%     structural reliability analysis. Structural Safety, vol. 27, no. 3.
%     DOI:10.1016/j.strusafe.2004.11.001
%   - Echard, B., N. Gayton, and M. Lemaire. (2011).
%     AK-MCS: An active learning reliability method combining kriging and
%     Monte Carlo simulation. Reliability Engineering and System Safety,
%     vol. 33, no. 2, pp. 145-154. DOI:10.1016/j.ress.2012.10.008
%
% lb = [-7 -7];
% ub = [7 7];
% https://uqworld.org/t/liquid-hydrogen-tank-problem/58
%% Check inputs
%
narginchk(1,2)
assert(size(X,2)==2,'only 2 input variables allowed')

%% Evaluate function
%
if nargin == 1 % Default scalar parameter
    P = 6;
    
    Y = min([3 + 0.1*(X(:,1)-X(:,2)).^2 - 1/sqrt(2)*(X(:,1)+X(:,2)),...
        3 + 0.1*(X(:,1)-X(:,2)).^2 + 1/sqrt(2)*(X(:,1)+X(:,2)),...
        (X(:,1)-X(:,2)) + 1/sqrt(2)*P,...
        (X(:,2)-X(:,1)) + 1/sqrt(2)*P],...
        [], 2);
end

if nargin == 2 % User-specified scalar parameter
    Y = min([3 + 0.1*(X(:,1)-X(:,2)).^2 - 1/sqrt(2)*(X(:,1)+X(:,2)),...
        3 + 0.1*(X(:,1)-X(:,2)).^2 + 1/sqrt(2)*(X(:,1)+X(:,2)),...
        (X(:,1)-X(:,2)) + 1/sqrt(2)*P,...
        (X(:,2)-X(:,1)) + 1/sqrt(2)*P],...
        [], 2);
end

end
