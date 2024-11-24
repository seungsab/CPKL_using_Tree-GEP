function [Y_levy] = levy2D(X)

% ===============================================================
% Dimensions: 2
% ===============================================================
% Input:
% X = X(i)~Uniform[10,10]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref:
% (1) Laguna, M., Martí, R., 2005. Experimental testing of advanced scatter search designs for global optimization of multimodal functions. 
%     Journal of Global Optimization 33(2) 235-255.
% (2) Global Optimization Test Functions Index. Retrieved March 2018, 
%     from http://infinity77.net/global_optimization/test_functions.html#test-functions-index.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y_levy = ((sin(3 * pi * X(1)))^2) + ((X(1) - 1)^2 * (1 + (sin(3 * pi * X(2)))^2)) + ((X(2) - 1)^2 * (1 + (sin(2 * pi * X(2)))^2));
end
