function [y] = Rosenbrock_fcnt(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Rosenbrock FUNCTION
%
% INPUTS:
%
% xx = [x1,...,xn]
% 
%         xn �� [-2, 3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = zeros(size(xx,1),1);

for ii=1:size(xx,2)-1
    y = y + (100*(xx(:,ii+1)-(xx(:,ii).^2)).^2+(1-xx(:,ii).^2));
end


end