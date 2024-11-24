function D = mypdist(X,Y,M,C)
%MYPDIST    Calculates different types of distances for GPLAB M3GP method
%   MYPDIST(X,Y,M,C) returns one of 4 different types of distance,
%   based on the MATLAB functions PDIST and PDIST2:
%   1. pdist(X)
%   2. pdist2(X,Y) where Y X is Observation matrix and Y is a vector of mean values
%   3. pdist(X,'mahalanobis',M)
%   4. pdist2(X,Y,'mahalanobis',M)
%
%   References:
%      - Ingalalli, Silva, Castelli, Vanneschi (2014). A Multi-dimensional
%      Genetic Programming Approach for Multi-class Classification
%      Problems, EuroGP-2014.
%      - Muñoz, Silva, Trujillo (2015). M3GP – Multiclass Classification
%      with GP, EuroGP-2015.
%
%   See also PDIST, PDIST2 (MATLAB functions)
%
%   Created by Vijay Ingalalli
%   Copyright (C) 2014-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox


[n,~] = size(X);

if nargin == 1 %Euclidean
	length = n*(n-1)/2;
	D = zeros(length,1);
	k=1;
	for i=1:n-1
		for j=i+1:n
			D(k) = sqrt((X(i,:)-X(j,:))*(X(i,:)-X(j,:))');
			k = k+1;
		end
	end

elseif nargin == 2 %Euclidean between a matrix and a vector
	D = zeros(n,1);
	for i = 1:n
		D(i) = sqrt((X(i,:)-Y)*(X(i,:)-Y)');
	end

elseif nargin == 3	%Mahalanobis
	length = n*(n-1)/2;
	D = zeros(length,1);
	k=1;
	for i=1:n-1
		for j=i+1:n
			D(k) = sqrt((X(i,:)-X(j,:))*inv(M)*(X(i,:)-X(j,:))');
			k = k+1;
		end
	end

elseif nargin == 4 %Mahalanobis	between a matrix and a vector
	D = zeros(n,1);
	for i = 1:n
		D(i) = sqrt((X(i,:)-Y)*inv(C)*(X(i,:)-Y)');
	end
end	
