function ind=M3GPaccuracy(ind,params,data,terminals,varsvals)
%M3GPACCURACY    Measures the fitness of a GPLAB individual according to M3GP.
%   M3GPACCURACY(INDIVIDUAL,PARAMS,DATA,TERMINALS,VARSVALS) returns
%   the M3GP fitness of INDIVIDUAL, measured as the accuracy obtained on the
%   DATA dataset, and also returns the class obtained in each fitness case.
%
%   Input arguments:
%      INDIVIDUAL - the individual whose fitness is to measure (struct)
%      PARAMS - the current running parameters (struct)
%      DATA - the dataset on which to measure the fitness (struct)
%      TERMINALS - the variables to set with the input dataset (cell array)
%      VARSVALS - the string of the variables of the fitness cases (string)
%   Output arguments:
%      INDIVIDUAL - the individual whose fitness was measured (struct)
%
%   See also CALCFITNESS, ANTFITNESS
%
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% map the data on the ndimensional space with the hyperfeatures specified
% by the branches of the tree:

X=data.example;
nsamples=size(X,1);
nclasses=length(data.classes);
ndimensions=length(ind.tree.kids);
ndimensionalres=zeros(nsamples,ndimensions);

for d=1:ndimensions
    if ndimensions==1
        outstr=ind.str(6:end-1); % remove the initial "root(" and final ")" parts
    else
        outstr=tree2str(ind.tree.kids{d});
    end
    for i=params.numvars:-1:1
        %outstr=strrep(outstr,strcat('X',num2str(i)),strcat('X(:,',num2str(i),')'));
        outstr=strrep(outstr,['X',sprintf('%d',i)],['X(:,',sprintf('%d',i),')']);
    end
    try
        res=eval(outstr);
    catch % because of the "nesting 32" error of matlab
        res=str2num(evaluate_tree(ind.tree.kids{d},X));
    end
    if length(res)<nsamples
        res=res*ones(nsamples,1);
    end
    ndimensionalres(:,d)=res;
end

if isempty(ind.mapping)
    % training, save everything:
    ind.mapping=ndimensionalres;
    [ind.covariancematrix,ind.inversematrix,ind.centroids]=calcCovarianceCentroids(ind.mapping,data);
    samples2classify=ind.mapping;
else
    % otherwise, testing, do not save anything
    samples2classify=ndimensionalres;
end


% mapping done, covariance matrix and centroids calculated (and inverse when possible), now fitness:

%build matrix of distances between each sample and each centroid:
%nsamples2classify=size(samples2classify,1);
%nmappedsamples=size(ind.mapping,1);
distancematrix=zeros(nsamples,nclasses);

for c=1:nclasses
    repcentroid=repmat(ind.centroids(:,c)',nsamples,1);
    % mahalanobis distance (euclidean when inverse is identity):

    %'option1: sample by sample'
     for s=1:nsamples
     	distancematrix(s,c)=sqrt((samples2classify(s,:)-repcentroid(s,:))*ind.inversematrix{c}*(samples2classify(s,:)-repcentroid(s,:))');
     end
     
     %'option2: all samples at once and take only diagonal of result' ---- which one is faster???
     % this one may be faster, but runs out of memory with big datasets
      %distancematrix(:,c)=sqrt(diag((samples2classify-repcentroid)*ind.inversematrix{c}*(samples2classify-repcentroid)'));
end

% find predicted classes, based on distances to class centroids:
i=repmat(1:nclasses,nsamples,1); % matrix of indices: [1 2 3; 1 2 3; etc]
fi=(distancematrix==min(distancematrix')').*i; % matrix of zeros and
% indices where the minimum distance was found:
% [0 0 3; 1 0 0; 0 2 0; 0 2 3; etc] where last case means that a sample is
% at the same distance to two different centroids (of 2nd and 3rd classes)
fi(fi==0)=nclasses+1; % now the zeros became the maximum, nclasses+1
ind.result=data.classes(min(fi')); % class is given by the minimum index

ind.fitness=sum(ind.result==data.result)/nsamples; % fitness is accuracy
ind.adjustedfitness=ind.fitness;

% Maybe use this to calculate certainty, so save it for now:
%ind.distancematrix=distancematrix;

