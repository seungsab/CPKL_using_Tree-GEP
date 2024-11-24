function newinds=brood(pop,params,state,data,indlist)
%BROOD    Creates new individuals for GPLAB by brood recombination.
%   NEWINDS=BROOD(POPULATION,PARAMS,STATE,DATA,PARENTS) returns the two
%   best new individuals among the larger brood created by multiple
%   applications of crossover.
%
%   Input arguments:
%      POPULATION - the population where the parents are (array)
%      PARAMS - the parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      DATA - the current dataset for the algorithm to run (array)
%      PARENTS - the indices of the parents in POPULATION (1x2 matrix)
%   Output arguments:
%      NEWINDS - the two newly created individuals (1x2 matrix)
%
%   See also CROSSOVER, MUTATION, APPLYOPERATOR
%
%   Copyright (C) 2012-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

tmpinds=[];

for b=1:params.broodpairs % create params.broodpairs pairs of offpring
    inds=crossover(pop,params,state,data,indlist);
    if isempty(tmpinds)
        tmpinds=inds;
    else
        tmpinds(end+1:end+2)=inds;
    end
end

% calculate fitnesses
for i=1:length(tmpinds)
    tmpinds(i)=calcfitness(tmpinds(i),params,data,state,0);
end
    
% sort by fitness
%(don't forget to use adjustedfitness instead of fitness, here and elsewhere)
if params.lowerisbetter
    [ans,i]=sort([tmpinds.adjustedfitness],2,'ascend');
else
    [ans,i]=sort([tmpinds.adjustedfitness],2,'descend');
end

% choose the best 2 offspring (this brood crossover does not maintain avg
% pop size because the 2 offspring may not be from the same operation)
newinds=tmpinds(i(1:2));
