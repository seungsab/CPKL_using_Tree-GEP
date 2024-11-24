function newind=pointmutation(pop,params,state,data,i)
%POINTMUTATION    Creates a new individual for GPLAB by point mutation.
%   NEWIND=POINTMUTATION(POPULATION,PARAMS,STATE,DATA,PARENT) returns a
%   new individual created by substituting 1 or n nodes of PARENT by new
%   random nodes with the same arity.
%
%   Input arguments:
%      POPULATION - the population where the parent is (array)
%      PARAMS - the parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      DATA - the current dataset for the algorithm to run (array)
%      PARENT - the index of the parent in POPULATION (integer)
%   Output arguments:
%      NEWIND - the newly created individual (struct)
%
%   See also CROSSOVER, APPLYOPERATOR
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

ind=pop(i);

% calculate number of nodes
if isempty(ind.nodes)
   ind.nodes=nodes(ind.tree);
end

if strcmp(params.pointmutationtype,'1point')
    x=intrand(1,ind.nodes); % 1 random mutation point
else
    rr=rand(1,ind.nodes); % one mutation probability per node
    x=find(rr<1/ind.nodes); % possibly several mutation points
end

% change all nodes in x to other functions of the same arity:
ind.tree=newnodefunc(ind.tree,x,state);

ind.id=[];
ind.origin=[params.pointmutationtype 'mutation'];
ind.parents=[pop(i).id];
ind.xsites=x;
ind.str=tree2str(ind.tree);
ind.fitness=[];
ind.adjustedfitness=[];
ind.result=[];
ind.testfitness=[];
ind.testadjustedfitness=[];
ind.nodes=ind.tree.nodes;
ind.introns=[];
ind.level=[];

newind=ind;
