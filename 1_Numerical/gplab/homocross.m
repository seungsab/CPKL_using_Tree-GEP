function newinds=homocross(pop,params,state,data,indlist)
%HOMOCROSS    Creates new GPLAB individuals by homologous crossover.
%   NEWINDS=HOMOCROSS(POPULATION,PARAMS,STATE,DATA,PARENTS) returns two
%   new individuals created by swapping node functions or entire subtrees
%   of the two PARENTS at random homologous points.
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
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

i1=indlist(1);
i2=indlist(2);

ind1=pop(i1);
ind2=pop(i2);

% weights W to decide whether the first or second parent contributes
% with genetic material (see Moraglio's refs)

if strncmp(params.homocrosstype,'gde',3)
    W1=params.gdeW1;
else
    W1=0.5; % pure homologous crossover
end
%W2=1-W1;

if isempty(ind1.nodes)
   ind1.nodes=nodes(ind1.tree);
end
if isempty(ind2.nodes)
   ind2.nodes=nodes(ind2.tree);
end

% Identify possible cross points (ids and arities) in each parent:
[cp_ids1,cp_ids2,cp_arities1,cp_arities2]=crosspoints(ind1.tree,ind2.tree,[],[],[],[],{},{},[],[],[]);
% (start the recursive function with empty lists)

% Store original cross point ids:
ocp_ids1=cp_ids1;
ocp_ids2=cp_ids2;

% For each point decide whether to swap or not:
% (flip as many coins as possible cross points)
ncp=length(cp_ids1);
r=rand(1,ncp);
for i=1:ncp
    if r(i)<W1
        % swap:
        if cp_arities1(i)==cp_arities2(i)
            [ind1.tree,ind2.tree]=swapnodefunc(ind1.tree,ind2.tree,cp_ids1(i),cp_ids2(i),1);
            % (1 means both trees have the same shape until this point)
        else
            [nind1tree,nind2tree]=swapnodes(ind1.tree,ind2.tree,cp_ids1(i),cp_ids2(i));
            % update the remaining ids to swap, because of the size change
            d=nind1tree.nodes-ind1.tree.nodes;
            cp_ids1(i+1:end)=cp_ids1(i+1:end)+d;
            cp_ids2(i+1:end)=cp_ids2(i+1:end)-d;
            ind1.tree=nind1tree;
            ind2.tree=nind2tree;
        end
    end
end


ind1.str=tree2str(ind1.tree);
ind1.id=[];
ind1.origin=params.homocrosstype;
ind1.parents=[pop(i1).id,pop(i2).id];
ind1.xsites=[ocp_ids1',ocp_ids2'];
ind1.fitness=[];
ind1.adjustedfitness=[];
ind1.result=[];
ind1.testfitness=[];
ind1.testadjustedfitness=[];
ind1.level=[];
ind1.nodes=ind1.tree.nodes;
ind1.introns=[];

ind1.dimensions=[];
ind1.pruned=0;
ind1.mapping=[];
ind1.covariancematrix={};
ind1.inversematrix={};
ind1.centroids=[];
   
if strcmp(params.homocrosstype,'homocross + pointmutation')
    tmppop(1)=ind2;
    ind2=pointmutation(tmppop,params,state,data,1);
    ind2.origin=['homocross + ' params.pointmutationtype 'mutation'];
    ind2.xsites(end+1:length(ocp_ids2))=0;
    ocp_ids1(end+1:length(ind2.xsites))=0;
    ocp_ids2(end+1:length(ind2.xsites))=0;
    ind2.xsites=[ocp_ids2',ocp_ids1',ind2.xsites'];
else
    ind2.origin=params.homocrosstype;
    ind2.xsites=[ocp_ids2',ocp_ids1'];
end

ind2.str=tree2str(ind2.tree);      	
ind2.id=[];
ind2.parents=[pop(i2).id,pop(i1).id];
ind2.fitness=[];
ind2.adjustedfitness=[];
ind2.result=[];
ind2.testfitness=[];
ind2.testadjustedfitness=[];
ind2.level=[];
ind2.nodes=ind2.tree.nodes;
ind2.introns=[];   

ind2.dimensions=[];
ind2.pruned=0;
ind2.mapping=[];
ind2.covariancematrix={};
ind2.inversematrix={};
ind2.centroids=[];

newind1=ind1;
newind2=ind2;

%newinds=[newind1,newind2];
newinds=newind2; % ignore first individual
% (don't change this when you're using GDE!!!)
% (if you need 2 individuals, GDE should always get the second!)
% (if you change this to 2 individuals, make sure to mutate the first too)
