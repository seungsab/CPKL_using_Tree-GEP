function newinds=crossover(pop,params,state,data,indlist)
%CROSSOVER    Creates new individuals for the GPLAB algorithm by crossover.
%   NEWINDS=CROSSOVER(POPULATION,PARAMS,STATE,DATA,PARENTS) returns two new
%   individuals created by swaping subtrees of the two PARENTS at random
%   points.
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
%   See also MUTATION, APPLYOPERATOR
%
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

i1=indlist(1);
i2=indlist(2);

ind1=pop(i1);
ind2=pop(i2);

if isempty(ind1.nodes)
   ind1.nodes=nodes(ind1.tree);
end
if isempty(ind2.nodes)
   ind2.nodes=nodes(ind2.tree);
end

if params.M3GP
    % decide which crossover:
    r=rand;
    if r<0.5 % normal crossover, just do not touch the root node
        x1=intrand(2,ind1.nodes); % crosspoint1
        x2=intrand(2,ind2.nodes); % crosspoint2
        indorigin='M3GPcrossover_normal';
    else % choose a node from the list of dimension-root nodes
        nkids1=length(ind1.tree.kids);
        listnodes1=zeros(1,nkids1);
        for k=1:nkids1
            listnodes1(k)=ind1.tree.kids{k}.nodeid;
        end
        nkids2=length(ind2.tree.kids);
        listnodes2=zeros(1,nkids2);
        for k=1:nkids2
            listnodes2(k)=ind2.tree.kids{k}.nodeid;
        end
        x1=listnodes1(intrand(1,nkids1)); % crosspoint1, a dimension root
        x2=listnodes2(intrand(1,nkids2)); % crosspoint2, a dimension root
        indorigin='M3GPcrossover_swapdimensions';
    end
else
    x1=intrand(1,ind1.nodes); % crosspoint1
    x2=intrand(1,ind2.nodes); % crosspoint2
    indorigin='crossover';
end

% swap nodes in only one step:
[ind1.tree,ind2.tree]=swapnodes(ind1.tree,ind2.tree,x1,x2);

%xnode1=findnode(ind1.tree,x1); % node to swap1
%xnode2=findnode(ind2.tree,x2); % node to swap2
%ind1.tree=swapnode(ind1.tree,x1,xnode2); % substitute crosspoint1 by xnode2
%ind2.tree=swapnode(ind2.tree,x2,xnode1); % substitute crosspoint2 by xnode1

ind1.str=tree2str(ind1.tree);
ind1.id=[];
ind1.origin=indorigin;
ind1.parents=[pop(i1).id,pop(i2).id];
ind1.xsites=[x1,x2];
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

ind2.str=tree2str(ind2.tree);      	
ind2.id=[];
ind2.origin=indorigin;
ind2.parents=[pop(i2).id,pop(i1).id];
ind2.xsites=[x2,x1];
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

newinds=[newind1,newind2];
