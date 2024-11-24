function newinds=gdeCM(pop,params,state,data,indlist)
%GDECM    Creates new GPLAB individual by determining Center of Mass.
%   NEWIND=GDECM(POPULATION,PARAMS,STATE,DATA,PARENTS) returns a new
%   individual that is the Center of Mass of the PARENTS.
%
%   Input arguments:
%      POPULATION - the population where the parents are (array)
%      PARAMS - the parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      DATA - the current dataset for the algorithm to run (array)
%      PARENTS - the indices of the parents in POPULATION (1x2 matrix)
%   Output arguments:
%      NEWIND - the newly created individual (struct)
%
%   Notes:
%      GDECM is NOT currently used by GPLAB - it was meant for a NM implementation.
%
%   References:
%   Moraglio A, Silva S (2011). Geometric Nelder-Mead Algorithm on the Space
%   of Genetic Programs. In Proceedings of the Genetic and Evolutionary Computation
%   Conference (GECCO 2011), 1307–1314.
%
%   See also CENTERMASS, CROSSOVER
%
%   Copyright (C) 2011-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

for i=1:length(indlist)
    ix(i)=indlist(i);
    indx(i)=pop(ix(i));
    treex(i)=indx(i).tree;
end    
ind1=indx(1); % just to use the struct fields, as in crossover and others

[unique_ops,unique_i]=unique(state.functions(:,1),'legacy'); 
% list of possible operators and their arities:
possiblenodes=cell(length(unique_ops),2);
possiblenodes(:,1)=unique_ops;
possiblenodes(:,2)=state.functions(unique_i,2);

ind1.tree=centermass(treex,possiblenodes,0); % 0 = id of last node created

ind1.str=tree2str(ind1.tree);
ind1.id=[];
ind1.origin='gdeCM';
ind1.parents=[pop(ix).id];
ind1.xsites=[];
ind1.fitness=[];
ind1.adjustedfitness=[];
ind1.result=[];
ind1.testfitness=[];
ind1.testadjustedfitness=[];
ind1.level=[];
ind1.nodes=ind1.tree.nodes;
ind1.introns=[];
   
newinds=[ind1];
