function newind=lightmutation(pop,params,state,data,i)
%LIGHTMUTATION    Creates a new individual for GPLAB by light mutation.
%   NEWIND=LIGHTMUTATION(POPULATION,PARAMS,STATE,DATA,PARENT) returns a
%   new individual created by substituting the minimum amount of genetic
%   material in POPULATION(PARENT). The choice of the branch to change
%   can be directed by information contained in STATE.TMP.TARGETSIZE.
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
%   See also CROSSOVER, MUTATION, APPLYOPERATOR
%
%   References:
%   Vanneschi et al. Fitness Distance Correlation in Structural
%   Mutation Genetic Programming. EuroGP-2003.
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

ind=pop(i);

% current size:
if isempty(ind.nodes)
   ind.nodes=nodes(ind.tree);
end

% difference needed to achieve minimum allowed size on the intended bin:
d=abs(state.tmp.targetsize.min-ind.nodes);

% if individual is too large to fit bin, shrink it:
% (d+1 nodes will be removed, and one node will be added)
if ind.nodes>state.tmp.targetsize.min
    s_remove=d+1;
else
    % if individual is too small to fit bin, enlarge it:
    % (1 node will be removed, and one terminal branch will be added)
    s_remove=1;
end

% NOTE:
% regardless of the difference needed to achieve the intended size,
% only terminal branches will be inserted or removed in each mutation step 

% collect list of branches compatible with intended removal size:
% (only terminals or terminal branches will be selected)
% (the list mysizelist will contain the node ids of compatible branches)
mysizelist=findbranches(ind.tree,s_remove);
if isempty(mysizelist)
    newind=ind;
    return % do nothing, mutation is not possible
end
% then choose one branch at random, and it becomes the node to swap:
x=mysizelist(intrand(1,length(mysizelist)));

if s_remove==1
    % enlarge
    % this is done in several steps:
    
    % step 1: create a new terminal branch of depth 2:
    newtree=maketree(2,state.functions,state.arity,1,'1',[],x-1);
    % (2 is the depth we want for the new branch)
    % (1 means we aim at the exact depth)
    % ('1' means we are looking at depth, not size, for the new branch)
    % (x-1 is the last node that will remain unchanged in the tree that
    %  will be mutated. the root of the new branch will have id x)
    % (see MAKETREE)
    
    % step 2: one of the terminals of the new branch is replaced by the
    % terminal chosen on the original tree:
    % (just choose the first, it was randomly generated anyway)
    nind.tree=newtree;
    nind.tree=swapnodes(nind.tree,ind.tree,2,x);
    % (2 is the first terminal of the new branch -> 2nd node, not nodeid=2)
    % (this is mutation like other mutations, i.e., like crossover)
    
    % step 3: finally swap new modified branch with previously selected
    % node of the original tree:
    ind.tree=swapnodes(ind.tree,nind.tree,x,1);
    
else
    % shrink
    % this is done in several steps:
    
    %step 1: find the branch to remove:
    nodei=findnode(ind.tree,x);
    
    % step 2: choose one of its leaves randomly:
    nodei=nodei.kids{intrand(1,length(nodei.kids))};
    
    % step 3: replace the branch to remove with the chosen leaf:
    ind.tree=swapnodes(ind.tree,nodei,x,1);
    
end

ind.id=[];
ind.origin=[ind.origin '+lightmutation'];
%ind.parents=[pop(i).id]; % do not overwrite original origin
%ind.xsites=[x]; % do not overwrite original origin
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
