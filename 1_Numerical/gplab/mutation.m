function newind=mutation(pop,params,state,data,i)
%MUTATION    Creates a new individual for GPLAB by mutation.
%   NEWIND=MUTATION(POPULATION,PARAMS,STATE,DATA,PARENT) returns a new
%   individual created by substituting a random subtree of PARENT by a
%   new randomly created tree, with the same depth/size restrictions
%   as the initial random trees.
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
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

ind=pop(i);

% calculate number of nodes (we need it to pick a random branch)
if isempty(ind.nodes)
   ind.nodes=nodes(ind.tree);
end

if params.M3GP
    % decide which mutation:
    r=rand;
    if r<1/3 % normal mutation, just do not touch the root node
        x=intrand(2,ind.nodes); % mutation point
        indorigin='M3GPmutation_normal';
    elseif r<2/3 % add a dimension as a new last branch of the root node
        newtree=maketree(state.iniclevel,state.functions,state.arity,0,params.depthnodes,[],ind.tree.nodes);
        if existsM3GPbranch(newtree,ind.tree)
            x=0; % branch already exists, this will signal that nothing was done
        else
            %newtree=updatenodeids(newtree,newtree.nodes);
            x=1; % only useful to save mutation point, in the end
            ind.tree.kids{end+1}=newtree;
            ind.tree.nodes=ind.tree.nodes+newtree.nodes;
            ind.tree.maxid=ind.tree.maxid+newtree.nodes;
        end
        indorigin='M3GPmutation_adddimension';
    else % remove a random dimension, by removing a branch of the root node
        nkids=length(ind.tree.kids);
        if nkids>1
            listnodes=zeros(1,nkids);
            for k=1:nkids
                listnodes(k)=ind.tree.kids{k}.nodeid;
            end
            x=listnodes(intrand(1,nkids)); % removal point, a dimension root
            ind.tree=pruneM3GPdimension(ind.tree,x);
        else % else do nothing, do not remove the only dimension
            x=0; % this will signal that nothing was done
        end 
        indorigin='M3GPmutation_removedimension';
    end
else
    x=intrand(1,ind.nodes); % mutation point
    indorigin='mutation';
end

if or(~params.M3GP,strcmp(indorigin,'M3GPmutation_normal'))

    newtree=maketree(state.iniclevel,state.functions,state.arity,0,params.depthnodes,[],x-1);
    % (the maximum size of the new branch is the same as the initial random trees)
    % (0 means no exact level)

    % swap old branch with new branch in only one step, as if this were
    % crossover (but discard the resulting nind):
    nind.tree=newtree;
    ind.tree=swapnodes(ind.tree,nind.tree,x,1);
    %ind.tree=swapnode(ind.tree,x,newtree);
end

ind.id=[];
ind.origin=indorigin;
ind.parents=[pop(i).id];
ind.xsites=[x];
ind.str=tree2str(ind.tree);
ind.fitness=[];
ind.adjustedfitness=[];
ind.result=[];
ind.testfitness=[];
ind.testadjustedfitness=[];
ind.nodes=ind.tree.nodes;
ind.introns=[];
ind.level=[];

ind.dimensions=[];
ind.pruned=0;
ind.mapping=[];
ind.covariancematrix={};
ind.inversematrix={};
ind.centroids=[];

ind.distancematrix=[];

ind.countclusters=[];
ind.mapclusters=[];

newind=ind;
