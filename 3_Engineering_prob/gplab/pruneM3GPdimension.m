function prunedtree=pruneM3GPdimension(tree,x)
%PRUNEM3GPDIMENSION    Prune entire branch from the root of a M3GP GPLAB tree.
%   PRUNEDTREE=PRUNEM3GPDIMENSION(TREE,X) returns the new tree (TREE) resulting
%   from removing the branch with root X. Nodes are numbered depth-first.
%
%   Input arguments:
%      TREE - the tree to be pruned (struct)
%      X - the number of the node to be remove with its branch (integer)
%   Output arguments:
%      PRUNEDTREE - the resulting pruned tree (struct)
%
%   Notes: Not a recursive function. Cannot prune internal branches.
%
%   See also SWAPNODES, MUTATION
%
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

found=0;
sizebranch=0;
prunedtree=tree;
prunedtree.kids={};

for k=1:length(tree.kids)
    if found % branches after pruned branch
        prunedtree.kids{k-1}=updatenodeids(tree.kids{k},-sizebranch);
        % k-1 because there is one less branch, updatenodeids to reflect it
    else
        if tree.kids{k}.nodeid==x % here is the branch to prune
            found=1;
            sizebranch=tree.kids{k}.nodes;
            % do not copy this branch
            prunedtree.nodes=tree.nodes-sizebranch; % size reduction
            prunedtree.maxid=tree.maxid-sizebranch; % size reduction
        else % branches before branch to prune
            prunedtree.kids{k}=tree.kids{k}; % copy without changes
        end
    end
end
