function [cmtree,lastnode]=centermass(ptrees,possiblenodes,lastnode)
%CENTERMASS    Creates new GPLAB tree by calculating the Center of Mass of other trees.
%   CENTERMASS(PTREES,NODES,LASTNODE) creates a new tree that represents the Center of
%   Mass of other trees. For each node to create on the new tree, NODES stores the
%   operators and respective arities found on the homologous positions on the other trees.
%
%   Additional input parameter (not set in the first call)
%   is LASTNODE, the id of the last node created in the tree.
%
%   Additional output argument (essential in recursiveness)
%   is LASTNODE, the id of the last node created in the tree.
%
%   Input arguments:
%      PTREES - parent trees (1xn matrix)
%      NODES - column 1 contains possible operators, column 2 their arities (nx2 matrix)
%      LASTNODE - the id of the last node created in the tree (integer)
%   Output arguments:
%      CMTREE - the new tree (struct)
%      LASTNODE - the id of the last node created in the tree (integer)
%
%   Notes:
%      CENTERMASS is a recursive function.
%      CENTERMASS is NOT currently used by GPLAB - it was meant for a NM implementation.
%
%   References:
%   Moraglio A, Silva S (2011). Geometric Nelder-Mead Algorithm on the Space
%   of Genetic Programs. In Proceedings of the Genetic and Evolutionary Computation
%   Conference (GECCO 2011), 1307–1314.
%
%   See also GDECM
%
%   Copyright (C) 2011-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

thisnode=lastnode+1;

% find most common node:

% initialize counting:
countpossiblenodes=zeros(length(possiblenodes(:,1)),1);

% count occurrences of each operator, with a for cycle (cell arrays...)
for t=1:length(ptrees)
    f=find(strcmp(possiblenodes(:,1),ptrees(t).op));
    countpossiblenodes(f)=countpossiblenodes(f)+1;
end

% use node with maximum count:
% (f identifies the index where it can be found in possiblenodes)
f=find(countpossiblenodes==max(countpossiblenodes));
f=f(intrand(1,length(f))); % because there may be more than one


cmtree.op=possiblenodes{f,1};
cmtree.kids=[];
cmtree.nodeid=thisnode;

% arity of the chosen node:
k=possiblenodes{f,2};

if k~=0
    % generate k branches:
    for i=1:k
        % build the list of "parent subtrees" for the ith subtree:
        if exist('kptrees','var')
            clear kptrees
        end
        for p=1:length(ptrees)
            % if this tree has the right arity (k),
            % add its ith subtree to the list
            if length(ptrees(p).kids)==k
                if ~exist('kptrees','var')
                    kptrees(1)=ptrees(p).kids{i};
                else
                    kptrees(end+1)=ptrees(p).kids{i};
                end
            end
        end
        % calculate the centre of mass of the list of parent subtrees
        [kcmtree,lastnode]=centermass(kptrees,possiblenodes,thisnode);
        cmtree.kids{i}=kcmtree;
        thisnode=lastnode+1;
    end 
end

cmtree.nodes=thisnode-cmtree.nodeid+1;
cmtree.maxid=lastnode+1;
