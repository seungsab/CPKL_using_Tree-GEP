function mysizelist=findbranches(tree,s,mysizelist)
%FINDBRANCHES    Collects list of terminal branches from a GPLAB individual.
%   MYSIZELIST=FINDBRANCHES(TREE,SIZE,MYSIZELIST) returns a list of
%   node ids of terminal branches (terminals or branches of depth 2)
%   no larger than SIZE, where size is the number of nodes.
%
%   Input parameter MYSIZELIST is not set (is empty) on first call.
%
%   Input arguments:
%      TREE - the tree where to look for terminal branches (struct)
%      SIZE - maximum size of each terminal branch collected (integer)
%      MYSIZELIST - current terminal branch list (array)
%   Output arguments:
%      MYSIZELIST - final terminal branch list (array)
%
%   Notes:
%      FINDBRANCHES is a recursive function.
%
%   See also LIGHTMUTATION, MUTOPEQ
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

if ~exist('mysizelist')
   mysizelist=[]; 
end

% allow the selection of terminals only if the removal size is exactly 1
if s==1 
    allowterminals=1;
else
    allowterminals=0;
end

% is this (tree) a terminal branch?
% (a branch is terminal if all its kids are terminal)
% (a terminal node is also a terminal branch)
nonterminalkids={};
for i=1:length(tree.kids)
    if ~isempty(tree.kids{i}.kids)
        % save all non-terminal kids to visit later:
        nonterminalkids{end+1}=tree.kids{i};
    end
end

if isempty(nonterminalkids)
    % this (tree) is a terminal branch - collect it if its size is good:
    if (tree.nodes==1 && allowterminals) || (tree.nodes>1 && tree.nodes<=s) 
        mysizelist(end+1)=tree.nodeid;
    end
else
    % this (tree) is not a terminal branch, visit all its kids:
    for i=1:length(nonterminalkids)
        mysizelist=findbranches(nonterminalkids{i},s,mysizelist);
    end
end
