function found=existsM3GPbranch(mybranch,tree)
%EXISTSM3GPBRANCH    Check if a branch already exist in a M3GP GPLAB tree.
%   FOUND=EXISTSM3GPBRANCH(BRANCH,TREE) returns true (1) if BRANCH is one
%   of the branches of the root of TREE.
%
%   Input arguments:
%      BRANCH - the branch to search in TREE (struct)
%      TREE - the tree where to look for BRANCH (struct)
%   Output arguments:
%      FOUND - whether the branch exists in tree (boolean)
%
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

found=false;
for k=1:length(tree.kids)
    if found
        return
    end
    if mybranch.nodes==tree.kids{k}.nodes % number of nodes is equal
        if strcmp(mybranch.op,tree.kids{k}.op) % main operator is equal
            if and(isempty(mybranch.kids),isempty(tree.kids{k}.kids)) % kids are both empty
                found=true;
            else % kids not empty, check them searching for differences
                found=true; % assume all kids are equal, change if a difference is found
                for bk=1:length(mybranch.kids)
                    if ~existsM3GPbranch(mybranch.kids{bk},tree.kids{k})
                        found=false;
                        break % as soon as one is different, no need to continue searching
                    end
                end
            end
        end
    end
end
