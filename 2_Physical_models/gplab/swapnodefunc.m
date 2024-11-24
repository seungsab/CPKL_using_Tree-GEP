function [tree1,tree2]=swapnodefunc(tree1,tree2,x1,x2,sameshape)
%SWAPNODEFUNC    Swaps a node in one GPLAB tree with a node in another GPLAB tree.
%   [TREE1,TREE2]=SWAPNODEFUNC(TREE1,TREE2,X1,X2,SAMESHAPE) returns TREE1 and TREE2
%   where the original respective nodes X1 and X2 have been swapped with one another.
%   
%   Input arguments:
%      TREE1 - the first tree where swapping will be made (struct)
%      TREE2 - the second tree where swapping will be made (struct)
%      X1 - the id of the node to be swapped in TREE1 (integer)
%      X2 - the id of the node to be swapped in TREE2 (integer)
%      SAMESHAPE -  always 1 for now (boolean)
%   Output arguments:
%      TREE1 - the first tree with a swapped node (struct)
%      TREE2 - the second tree with a swapped node (struct)
%
%   Note: This is a recursive function.
%
%   See also NEWNODEFUNC, HOMOCROSS
%
%   Copyright (C) 2010-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% search node in tree1:
% (in tree2 the node will be in the same position (sameshape until this point)
%  or tree2 will be a single node)
if tree1.nodeid==x1
    tmp_op=tree1.op;
    tree1.op=tree2.op;
    tree2.op=tmp_op;
else
    for i=1:length(tree1.kids)
        if tree1.kids{i}.maxid>=x1
            %this is the kid where the swap will occur:
            if sameshape
                [tree1.kids{i},tree2.kids{i}]=swapnodefunc(tree1.kids{i},tree2.kids{i},x1,x2,sameshape);
            else
                [tree1.kids{i},tree2]=swapnodefunc(tree1.kids{i},tree2,x1,1,sameshape);
            end
            break
        end
    end
end