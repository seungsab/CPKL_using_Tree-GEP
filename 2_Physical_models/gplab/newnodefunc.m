function tree=newnodefunc(tree,x,state)
%NEWNODEFUNC    Replaces one node in a GPLAB tree with a node of the same arity.
%   TREE=NEWNODEFUNC(TREE,X,STATE) returns TREE with node X replaced
%   by a new node with the same arity, from the ones available in STATE.
%   
%   Input arguments:
%      TREE - the original tree to replace a node (struct)
%      X - the id of the node to be replaced (integer)
%      STATE - the current state of the algorithm (struct)
%   Output arguments:
%      TREE - the tree with the replaced node (struct)
%
%   Note: This is a recursive function.
%
%   See also SWAPNODEFUNC, GDEER
%
%   Copyright (C) 2010-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox


% search node in tree. when found, search for function of same arity:
if tree.nodeid==x
    f_funcs=state.functions(:,1);
    f_arities=state.arity;
    f_funcs=f_funcs(find(f_arities==length(tree.kids)));
    rr=intrand(1,length(f_funcs));
    f=cell2mat(f_funcs(rr));
    tree.op=f;
else
    for i=1:length(tree.kids)
        if tree.kids{i}.maxid>=x
            %this is the kid where the new func will be inserted:
            tree.kids{i}=newnodefunc(tree.kids{i},x,state);
            break
        end
    end
end