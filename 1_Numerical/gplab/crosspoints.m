function [ids1,ids2,arities1,arities2,ops1,ops2,nodes1,nodes2,equaltrees]=crosspoints(tree1,tree2,ids1,ids2,arities1,arities2,ops1,ops2,nodes1,nodes2,equaltrees)
%CROSSPOINTS    Identifies possible homologous crossing nodes between two GPLAB trees.
%   [IDS1,IDS2,ARS1,ARS2,OPS1,OPS2,NS1,NS2,EQTREES]=CROSSPOINTS(T1,T2,IDS1,IDS2,ARS1,ARS2,OPS1,OPS2,NS1,NS2,EQTREES)
%   returns the following information regarding all possible crossing points between
%   two trees: the node identifiers, the respective arities, the respective operators,
%   the respective branch sizes, and whether both trees are equal until this crossing
%   point. The crossing points must be homologous. The info is returned for each tree
%   separately (except the info on whether they are equal), in arrays and cell arrays.
%
%   Input arguments:
%      T1 - the first tree where to look for crossing points (struct)
%      T2 - the second tree where to look for crossing points (struct)
%      IDS1 - the ids of nodes ASPCP1 (see notes below) (1xn array)
%      IDS2 - the ids of nodes ASPCP2 (see notes below) (1xn array)
%      ARS1 - the arties of nodes ASPCP1 (see notes below) (1xn array)
%      ARS2 - the arities of nodes ASPCP2 (see notes below) (1xn array)
%      OPS1 - the operators of nodes ASPCP1 (see notes below) (1xn cell array)
%      OPS2 - the operators of nodes ASPCP2 (see notes below) (1xn cell array)
%      NS1 - the size of the branches with root ASPCP1 (see notes below) (1xn cell array)
%      NS1 - the size of the branches with root ASPCP2 (see notes below) (1xn cell array)
%      EQTREES - whether both trees are equal until this point (1xn array)
%   Output arguments:
%      Same as input arguments except T1 and T2
%
%   Notes:
%      ASPCP1 = already set as possible crossing points in T1
%      ASPCP2 = already set as possible crossing points in T2
%      This is a recursive function.
%      On the first call all parameters except T1 and T2 are empty.
%
%   See also HOMOCROSS, GDEER
%
%   Copyright (C) 2010-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox


% ids:
ids1(end+1)=tree1.nodeid;
ids2(end+1)=tree2.nodeid;
% arities:
arities1(end+1)=length(tree1.kids);
arities2(end+1)=length(tree2.kids);
% operators:
ops1{end+1}=tree1.op;
ops2{end+1}=tree2.op;
% nodes:
nodes1(end+1)=tree1.nodes;
nodes2(end+1)=tree2.nodes;
% equal trees:
equaltrees(end+1)=strcmp(ops1{end},ops2{end});
l=length(equaltrees);

if arities1(end)==arities2(end)
    for i=1:length(tree1.kids)
        cl=length(equaltrees); % current length
        [ids1,ids2,arities1,arities2,ops1,ops2,nodes1,nodes2,equaltrees]=crosspoints(tree1.kids{i},tree2.kids{i},ids1,ids2,arities1,arities2,ops1,ops2,nodes1,nodes2,equaltrees);
        equaltrees(l)=and(equaltrees(l),equaltrees(cl+1));
    end 
end