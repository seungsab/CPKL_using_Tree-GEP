function newinds=gdeER(pop,params,state,data,indlist)
%GDEER    Creates a new GPLAB individual by GDE extension ray.
%   NEWIND=GDEER(POPULATION,PARAMS,STATE,DATA,PARENTS) returns one new
%   individual created by GDE extension ray of two parents.
%   (GDE = Geometric Differential Evolution)
%
%   Input arguments:
%      POPULATION - the population where the parents are (array)
%      PARAMS - the parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      DATA - the dataset on which the algorithm runs (struct)
%      PARENTS - the indices of the parents in POPULATION (integer)
%   Output arguments:
%      NEWIND - the newly created individual (struct)
%
%   References:
%   Moraglio A, Silva S (2010). Geometric Differential Evolution on
%   the Space of Genetic Programs, Proceedings of EuroGP-2010, 171-183.
%
%   See also GDE1,GDE2
%
%   Copyright (C) 2010-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

i1=indlist(1);
i2=indlist(2);

ind1=pop(i1);
ind2=pop(i2);


% weights W to decide whether the first or second parent contributes
% with genetic material (see Moraglio's refs)
%W1=1/(1+params.gdeF);
%W2=1-W1;
% this is decided in the call to this function, in params.gdeW1 and gdeW2

% probability P (see Moraglio's refs)
d_ab=shd(ind1.tree,ind2.tree);
if d_ab==1
    P=1;
else
    d_bc=d_ab*params.gdeW1/params.gdeW2;
    P=d_bc/(1-d_ab);
end

if isempty(ind1.nodes)
   ind1.nodes=nodes(ind1.tree);
end
if isempty(ind2.nodes)
   ind2.nodes=nodes(ind2.tree);
end

% Identify possible cross points (ids, arities, etc) in each parent:
[cp_ids1,cp_ids2,cp_arities1,cp_arities2,cp_ops1,cp_ops2,cp_nodes1,cp_nodes2,cp_equaltrees]=crosspoints(ind1.tree,ind2.tree,[],[],[],[],{},{},[],[],[]);
% (start the recursive function with empty lists)

% Store original cross point ids:
ocp_ids1=cp_ids1;
ocp_ids2=cp_ids2;

% To begin with, the only offspring is equal to the second parent:
% (so the only modified tree will be ind2.tree)

% For each point decide whether to change material in the offspring:
% (flip as many coins as possible cross points)
ncp=length(cp_ids1);
r=rand(1,ncp);
i=1;
while i<=ncp
    if r(i)<P
        % change:
        if strcmp(params.homocrosstype,'gdeER2') && cp_equaltrees(i) %er2
            % if it's a terminal, simply swap with newnodefunc:
            if cp_arities1(i)==0
                ind2.tree=newnodefunc(ind2.tree,cp_ids2(i),state);
                % (this will place a new function node with the same arity)
            else
                % make new random tree of same size:
                newtree=maketree(cp_nodes2(i),state.functions,state.arity,1,'2',[],cp_ids2(i)-1);
                %(1 means exact level, '2' means restrict nodes, not depth)
                nind2tree=swapnodes(ind2.tree,newtree,cp_ids2(i),1);
                % update the lists of crossing points, to remove the ones
                % that were inside the replaced subtree (a and b delimit them):
                %a=i+1;
                %b=i+cp_nodes1(i)-1;
                %ff=find(or(cp_ids1<a,cp_ids1>b));
                ff=[1:i i+cp_nodes1(i):ncp];
                ncp=length(ff);
                cp_ids1=cp_ids1(ff);
                cp_ids2=cp_ids2(ff);
                cp_arities1=cp_arities1(ff);
                cp_arities2=cp_arities2(ff);
                cp_ops1=cp_ops1(ff);
                cp_ops2=cp_ops2(ff);
                cp_nodes1=cp_nodes1(ff);
                cp_nodes2=cp_nodes2(ff);
                cp_equaltrees=cp_equaltrees(ff);
                % update the remaining ids to swap, because of possible size change
                % (maketree cannot guarantee exact size)
                df=nind2tree.nodes-ind2.tree.nodes;
                cp_ids2(i+1:end)=cp_ids2(i+1:end)+df;
                ind2.tree=nind2tree;
            end
        elseif cp_arities1(i)==cp_arities2(i) % type 1
            if strcmp(cp_ops1{i},cp_ops2{i})
                % create random node to swap, then swap
                % select random function:
                f_funcs=state.functions(:,1);
                f_arities=state.arity;
                f_funcs=f_funcs(find(f_arities==cp_arities2(i)));
                rr=intrand(1,length(f_funcs));
                f=cell2mat(f_funcs(rr));
                % swap:
                ind2.tree=newnodefunc(ind2.tree,cp_ids2(i),state);
            end
        else % type 2 % impossible case of different arity and same ops
            % IF THIS CASE BECOMES POSSIBLE, CHECK CODE!!!!
            if strcmp(cp_ops1{i},cp_ops2{i})
                % create random branch to swap, then swap
                % make new branch, no deeper than average depth:
                dp=ceil((cp_nodes1(i)+cp_nodes2(i))/2);
                % (ceil compensates the fact that maketree may not achieve
                % the exact number of intended nodes, just close)
                newtree=maketree(dp,state.functions,state.arity,1,'2',[],cp_ids2(i)-1);
                %(1 means exact level, '2' means restrict nodes, not depth)
                nind2tree=swapnodes(ind2.tree,newtree,cp_ids2(i),1);
                % update the remaining ids to swap, because of the size change
                df=nind2tree.nodes-ind2.tree.nodes;
                cp_ids2(i+1:end)=cp_ids2(i+1:end)+df;
                ind2.tree=nind2tree;
            end
        end
    end
    i=i+1;
end

   
ind2.str=tree2str(ind2.tree);      	
ind2.id=[];
ind2.origin=params.homocrosstype;
ind2.parents=[pop(i2).id,pop(i1).id];
ind2.xsites=[ocp_ids2',ocp_ids1'];
ind2.fitness=[];
ind2.adjustedfitness=[];
ind2.result=[];
ind2.testfitness=[];
ind2.testadjustedfitness=[];
ind2.level=[];
ind2.nodes=ind2.tree.nodes;
ind2.introns=[];   

newind2=ind2;

newinds=newind2;
