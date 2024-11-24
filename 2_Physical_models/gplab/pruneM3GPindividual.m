function ind=pruneM3GPindividual(ind,params,data,state)
%PRUNEM3GPINDIVIDUAL    Prune useless branches of a M3GP GPLAB individual.
%   PRUNEDIND=PRUNEM3GPINDIVIDUAL(INDIVIDUAL) returns the new individual
%   resulting from removing the tree branches that do not help fitness.
%
%   Input arguments:
%      INDIVIDUAL - the individual whose tree is to be pruned (struct)
%   Output arguments:
%      PRUNEDIND - the resulting pruned tree (struct)
%
%   See also PRUNEM3GPDIMENSION
%
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox
%return
npruned=0;
for k=1:length(ind.tree.kids)
    prunedind=ind;
    prunedtree=pruneM3GPdimension(ind.tree,ind.tree.kids{k-npruned}.nodeid);
    prunedind.tree=prunedtree;
    prunedind.str=tree2str(prunedind.tree);
    prunedind.nodes=prunedtree.nodes;
    % clear some fields before calculating fitness: % mapping is enough
    prunedind.mapping=[];
    [tmpind,tmpstate]=calcfitness(prunedind,params,data,state,0);
    if or(and(params.lowerisbetter,tmpind.adjustedfitness<=ind.adjustedfitness),and(~params.lowerisbetter,tmpind.adjustedfitness>=ind.adjustedfitness))
        ind=tmpind;
        npruned=npruned+1;
    end
end

ind.pruned=1;

