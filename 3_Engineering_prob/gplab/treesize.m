function ind=treesize(ind,params,data,terminals,varsvals)
%TREESIZE    Returns the tree size of a GPLAB individual.
%   TREESIZE(INDIVIDUAL,PARAMS,DATA,TERMINALS,VARSVALS) returns
%   the tree size of INDIVIDUAL, measured as the number of nodes.
%
%   Input arguments:
%      INDIVIDUAL - the individual whose fitness is to measure (struct)
%      PARAMS - the current running parameters (struct)
%      DATA - the dataset on which to measure the fitness (struct)
%      TERMINALS - the variables to set with the input dataset (cell array)
%      VARSVALS - the string of the variables of the fitness cases (string)
%   Output arguments:
%      INDIVIDUAL - the individual whose tree size was measured (struct)
%
%   See also NODES
%
%   Copyright (C) 2011-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

if isempty(ind.nodes)
    ind.nodes=nodes(ind.tree);
end

