function [validind,state]=DynOpEq(ind,pop,params,state,data,parentindices)
%DYNOPEQ    Applies operator equalisation filters to a new GPLAB individual.
%   [VALIDIND,STATE]=DYNOPEQ(IND,POP,PARAMS,STATE,DATA,PARENTS) decides
%   whether an individual (IND) can be accepted, considering the
%   current target distribution of length classes and fitness statistics
%   for each class. If not accepted, returns empty [].
%
%   Input arguments:
%      IND - individual to be validated (array)
%      POPULATION - the current population of the algorithm (array)
%      PARAMS - the running parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      DATA - the dataset for use in the algorithm (struct)
%      PARENTS - the indices of the parents of IND (matrix)
%   Output arguments:
%      VALIDIND - valid individual (IND or one of its parents) (array)
%      STATE - the updated state of the algorithm (struct)
%
%   References:
%      Silva S, Dignum S, "Extending Operator Equalisation: Fitness Based
%      Self Adaptive Length Distribution for Bloat Free GP", EuroGP-2009.
%
%   See also VALIDATEINDS, STRICTDEPTH, DYNNODES, ... (the other filters)
%
%   Copyright (C) 2008-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% tree size is needed:
if isempty(ind.nodes)
	ind.nodes=nodes(ind.tree);
end

% fitness is needed:
if isempty(ind.fitness)
    [ind,state]=calcfitness(ind,params,data,state,0);
end

% determine appropriate bin for this tree size:
bin=opeq_getbin(ind.nodes,params.opeq_binsize);

% reject is the default:
accept=0;

% if bin exists...
if length(state.opeq.exists)>=bin
    

    % if bin is not full, accept unconditionaly
    if ~state.opeq.full(bin) 
        accept=1;

    % if bin is full but this individual will improve mean/best fitness...
    % (current implementation: best)
    % or bin has no individuals...
    % accept
    %elseif state.opeq.current_occupancy(bin)==0 || ((params.lowerisbetter && mean(state.opeq.fitnesses{bin})<mean([state.opeq.fitnesses{bin} ind.fitness])) || (~params.lowerisbetter && mean(state.opeq.fitnesses{bin})>mean([state.opeq.fitnesses{bin} ind.fitness])))
    elseif state.opeq.current_occupancy(bin)==0 || ((params.lowerisbetter && ind.fitness<state.opeq.bestfitness(bin)) || (~params.lowerisbetter && ind.fitness>state.opeq.bestfitness(bin)))
        accept=1;
    end

% if bin does not exist...
else
    % if bin did not exist before... accept with restriction:
    % new individual must be better than the current best
    if (params.lowerisbetter && ind.fitness<min(state.opeq.bestfitness)) || (~params.lowerisbetter && ind.fitness>max(state.opeq.bestfitness))
        accept=1;
        state.opeq.exists(1:bin)=1;
        state.opeq.target_occupancy(bin)=1;
        state.opeq.current_occupancy(bin)=0;
    end
end

if accept
    state.opeq.current_occupancy(bin)=state.opeq.current_occupancy(bin)+1;
    if length(state.opeq.fitnesses)<bin
        state.opeq.fitnesses{bin}=[];
    end
    state.opeq.fitnesses{bin}(end+1)=ind.fitness;
    if state.opeq.current_occupancy(bin)>=state.opeq.target_occupancy(bin)
        state.opeq.full(bin)=1;
    end
    if params.lowerisbetter
        state.opeq.bestfitness(bin)=min(state.opeq.fitnesses{bin});
        state.opeq.worstfitness(bin)=max(state.opeq.fitnesses{bin});
    else
        state.opeq.bestfitness(bin)=max(state.opeq.fitnesses{bin});
        state.opeq.worstfitness(bin)=min(state.opeq.fitnesses{bin});
    end
    
    validind=ind;

else

    if length(state.opeq.rejected)<bin
        state.opeq.rejected(bin)=0;
        state.opeq.rejectedsizes(bin)=0;
    end
    state.opeq.rejected(bin)=state.opeq.rejected(bin)+1;
    state.opeq.rejectedsizes(bin)=state.opeq.rejectedsizes(bin)+ind.nodes;
    
    % return empty
    validind=[];
end        

