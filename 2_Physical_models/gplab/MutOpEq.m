function [validind,state]=MutOpEq(ind,pop,params,state,data,parentindices)
%MUTOPEQ    Applies operator equalisation filters to a new GPLAB individual.
%   [VALIDIND,STATE]=MUTOPEQ(IND,POP,PARAMS,STATE,DATA,PARENTS) accepts
%   an individual (IND), considering the current target distribution of
%   length classes and fitness statistics for each class. Or, it adapts
%   the individual (using light mutations) to fit the target.
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
%      Silva S, Vanneschi L, "Operator Equalisation, Bloat and Overfitting 
%      - A Study on Human Oral Bioavailability Prediction", GECCO-2009.
%
%   See also VALIDATEINDS, STRICTDEPTH, DYNNODES, ... (the other filters)
%
%   Copyright (C) 2008-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% tree size is needed:
if isempty(ind.nodes)
	ind.nodes=nodes(ind.tree);
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

    % do not accept more individuals in a bin that is full
    % (accept if bin has no individuals)
    %elseif state.opeq.current_occupancy(bin)==0 || ((params.lowerisbetter && mean(state.opeq.fitnesses{bin})<mean([state.opeq.fitnesses{bin} ind.fitness])) || (~params.lowerisbetter && mean(state.opeq.fitnesses{bin})>mean([state.opeq.fitnesses{bin} ind.fitness])))
    %elseif state.opeq.current_occupancy(bin)==0 || ((params.lowerisbetter && ind.fitness<state.opeq.bestfitness(bin)) || (~params.lowerisbetter && ind.fitness>state.opeq.bestfitness(bin)))
    elseif state.opeq.current_occupancy(bin)==0
        accept=1;
    end

% if bin does not exist...
else
    % if bin did not exist before... accept with restriction:
    % new individual must be better than the current best
    if isempty(ind.fitness)
        [ind,state]=calcfitness(ind,params,data,state,0);
    end
    if (params.lowerisbetter && ind.fitness<min(state.opeq.bestfitness)) || (~params.lowerisbetter && ind.fitness>max(state.opeq.bestfitness))
        accept=1;
        state.opeq.exists(1:bin)=1;
        state.opeq.target_occupancy(bin)=1;
        state.opeq.current_occupancy(bin)=0;
        fprintf('     (created new bin for new best individual of size %d)\n',ind.nodes);
    end
end

% we adapt the individual using weak mutations, until acceptance
% --------------------------------------------------------------

if ~accept

    % pick closest non-full bin:
    % (find non-full bins that have target_occupancy>0)
    f=find(~state.opeq.full & state.opeq.target_occupancy>0);
    if isempty(f)
        % it may happen that the target does not sum to popsize, because
        % of rounding in opeq_calctarget.m (or it may happen that we create
        % more than popsize individuals, because of the number of children
        % created by each genetic operator).
        % so, here we accept additional individuals when the target is full.
        % (find any bins that have target_occupancy>0)
        f=find(state.opeq.target_occupancy>0);
    end
    df=abs(f-bin); % distance of target bins to bin
    cb=find(df==min(df)); % find the smallest distance
    cb=f(cb(1)); % choose the smallest bin, in case there were several
        
    % calculate min and max size for the individual to fit this bin:
    % (save it in state.tmp to pass to the light mutation operator)
    targetsize.min=(cb-1)*params.opeq_binsize+1;
    targetsize.max=(cb-1)*params.opeq_binsize+params.opeq_binsize;
    state.tmp.targetsize=targetsize;
    
    % the next cycle can call the light mutation several times, as many as
    % necessary till acceptance. one acceptance criteria is, however, the
    % size not changing after the mutation, meaning it was not possible to
    % reach the intended size:
    while ~accept        
        % apply light mutation (to shrink or enlarge the individual):
        % (we pass only 'ind' instead of 'pop'. the index has to be 1)
        newind=lightmutation(ind,params,state,data,1);
        % the individual is accepted if its size is appropriate for bin,
        % or if its size has not changed at all
        if (newind.nodes>=targetsize.min && newind.nodes<=targetsize.max) || (newind.nodes==ind.nodes)
            accept=1;
        end
        ind=newind;
    end

end

% NOTE: we always reach this point with accept=1,
% but we may still have to mutate one more time:
% when the size of the individual has not changed because it cannot reach
% the correct size, but it is still outside the boundaries of the target

% determine the bin again, because it may have changed:
bin=opeq_getbin(ind.nodes,params.opeq_binsize);
    
% check if this bin is within the target:
if bin>length(state.opeq.current_occupancy)
    % it is not, so we must shrink the individual a bit more:
    state.tmp.targetsize.min=1;
    state.tmp.targetsize.max=length(state.opeq.current_occupancy);
    newind=lightmutation(ind,params,state,data,1);
    % in this case we gave light mutation total freedom to guarantee the
    % individual is now within the target, but the mutation is always light
    ind=newind;
    % determine the bin again, because it certainly changed:
    bin=opeq_getbin(ind.nodes,params.opeq_binsize);
end

if accept % always happens in MutOpEq
    
    % fitness is needed:
    if isempty(ind.fitness)
        [ind,state]=calcfitness(ind,params,data,state,0);
    end

    % and now fill the opeq stats arrays:
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
    
end        

validind=ind;
