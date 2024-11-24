function [validind,state]=dynnodes(ind,pop,params,state,data,parentindices)
%DYNNODES    Applies dynamic size filters to a new GPLAB individual.
%   [VALIDIND,STATE]=DYNNODES(IND,POP,PARAMS,STATE,DATA,PARENTS)
%   tests if an individual (IND) conforms to the dynamic maximum
%   size (number of nodes) rules. If not, returns one of its parents.
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
%      Silva S. and Almeida J. (2003) "Dynamic maximum tree depth - a
%      simple technique for avoiding bloat in tree-based GP". GECCO-2003.
%      Silva S. and Costa E. (2004) "Dynamic limits for bloat control -
%      variations on size and depth". GECCO-2004.
%
%   See also VALIDATEINDS, DYNDEPTH, STRICTNODES, ... (the other filters)
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% tree nodes is needed:
if isempty(ind.nodes)
	ind.nodes=nodes(ind.tree);
end
      
% if more nodes than state.maxlevel, measure fitness:
if ind.nodes>state.maxlevel
   [ind,state]=calcfitness(ind,params,data,state,0);
   % if this individual is the best so far, increase dynamic nodes:
   if ~isempty(state.bestsofar) && ((params.lowerisbetter && ind.fitness<state.bestsofar.fitness) || (~params.lowerisbetter && ind.fitness>state.bestsofar.fitness))
      state.lastid=state.lastid+1;
		ind.id=state.lastid;
      state.maxlevel=ind.nodes;
      state.bestsofar=ind;
      if ~strcmp(params.output,'silent')
      	displaystatus('     (increasing maximum nodes to %d)\n',state.maxlevel);
      end
        
    % else if this individual is not the best so far, check the size of the parents:  
   else
      mnodes=max([pop(parentindices).nodes]);
      % if both parents within limits (mnodes<=state.maxlevel), accept
      % offpring unconditionally (allow exploration of search space)
      % if parents were first generation (isempty(mnodes)) or this ind is
      % no bigger than biggest parent (ind.nodes<=mnodes), also accept
      if isempty(mnodes) || ind.nodes<=mnodes || mnodes<=state.maxlevel
          state.lastid=state.lastid+1;
          ind.id=state.lastid;
      else % else reject (substitute by a parent):
         ind=pop(parentindices(1));
         if length(state.rejected)<ind.nodes
		state.rejected(ind.nodes)=1;
	 else
	        state.rejected(ind.nodes)=state.rejected(ind.nodes)+1;
	 end
      end
   end        
           
end

validind=ind;
