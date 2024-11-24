function [pop,state]=generation(pop,params,state,data,varargin)
%GENERATION    Creates a new generation for the GPLAB algorithm.
%   GENERATION(POPULATION,PARAMETERS,STATE,DATA) returns a new
%   population of individuals created by applying the genetic
%   operators to the current population. Also implements GDE (see refs).
%
%   GENERATION(POPULATION,PARAMETERS,STATE,DATA,NEWPOPSIZE)
%   creates a new population with NEWPOPSIZE individuals, which
%   can be different from the size of the current population.
%
%   [POPULATION,STATE]=GENERATION(...) also returns the updated
%   state of the algorithm.
%
%   Input arguments:
%      POPULATION - the current population (array)
%      PARAMETERS - the algorithm running parameters (struct)
%      STATE - the algorithm current state (struct)
%      DATA - the dataset on which the algorithm runs (struct)
%      NEWPOPSIZE - the number of individuals to create (integer)
%   Output arguments:
%      POPULATION - the new population (array)
%      STATE - the algorithm current state, now updated (struct)
%
%   See also GENPOP
%
%   References:
%   Moraglio A, Silva S (2010). Geometric Differential Evolution on
%   the Space of Genetic Programs, Proceedings of EuroGP-2010, 171-183.
%
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

if nargin==5
   newpopsize=varargin{1};
else
   newpopsize=state.popsize;
end

state.generation=state.generation+1;
if ~strcmp(params.output,'silent')
	displaystatus('- Creating generation %d -\n',state.generation);
end

% create a new generation:
      
% first set a temporary empty population, with newpopsize individuals...
tmppop=initpop(newpopsize); % tmppop is an empty population
currentsize=0; % number of non empty individuals in tmppop
      
% ...then apply genetic operators to pop and fill tmppop with the new individuals...
% (and eventually calculate stuff for variable operator probabilities)
parentspool=[];
while currentsize<newpopsize
   
   % first choose between reproduction and genetic operator:
   rrate=params.reproduction;
   if (rrate>0 && rand<rrate)
      % reproduction
      opnum=0;
   else
      % genetic operator - if there's only one don't bother to randomize
      if length(params.operatornames)==1
          opnum=1;
      else
          opnum=pickoperator(state); % pickoperator gives the index of the operator chosen
      end
   end
   [tmppop,newsize,parentspool,state]=applyoperator(pop,params,state,data,tmppop,currentsize,opnum,parentspool);
   
   % if operator probabilities are variable, call procedure to adapt them:
   if strcmp(params.operatorprobstype,'variable')
      [state,tmppop]=automaticoperatorprobs(tmppop,pop,params,state,data,currentsize,newsize);
   end
   
   currentsize=newsize;
end

% note that tmppop may have more individuals than needed, because the last operator applied
% may have produced more offspring than needed to fill the population. tmppop will be
% subject to applysurvival, where the worst individuals will be discarded. (even if there is
% no elitism, if tmppop has more individuals than needed, the worst will be discarded.)

% ...measure fitness on the tmppop individuals with empty fitness...
[tmppop,state]=calcpopfitness(tmppop,params,data,state);

% ...and finally apply survival to choose individuals from both pop and tmppop,
% creating a definite new population

if params.gde
    % (if this is GDE, just use tmppop)
    pop=tmppop;
else
    [pop,state]=applysurvival(pop,params,state,tmppop);
    if params.M3GP
        % prune the first individual on the list, which is also the best
        % (applysurvival has ordered them by fitness)
        if and(~pop(1).pruned,length(pop(1).tree.kids)>1)
            pop(1)=pruneM3GPindividual(pop(1),params,data,state);
        end
    end
end

% set state measures (fitness, rank, level history):
[state,pop]=updatestate(params,state,data,pop);
