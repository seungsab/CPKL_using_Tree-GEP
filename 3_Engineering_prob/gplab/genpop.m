function [pop,state]=genpop(params,state,data,n,level)
%GENPOP    Creates the initial generation for the GPLAB algorithm.
%   GENPOP(PARAMETERS,STATE,DATA,POPSIZE,MAXLEVEL) returns a
%   population of POPSIZE new individuals for the GPLAB algorithm,
%   with tree representations no deeper than MAXLEVEL, already
%   containing the fitness values as measured in DATA.
%
%   [POPULATION,STATE]=GENPOP(PARAMETERS,STATE,DATA,POPSIZE,MAXLEVEL)
%   also returns the updated state after creating the population,
%   containing several fitness measures and, depending on the
%   algorithm PARAMETERS, data for the procedure of automatic operator
%   probabilities setting (Davis 89).
%
%   Input arguments:
%      PARAMETERS - the algorithm running parameters (struct)
%      STATE - the algorithm current state (struct)
%      DATA - the dataset on which the algorithm runs (struct)
%      POPSIZE - the number of individuals to create (integer)
%      MAXLEVEL - the maximum depth of the new individuals (integer)
%   Output arguments:
%      POPULATION - the population of new individuals (array)
%      STATE - the algorithm current state, now updated (struct)
%
%   See also INITPOP, GENERATION
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

if ~strcmp(params.output,'silent')
   displaystatus('\n- Creating initial generation -\n');
end

% create a new population:
if params.M3GP
    [pop,state.lastid]=initpop(n,state.lastid,params.inicmaxlevel,state.functions,state.arity,params.initpoptype,params.depthnodes,1);
    % the last parameters is the number of dimensions of the initial
    % individuals. by default, the initial individuals have only one
    % dimension. If not M3GP, dimensions=[]
else
    [pop,state.lastid]=initpop(n,state.lastid,params.inicmaxlevel,state.functions,state.arity,params.initpoptype,params.depthnodes,[]);
end

% if variable probabilities, update adaptwindow and op history:
if strcmp(params.operatorprobstype,'variable')
   state=automaticoperatorprobs(pop,[],params,state,data,0,n);
   % (the second input parameter should be the population of
   %  the previous generation, but it still does not exist)
end

[state,pop]=updatestate(params,state,data,pop);

% if "populations dynamiques", set pivot:
if strcmp(params.survival,'pivotfixe')
    state.pivot=state.maxfitness/state.maxgen;
end

