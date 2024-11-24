function [indid,indindex,popexpected,popnormfitness]=sampleall(pop,params,state,nsample,toavoid)
%TOURNAMENT    Sampling of GPLAB individuals by the tournament method.
%   TOURNAMENT(POP,PARAMS,STATE,NSAMPLE,TOAVOID) returns NSAMPLE ids of
%   individuals chosen from POP using the tournament method. The ids in
%   TOAVOID are not chosen.
%
%   [IDS,INDICES]=TOURNAMENT(POP,PARAMS,STATE,NSAMPLE,TOAVOID) also
%   returns the indices in POP of the chosen individuals.
%
%   Input arguments:
%      POPULATION - the current population of the algorithm (array)
%      PARAMS - the running parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      NSAMPLE - the number of individuals to draw (integer)
%      TOAVOID - the ids of the individuals to avoid drawing (1xN matrix)
%   Output arguments:
%      IDS - the ids of the individuals chosen (1xN matrix)
%      INDICES - the indices of the individuals chosen (1xN matrix)
%
%   Note:
%      The two last output arguments are not referred because they are
%      not used in this function. They are present only for compatibility
%      with the other functions for sampling individuals.
%
%   See also ROULETTE, SUS, LEXICTOUR, DOUBLETOUR, TOURBEST, SAMPLING
%
%   Copyright (C) 2003-2007 Sara Silva (sara@dei.uc.pt)
%   This file is part of the GPLAB Toolbox

% we are not calculating these, but they are output arguments:
popexpected=[];
popnormfitness=[];

state.popexpected=ones(1,state.popsize);

indindex=1:state.popsize;
indid=[pop.id];

