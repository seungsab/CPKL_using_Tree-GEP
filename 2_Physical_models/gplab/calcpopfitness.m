function [pop,state]=calcpopfitness(pop,params,data,state)
%CALCPOPFITNESS    Calculate fitness values for a GPLAB population.
%   CALCPOPFITNESS(POPULATION,PARAMS,DATA,STATE) returns the
%   population with the fitness values for all individuals.
%
%   [POPULATION,STATE]=CALCPOPFITNESS(...) also returns the
%   updated state of the algorithm.
%
%   Input arguments:
%      POPULATION - the current population of individuals (array)
%      PARAMS - the running parameters of the algorithm (struct)
%      DATA - the dataset on which to measure the fitness (struct)
%      STATE - the current state of the algorithm (struct)
%   Output arguments:
%      POPULATION - the population updated with fitness (array)
%      STATE - the updated state of the algorithm (struct)
%
%   See also CALCFITNESS
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

if params.rst
    % modify data subset to use according to parameters of random sampling technique:
    state=calcRSTdata(data,state,params.rst_rss,params.rst_rss,params.rst_ri);
    data.rstdata=state.rstdata;
end

% before measuring fitness, prepare the evaluation to be assigned to the variables:
% (because the evaluation set may be variable, this should be done every generation)
for t=1:params.numvars
   % for all variables (which are first in list of inputs), ie, X1,X2,X3,...
   %state.varsvals{t}=mat2str(data.example(:,t));
   temp1=sprintf('%d;',data.example(:,t));
   temp2=['[',temp1,']'];
   temp2(end-1)=[];
   state.varsvals{t}=temp2;
   
   if params.rst
      %state.rstvarsvals{t}=mat2str(state.rstdata.example(:,t));
      temp1=sprintf('%d;',data.rstdata.example(:,t));
      temp2=['[',temp1,']'];
      temp2(end-1)=[];
      state.rstvarsvals{t}=temp2;
   end
   
   if params.usetestdata
      %state.testvarsvals{t}=mat2str(data.test.example(:,t));
      temp1=sprintf('%d;',data.test.example(:,t));
      temp2=['[',temp1,']'];
      temp2(end-1)=[];
      state.testvarsvals{t}=temp2;
   end
end

for i=1:length(pop)
   if isempty(pop(i).fitness)
      [pop(i),state2]=calcfitness(pop(i),params,data,state,0);
      % (0 = learning data, not testing)
   end
end

