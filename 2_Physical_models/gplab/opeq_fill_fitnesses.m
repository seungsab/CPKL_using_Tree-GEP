function [fitnesses,current]=opeq_fill_fitnesses(pop,binsize);
%OPEQ_FILL_FITNESSES    Store fitnesses of GPLAB individuals for each OpEq bin.
%   [FITNESSES,CURRENT]=OPEQ_FILL_FITNESSES(POP,BINSIZE) returns the fitnesses
%   of all the individuals in the current population (POP) divided by bin,
%   taking into consideration the bin width (BINSIZE). Returns also the
%   occupancy (number of individuals) of each bin.
%
%   Input arguments:
%      POP - the current population (array)
%      BINSIZE - width of the bins (integer)
%   Output arguments:
%      FITNESSES - fitnesses per bin (cell array)
%      CURRENT - occupancy per bin (array)
%
%   References:
%      Silva S, Vanneschi L, "Operator Equalisation, Bloat and Overfitting 
%      - A Study on Human Oral Bioavailability Prediction", GECCO-2009.
%
%   See also OPEQ_GETBIN, OPEQ_CALCTARGET
%
%   Copyright (C) 2008-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

fitnesses={};
current=[];
for i=1:length(pop)
    % make sure fitness and nodes fields are filled
    bin=opeq_getbin(pop(i).nodes,binsize);
    if length(fitnesses)<bin
        fitnesses{bin}=[];
        current(bin)=0;
    end
    fitnesses{bin}(end+1)=pop(i).fitness;
    current(bin)=current(bin)+1;
end