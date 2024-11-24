function target_occupancy=opeq_calctarget(popsize,target,fitnesses,bestfitness,lowerisbetter);
%OPEQ_CALCTARGET    Calculates the OpEq target for the next GPLAB population.
%   TARGET_OCCUPANCY=OPEQ_CALCTARGET(POPSIZE,TARGET,FITNESSES,BESTS,LBETTER)
%   returns the target occupancy for the next generation. It does so by
%   transforming the FITNESSES distribution into a proportional distribution
%   of individuals, thus defining each bin's capacity.
%
%   Input arguments:
%      POPSIZE - the size of the current population (integer)
%      TARGET - type of target (string)
%      FITNESSES - fitnesses of the individuals in each bin (cell array)
%      BESTS - best fitness of each bin (array)
%      LBETTER - lower fitness is better (boolean)
%   Output arguments:
%      TARGET_OCCUPANCY - target occupancy per bin (array)
%
%   References:
%      Silva S, Vanneschi L, "Operator Equalisation, Bloat and Overfitting 
%      - A Study on Human Oral Bioavailability Prediction", GECCO-2009.
%
%   See also OPEQ_GETBIN, OPEQ_FILL_FITNESSES
%
%   Copyright (C) 2008-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% just to test a uniform distribution:
if strcmp(target,'uniform')
    target_occupancy(1:length(fitnesses))=round(popsize/length(fitnesses));
else

% this code should be vectorized somehow...

% calculate mean (or best) fitness for each bin (ignore empty bins):
% current implementation: mean. it allows a better differenciation in
% problems where possible fitness values are few (?)
m=[]; s=[];
for i=1:length(fitnesses)
    if ~isempty(fitnesses{i})
        m(end+1)=mean(fitnesses{i});
        %m(end+1)=bestfitness(i);
        s(end+1)=length(fitnesses{i});
    end
end

% normalize vector of mean/best fitnesses so its length is popsize

% remove infinite and nan values (give them the worst fitness):
m(abs(m)==Inf)=max(m(abs(m)~=Inf));
m(isnan(abs(m)))=max(m(~isnan(abs(m))));

mm=m;
if lowerisbetter
    m=normalize(-m,popsize);
else
    m=normalize(m,popsize);
end
if isnan(m(1))
    mm
end
% round occupancy to integer values:
t=round(m);

% set occupancy for all bins (zero for empty bins):
r=1;
max_m=0;
for i=1:length(fitnesses)
    if ~isempty(fitnesses{i})
        target_occupancy(i)=t(r);
        if m(r)>max_m
            bestbin=i;
            max_m=m(r);
        end
        r=r+1;
    else
        target_occupancy(i)=0;
    end
end

end