function newinds=gdeER1(pop,params,state,data,indlist)
%GDEER1    Creates a new GPLAB individual by GDE extension ray version 1.
%   NEWIND=GDEER1(POPULATION,PARAMS,STATE,DATA,PARENT) returns one new
%   individual created by GDE extension ray version 1 of two parents.
%   (GDE = Geometric Differential Evolution)
%
%   Input arguments:
%      POPULATION - the population where the parents are (array)
%      PARAMS - the parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      DATA - the dataset on which the algorithm runs (struct)
%      PARENT - the index of the parent in POPULATION (integer)
%   Output arguments:
%      NEWIND - the newly created individual (struct)
%
%   References:
%   Moraglio A, Silva S (2010). Geometric Differential Evolution on
%   the Space of Genetic Programs, Proceedings of EuroGP-2010, 171-183.
%
%   See also GDEER,GDEER2
%
%   Copyright (C) 2010-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

params.homocrosstype='gdeER1';
newinds=gdeER(pop,params,state,data,indlist);
