function newinds=gdeCX(pop,params,state,data,indlist)
%GDECX    Creates a new GPLAB individual by GDE crossover.
%   NEWINDS=GDECX(POPULATION,PARAMS,STATE,DATA,PARENTS) returns one new
%   individual created by swapping node functions or entire subtrees
%   of the two PARENTS at homologous points. At each point, each parent
%   has a certain probability of contributing genetic material.
%   (GDE = Geometric Differential Evolution)
%
%   Input arguments:
%      POPULATION - the population where the parents are (array)
%      PARAMS - the parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      DATA - the current dataset for the algorithm to run (array)
%      PARENTS - the indices of the parents in POPULATION (1x2 matrix)
%   Output arguments:
%      NEWIND - the newly created individual (1x2 matrix)
%
%   References:
%   Moraglio A, Silva S (2010). Geometric Differential Evolution on
%   the Space of Genetic Programs, Proceedings of EuroGP-2010, 171-183.
%
%   See also HOMOCROSS, CROSSOVER, MUTATION, APPLYOPERATOR
%
%   Copyright (C) 2010-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

params.homocrosstype='gdeCX';
newinds=homocross(pop,params,state,data,indlist);
newinds=newinds(1); % ignore second child
