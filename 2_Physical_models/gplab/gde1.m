function newinds=gde1(pop,params,state,data,indlist)
%GDE1    Creates a new GPLAB individual by GDE operators version 1.
%   NEWIND=GDE1(POPULATION,PARAMS,STATE,DATA,PARENT) returns one new
%   individual created by GDE operators (ER uses version 1, see ref)
%   using PARENT and three random "pseudo-parents".
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
%   See also GDE2
%
%   Copyright (C) 2010-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% from Moraglio's ref:

% pick at random 3 individuals to be "pseudo-parents":
r=intrand(1,state.popsize,1,3);

W=1/(1+params.gdeF);

% create intermediate ind E by doing gdeCX with pseudo-parents 1 and 3:
tmppop(1)=pop(r(1));
tmppop(2)=pop(r(3));
params.gdeW1=1-W;
params.gdeW2=W;
E=gdeCX(tmppop,params,state,data,[1,2]);


% create mutant ind U by doing gdeER with pseudo-parent 2 and ind E:
tmppop(1)=pop(r(2));
tmppop(2)=E;
params.gdeW1=W;
params.gdeW2=1-W;
U=gdeER1(tmppop,params,state,data,[1,2]);


% create candidate offspring by doing gdeCX with ind U and the real parent:
tmppop(1)=U;
tmppop(2)=pop(indlist);
params.gdeW1=params.gdeCR;
params.gdeW2=1-params.gdeCR;
newind=gdeCX(tmppop,params,state,data,[1,2]);


newind.str=tree2str(newind.tree);      	
newind.id=[];
newind.origin='gde1';
newind.parents=[pop(indlist).id];
newind.xsites=[];
newind.fitness=[];
newind.adjustedfitness=[];
newind.result=[];
newind.testfitness=[];
newind.testadjustedfitness=[];
newind.level=[];
newind.nodes=newind.tree.nodes;
newind.introns=[];

newinds=newind;
