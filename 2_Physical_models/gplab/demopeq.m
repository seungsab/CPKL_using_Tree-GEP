function [v,b]=demopeq
%DEMOPEQ    Demonstration function of the GPLAB toolbox
%
%   See also DEMO,DEMOPARITY,DEMOPLEXER,DEMOANT
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

fprintf('Running operator equalisation demo...');
p=resetparams;
p.opeq=1;
p.opeqtype='MutOpEq';
p.opeq_binsize=1;
%p.opeq_target='none';
%p.opeq_target='uniform';
p.opeq_target=[]; % dynamic

p=setoperators(p,'crossover',2,2);
p.operatorprobstype='fixed';
p.inicmaxlevel = 6;
p.initpoptype = 'rampedinit';
p.sampling = 'lexictour';
p.tournamentsize=10;

p.calcfitness='RMSE';

p = setfunctions(p, 'plus', 2, 'minus', 2, 'times', 2, 'mydivide', 2);

p.datafilex='dataBio_x.txt';
p.datafiley='dataBio_y.txt';

p.calccomplexity=0;
p.calcdiversity={};
p.graphics={};

p.fixedlevel=0;
%p.dynamiclevel='1';
p.dynamiclevel='0';

[v,b]=gplab(100,200,p);
