function [v,b]=demoM3GP
%DEMOM3GP    Demonstration function of the GPLAB toolbox.
%   DEMOM3GP runs a classification problem with the M3GP algorithm.
%
%   VARS=DEMOM3GP returns all the variables of the algorithm
%   needed to plot charts, continue runs, draw trees, etc.
%
%   [VARS,BEST]=DEMOM3GP also returns the best individual found
%   during the run (the same as 'vars.state.bestsofar').
%
%   See also DEMO,DEMOPARITY,DEMOANT,DEMOPLEXER
%
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox


fprintf('Running M3GP classification demo...');
p=resetparams;

p.M3GP=1;
p.calcfitness='M3GPaccuracy';
%p.calcfitness='M3GPaccuracy_1class_nclusters';

%p.calcfitness='KNNaccuracy';
%p.KNNvoting='KNNweightedvoting';
%p.KNNvoting='KNNbasicvoting';

%p.calcfitness='DTaccuracy';

p.lowerisbetter=0;


p.rst=0;
p.rst_rss=0.5; % use 50% of samples each time
p.rst_rsr=1; % change every generation
p.rst_ri=0; % probability of choosing to use the entire data set is null

p=setoperators(p,'crossover',2,2,'mutation',1,1);
p.operatorprobstype='fixed';
p.minprob=0;

p=setfunctions(p,'plus',2,'minus',2,'times',2,'kozadivide',2);
p=setterminals(p,'rand');

p.datafilex='heart_x.txt';
p.datafiley='heart_y.txt';

p.usetestdata=0;
%p.testdatafilex='gee_coords2_x.txt';
%p.testdatafiley='gee_coords2_y.txt';

p.calcdiversity={};
p.calccomplexity=0;
%p.graphics={'plotfitness','plotdiversity','plotcomplexity','plotoperators'};
p.graphics={};
p.fixedlevel=0;
p.dynamiclevel=0;

%p.savetofile='every10';
%p.savedir='C:\Users\SaraSilva\Desktop\GPLAB\GP_Data_test';


[v,b]=gplab(50,100,p);


drawtree(b.tree);