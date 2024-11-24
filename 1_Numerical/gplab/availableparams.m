function [myparams,defaults]=availableparams
%AVAILABLEPARAMS    Describes the GPLAB algorithm parameter variables.
%   AVAILABLEPARAMS returns a structured variable containing the
%   possible values of each parameter variable or simply hints on
%   validating each element.
%   
%   [PARAMVARS,DEFAULTS]=AVAILABLEPARAMS also returns the default
%   initial values for each element.
%
%   Output arguments:
%      PARAMVARS - the structure of the parameters variable (struct)
%      DEFAULTS - the default initial parameters (struct)
%
%   See also AVAILABLESTATE, RESETPARAMS, SETPARAMS, RESETSTATE
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% available parameter elements:

myparams=struct(...
   'maxresources',[],...
   'functions',[],...
   'terminals',[],...
   'numvars',[],...
   'autovars',[],...
   'initpoptype',[],...
   'expected',[],...
   'elitism',[],...
   'survival',[],...
   'dynamicresources',[],...
   'resourcespopsize',[],...
   'resourcesfitness',[],...
   'periode',[],...
   'ajout',[],...
   'sampling',[],...
   'drawperspin',[],...
   'calcfitness',[],...
   'adjustfitness',[],...
   'precision',[],...
   'lowerisbetter',[],...
   'gengap',[],...
   'savetofile',[],...
   'operatorprobstype',[],...
   'initialprobstype',[],...
   'reproduction',[],...
   'operatornames',[],...
   'operatornparents',[],...
   'operatornchildren',[],...
   'initialfixedprobs',[],...
   'initialvarprobs',[],...
   'homocrosstype',[],...
   'pointmutationtype',[],...
   'broodpairs',[],...
   'gde',[],...
   'gdeF',[],...
   'gdeCR',[],...
   'shd',[],...
   'savedir',[],...
   'savename',[],...
   'datafilex',[],...
   'datafiley',[],...
   'usetestdata',[],...
   'testdatafilex',[],...
   'testdatafiley',[],...
   'files2data',[],...
   'numbackgen',[],...
   'percentback',[],...
   'percentchange',[],...
   'adaptinterval',[],...
   'adaptwindowsize',[],...
   'smalldifference',[],...
   'minprob',[],...
   'tournamentsize',[],...
   'depthnodes',[],...
   'fixedlevel',[],...
   'dynamiclevel',[],...
   'inicmaxlevel',[],...
   'inicdynlevel',[],...
   'realmaxlevel',[],...
   'veryheavy',[],...
   'calccomplexity',[],...
   'calcdiversity',[],...
   'hits',[],...
   'output',[],...
   'graphics',[],...
   'keepevalssize',[],...
   'opeq',[],...
   'opeqtype',[],...
   'opeq_binsize',[],...
   'opeq_target',[],...
   'rst',[],...
   'rst_rss',[],...
   'rst_rsr',[],...
   'rst_ri',[],...
   'M3GP',[],...
   'KNNvoting',[],...
   'filters',[]);


% possible values or some hints for validation:

myparams.maxresources='posint';
myparams.functions='list <- no need for validation';
myparams.terminals='list <- no need for validation';
myparams.numvars='posint';
myparams.autovars='boolean';
myparams.initpoptype={'growinit','fullinit','rampedinit'};
myparams.expected={'absolute','rank85','rank89'};
myparams.elitism={'keepbest','replace','halfelitism','totalelitism'};
myparams.survival={'fixedpopsize','resources','pivotfixe'};
myparams.dynamicresources={'0','1','2'}; % 0=no; 1=normal; 2=heavy
myparams.resourcespopsize={'low','steady','free'};
myparams.resourcesfitness={'normal_accept','light_accept'};
myparams.periode='posint';
myparams.ajout={'M1','M2'};
myparams.sampling={'roulette','sus','tournament','lexictour','doubletour','sampleall'};
myparams.drawperspin='posint'; % how many individuals can be drawned per wheel spin
myparams.calcfitness={'regfitness','antfitness','antfitness_lib','RMSE'};
myparams.adjustfitness={'linearppp'};
myparams.precision='posint'; % number of decimal places to use
myparams.lowerisbetter='boolean';
myparams.gengap='posnumeric';
myparams.savetofile={'never','firstlast','every10','every100','always'};
myparams.operatorprobstype={'fixed','variable'};
myparams.initialprobstype={'fixed','variable'};
myparams.reproduction='special_01float';
myparams.savedir='string <- no need for validation'; % ([] means GPLAB will ask the user)
myparams.savename='string <- no need for validation'; % ([] means GPLAB will use default name)
myparams.datafilex='string <- no need for validation'; % ([] means GPLAB will ask the user)
myparams.datafiley='string <- no need for validation'; % ([] means GPLAB will ask the user)
myparams.usetestdata='boolean';
myparams.testdatafilex='string <- no need for validation'; % ([] means GPLAB will ask the user)
myparams.testdatafiley='string <- no need for validation'; % ([] means GPLAB will ask the user)
myparams.files2data={'xy2inout','anttrail'};
myparams.operatornames='list <- no need for validation';
myparams.operatornparents='list <- no need for validation';
% crossover, mutation, shrinkmutation, swapmutation, pointmutation, brood, homocross
myparams.operatornchildren='list <- no need for validation';
myparams.initialfixedprobs='list <- no need for validation';
myparams.initialvarprobs='list <- no need for validation';
myparams.homocrosstype={'homocross','gdeCX','gdeER1','gdeER2'};
myparams.pointmutationtype={'1point','npoint'};
myparams.broodpairs='posint';
myparams.gde='boolean';
myparams.gdeF='special_posnumeric';
myparams.gdeCR='special_posnumeric';
myparams.shd='boolean';
myparams.numbackgen='posint';
myparams.percentback='special_01float'; % a float between 0 and 1, inclusive
myparams.percentchange='01float';
myparams.adaptinterval='posnumeric'; % ([] means GPLAB will set this parameter)
myparams.adaptwindowsize='posint'; % ([] means GPLAB will set this parameter)
myparams.smalldifference='01float'; % (the same)
myparams.minprob='01float';
myparams.tournamentsize='posnumeric';
myparams.depthnodes={'1','2'}; % 1=limit depth; 2=limit nodes
myparams.fixedlevel='boolean'; % use strict limit: 0=no, 1=yes
myparams.dynamiclevel={'0','1','2'}; % 0=no; 1=normal; 2=heavy
myparams.inicmaxlevel='posint';
myparams.inicdynlevel='posint';
myparams.realmaxlevel='posint';
myparams.veryheavy='boolean';
myparams.calccomplexity='boolean';
myparams.calcdiversity='list <- no need for validation';
myparams.output={'silent','normal','verbose'};
myparams.hits='matrix <- no need for validation';
myparams.graphics='list <- no need for validation';
% (available plots: plotfitness,plotdiversity,plotcomplexity,plotoperators)
myparams.keepevalssize='special_posint';
myparams.opeq='boolean';
myparams.opeqtype={'DynOpEq','MutOpEq'};
myparams.opeq_binsize='posint';
myparams.opeq_target={'none','uniform'}; % or empty empty OR uniform forget none
myparams.rst='boolean'; % random sampling technique
myparams.rst_rss='posnumeric'; % random subset size, proportion of entire set, 0 means one sample
myparams.rst_rsr='posint'; % random subset reset, number of generations with same subset
myparams.rst_ri='special_01float'; % random interleaved, probability of using entire training set
myparams.M3GP='boolean';
myparams.KNNvoting={'KNNweightedvoting'};
myparams.filters='list <- no need for validation';

% default parameters:

defaults.maxresources='[]'; % when empty, gplab will use the amount of the initial pop
defaults.initpoptype='''rampedinit''';
defaults.expected='''rank85''';
defaults.elitism='''replace''';
defaults.survival='''fixedpopsize''';
defaults.dynamicresources='''0''';
defaults.resourcespopsize='''steady''';
defaults.resourcesfitness='''normal_accept''';
defaults.periode='1';
defaults.ajout='''M1''';
defaults.sampling='''lexictour''';
defaults.drawperspin='[]'; % empty means the maximum possible value
defaults.calcfitness='''regfitness''';
defaults.adjustfitness='[]';
defaults.precision='12';
defaults.lowerisbetter='1';
defaults.savetofile='''never''';
defaults.operatorprobstype='''fixed''';
defaults.initialprobstype='''fixed''';
defaults.homocrosstype='''homocross''';
defaults.pointmutationtype='''1point''';
defaults.broodpairs='5';
defaults.gde='0';
defaults.gdeF='0.8';
defaults.gdeCR='0.9';
defaults.shd='0';
defaults.gengap='[]'; % when empty, this default will be filled in gplab 
defaults.savedir='[]'; % when empty, GPLAB will ask the user
defaults.savename='[]'; % when empty, a default name is given
defaults.datafilex='[]'; % when empty, GPLAB will ask the user
defaults.datafiley='[]'; % when empty, GPLAB will ask the user
defaults.usetestdata='0';
defaults.testdatafilex='[]'; % when empty, GPLAB will ask the user
defaults.testdatafiley='[]'; % when empty, GPLAB will ask the user
defaults.files2data='''xy2inout''';
defaults.operatornames='{''crossover'',''mutation''}';
defaults.operatornparents='[2 1]';
defaults.operatornchildren='[2 1]';
defaults.initialfixedprobs='[]'; % when empty, this default will be filled in gplab
defaults.initialvarprobs='[]'; % the same
defaults.reproduction='0.1';
defaults.numbackgen='3';
defaults.percentback='0.25';
defaults.percentchange='0.25';
defaults.adaptinterval='[]'; % when empty, this default will be filled in gplab
defaults.adaptwindowsize='[]'; % the same
defaults.smalldifference='[]'; % the same
defaults.minprob='0.1'; % no operator probability will be lower than minprob/number of operators
defaults.tournamentsize='[]'; % when empty, this default will be filled in gplab
defaults.depthnodes='''1'''; % 1=depth, 2=nodes
defaults.fixedlevel='1'; % 1=yes
defaults.dynamiclevel='''1'''; % 1=normal
defaults.inicmaxlevel='[]'; % when empty, filled automatically (consider depth or size)
defaults.inicdynlevel='[]'; % when empty, filled automatically (consider depth or size)
defaults.realmaxlevel='[]'; % when empty, filled automatically (consider depth or size)
defaults.veryheavy='0'; % 0=no
defaults.calccomplexity='0';
defaults.calcdiversity='{}';
defaults.hits='[100 0]';
defaults.output='''normal''';
defaults.graphics='{}';
% (order of plots is top-right, bottom-right, top-left, bottom-left)
defaults.functions='{''plus'' 2; ''minus'' 2; ''times'' 2; ''sin'' 1; ''cos'' 1; ''mylog'' 1}';
defaults.terminals='{}';
defaults.numvars='[]'; % when empty, it will be filled automatically depending on 'autovars'
defaults.autovars='1';
defaults.keepevalssize='[]'; % when empty, this default will be filled in gplab
defaults.opeq='0';
defaults.opeqtype='''DynOpEq''';
defaults.opeq_binsize='1'; % bin width
defaults.opeq_target='[]'; % when empty, target is calculated dynamically based on fitness
defaults.rst='0'; % by default do not use rst
defaults.rst_rss='1'; % if rst is used, by default use only one training sample
defaults.rst_rsr='1'; % if rst is used, by default change sample every generation
defaults.rst_ri='0.5'; % if rst is used, each generation probability of 0.5 of using entire training set
defaults.M3GP='0';
defaults.KNNvoting='''KNNweightedvoting''';
defaults.filters='{}';
% (filters are to be filled automatically depending on 'depthnodes', 'fixedlevel' and 'dynamiclevel')
% (possible filters: strictdepth, strictnodes, dyndepth, dynnodes, heavydyndepth, heavydynnodes)
% (additional filter for operator equalisation: opeq)

% each operator name must be available as a matlab function.
% the user chooses whatever functions from the list of available functions:
% any MATLAB functions that verify closure may be used.

%defaults.functions='{
   %'ceil' 1;
   %'floor' 1;
   %'times' 2;
   %'mydivide' 2;
   %'plus' 2;
   %'minus' 2;
   %'uminus' 1;
   %'mysqrt' 1;
   %'mylog2' 1;
   %'mylog10' 1;
   %'mylog' 1;
   %'mypower' 2;
   %'sin' 1;
   %'cos' 1;
   %'min' 2;
   %'max' 2;
   %'abs' 1;
   %'myif' 3;
   %'eq' 2;
   %'ne' 2;
   %'lt' 2;
   %'gt' 2;
   %'le' 2;
   %'ge' 2;
   %'and' 2;
   %'or' 2;
   %'not' 1;
   %'xor' 2
%}';

