% INITIALIZE GEP
p=resetparams;

%% Set Training data
% This is a trick to define the kernel basis as inputs of GE programming
% These values are generated artificially and they are not meaningful
n_kernel_basis = 8;
tempx = randn(10,n_kernel_basis);
tempy = randn(10,1);

save([pwd '\gplab\kerenlbasisx.txt'],'tempx','-ascii'); % type tempx.txt
save([pwd '\gplab\kerenlbasisy.txt'],'tempy','-ascii'); % type tempy.txt
p.datafilex='kerenlbasisx.txt';
p.datafiley='kerenlbasisy.txt';

%% Tree initialization
p.depthnodes='1';
p.initpoptype = 'growinit';

%% Set tree depth and size limits
p.dynamiclevel = '1';

%% Operator and terminals
p = setfunctions(p, 'plus', 2, 'times', 2);

%% Define fitness function
p.calcfitness = 'GP_Structure_Discovery_GEP';
