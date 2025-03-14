function [vars,best]=gplab(g,varargin)
%GPLAB    Runs the GPLAB genetic programming algorithm.
%   GPLAB(NGENS,POPSIZE) initializes a GPLAB algorithm using
%   POPSIZE individuals, and runs it for NGENS generations
%   using default parameter variables. If NGENS=0 only the
%   initializations are done. It returns the algorithm
%   variables after the run.
%
%   GPLAB(NGENS,POPSIZE,PARAMS) uses previously set algorithm
%   parameters, PARAMS, instead of the default ones.
%
%   GPLAB(NGENS,VARS) continues a GPLAB run for NGENS generations,
%   starting from the point defined by the algorithm variables VARS.
%
%   [VARS,BEST] = GPLAB(...) also returns the best individual found
%   during the run, which is already part of the algorithm variables.
%
%   Input arguments:
%      NGENS - the number of generations to run the algorithm (integer)
%      POPSIZE - the number of individuals in the population (integer)
%      PARAMS - the algorithm running parameters (struct)
%      VARS - the algorithm variables (struct)
%        VARS.POP - the current population
%        VARS.PARAMS - the algorithm running parameters = PARAMS
%        VARS.STATE - the current state of the algorithm
%   Output arguments:
%      VARS - the algorithm variables (struct) - see Input arguments
%      BEST - the best individual found in the run (struct)
%
%   See also SETPARAMS, RESETPARAMS, RESETSTATE
%
%   --------------------------
%   See demo functions DEMO*
%   --------------------------
%
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

if (nargin<2) || (nargin>3)
   error('GPLAB: Wrong number of input arguments. Use either gplab(ngens,vars) to continue a run, or gplab(ngens,popsize,[optional params]) to start a run')   
   
elseif isstruct(varargin{1})
   % argument 1: the number of additional generations to run
   % argument 2: the algorithm variables
   if ~(isvalid(g,'posint'))
      error('GPLAB: The first argument must be an integer greater than 0.')
   end
   start=0;
   continuing=1;
   vars=varargin{1};
   n=vars.state.popsize;
   level=vars.state.maxlevel;
   ginic=vars.state.generation+1; % start generation number
   gend=ginic-1+g; % end generation number
   
else
   % argument 1: the number of generations to run
   % argument 2: the number of individuals in the population
   % argument 3: (optional) the parameters of the algorithm
   if ~(isvalid(g,'special_posint') && isvalid(varargin{1},'posint') && varargin{1}>=2)
      error('GPLAB: The first two arguments must be integers, and the second > 1')
   end
   start=1;
   continuing=0;
   n=varargin{1};
   if nargin==3
      vars.params=varargin{2};
   else
      vars.params=[];
   end
   vars.state=[];
   vars.data=[];
   ginic=1; % start generation number
   gend=g; % end generation number
end

% check parameter variables:
vars.params=checkvarsparams(start,continuing,vars.params,n);

% check data variables:
[vars.data,vars.params]=checkvarsdata(start,continuing,vars.data,vars.params);

% check state variables:
[vars.state,vars.params]=checkvarsstate(start,continuing,vars.data,vars.params,vars.state,n,g);

% initialize random number generator (see help on RAND):
rand('state',sum(100*clock));

displaystatus('\nRunning algorithm...\n');

% initiate graphics:
% (if we're not going to run generations or draw history, don't initiate the graphics)
if ~isempty(vars.params.graphics) && (ginic<=gend || continuing) 
   gfxState=graphicsinit(vars.params);
end


% initial generation:

if start
   [vars.pop,vars.state]=genpop(vars.params,vars.state,vars.data,n);
   if strcmp(vars.params.savetofile,'firstlast') || strcmp(vars.params.savetofile,'every10') || strcmp(vars.params.savetofile,'every100') || strcmp(vars.params.savetofile,'always')
      saveall(vars);
   end
   if ~strcmp(vars.params.output,'silent')
      displaystats('     #Individuals:  %d\n',vars.state.popsize,vars.state.generation);
      if strcmp(vars.params.survival,'resources')
	displaystats('     MaxResources:  %d\n',vars.state.maxresources,vars.state.generation);
      end
      displaystats('     UsedResources: %d\n',vars.state.usedresources,vars.state.generation);
      if vars.params.opeq
        displaystats('     Highest Bin:   %d\n',length(vars.state.opeq.bestfitness),vars.state.generation);
      end
      displaystats('     Best so far:   %d\n',vars.state.bestsofar.id,vars.state.generation);
      displaystats('     Fitness:       %f\n',vars.state.bestsofar.fitness,vars.state.generation);
      if vars.params.usetestdata
         displaystats('     Test fitness:  %f\n',vars.state.bestsofar.testfitness,vars.state.generation);
      end
      displaystats('     Depth:         %d\n',vars.state.bestsofar.level,vars.state.generation);
      displaystats('     Nodes:         %d\n',vars.state.bestsofar.nodes,vars.state.generation);
      if vars.params.M3GP
          displaystats('     Dimensions:    %d\n',vars.state.bestsofar.dimensions,vars.state.generation);
      end
      displaystatus('\n');
   end
   % (if we're not going to run generations, don't start the graphics:)
   if ~isempty(vars.params.graphics) && ginic<=gend
      gfxState=graphicsstart(vars.params,vars.state,gfxState);
   end
end

if continuing
   if ~isempty(vars.params.graphics)
      gfxState=graphicscontinue(vars.params,vars.state,gfxState);
   end
end

sc=0;

% generations:
   
for i=ginic:gend
   
   % stop condition?
   sc=stopcondition(vars.params,vars.state,vars.data);
   if sc
      % unless the option is to never save, save the algorithm variables now:
      if (~strcmp(vars.params.savetofile,'never'))
         saveall(vars);
      end
      break % if a stop condition has been reached, skip the for cycle
   end
   

   [vars.pop,vars.state]=generation(vars.pop,vars.params,vars.state,vars.data,vars.params.gengap);
   
   % save to file?
   if (strcmp(vars.params.savetofile,'firstlast') && i==g) || (strcmp(vars.params.savetofile,'every10') && rem(i,10)==0) || (strcmp(vars.params.savetofile,'every100') && rem(i,100)==0) || strcmp(vars.params.savetofile,'always')
      saveall(vars);
   end
   
   % textual output:
   if ~strcmp(vars.params.output,'silent')
      displaystats('     #Individuals:  %d\n',vars.state.popsize,vars.state.generation);
      if strcmp(vars.params.survival,'resources')
	displaystats('     MaxResources:  %d\n',vars.state.maxresources,vars.state.generation);
      end
      displaystats('     UsedResources: %d\n',vars.state.usedresources,vars.state.generation);
      if vars.params.opeq
        displaystats('     Highest Bin:   %d\n',length(vars.state.opeq.bestfitness),vars.state.generation);
      end
      displaystats('     Best so far:   %d\n',vars.state.bestsofar.id,vars.state.generation);
      displaystats('     Fitness:       %f\n',vars.state.bestsofar.fitness,vars.state.generation);
      if vars.params.usetestdata
         displaystats('     Test fitness:  %f\n',vars.state.bestsofar.testfitness,vars.state.generation);
      end
      displaystats('     Depth:         %d\n',vars.state.bestsofar.level,vars.state.generation);
      displaystats('     Nodes:         %d\n',vars.state.bestsofar.nodes,vars.state.generation);
      if vars.params.M3GP
          displaystats('     Dimensions:    %d\n',vars.state.bestsofar.dimensions,vars.state.generation);
      end
      displaystatus('\n');
   end
   
   % plots:
   if ~isempty(vars.params.graphics)
      gfxState=graphicsgenerations(vars.params,vars.state,gfxState);
   end 
      
end % for i=ginic:gend


% messages regarding the stop condition reached:

if sc
   if vars.state.generation==0
      displaystatus('\nStop condition #%d was reached after initial generation.\n',sc);      
   else
      displaystatus('\nStop condition #%d was reached after generation %d.\n',sc,vars.state.generation);
   end
else
   displaystatus('\nMaximum generation %d was reached.\n',vars.state.generation);
end      

best=vars.state.bestsofar;
vars.state.keepevals=[]; % clear memory, we don't want to save all this!

displaystatus('\nDone!\n\n');
