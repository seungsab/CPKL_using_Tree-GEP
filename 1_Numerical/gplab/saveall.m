function saveall(vars)
%SAVEALL    Saves all the GPLAB algorithm variables to disk.
%   SAVEALL(VARIABLES) saves all the data stored in VARIABLES
%   into a file identified with the number of the current
%   generation, in the directory indicated in the algorithm
%   parameters.
%
%   Input arguments:
%      VARIABLES - contains: POP, PARAMS, STATE, DATA (struct)
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   Acknowledgements: Matthew Clifton (matthew_clifton@hotmail.com)
%   This file is part of the GPLAB Toolbox

%   Modified May 2007:
%   - added '.mat' to file name for compatibility with Octave

if isempty(vars.params.savename)
    filename=[vars.params.savedir '/' num2str(vars.state.generation) '.mat'];
else
    filename=[vars.params.savedir '/' vars.params.savename num2str(vars.state.generation) '.mat'];
end
save(filename,'vars');
