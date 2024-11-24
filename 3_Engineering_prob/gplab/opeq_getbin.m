function bin=opeq_getbin(size,binsize)
%OPEQ_GETBIN    Calculates which OpEq bin the GPLAB individual fits.
%   BIN=OPEQ_GETBIN(SIZE,BINSIZE) returns the number of the bin where
%   an individual of size SIZE fits, considering the bin width is BINSIZE.
%
%   Input arguments:
%      SIZE - number of nodes of the individual (integer)
%      BINSIZE - width of the bins (integer)
%   Output arguments:
%      BIN - number of the bin where the individual fits (integer)
%
%   References:
%      Silva S, Vanneschi L, "Operator Equalisation, Bloat and Overfitting 
%      - A Study on Human Oral Bioavailability Prediction", GECCO-2009.
%
%   See also OPEQ_FILL_FITNESSES, OPEQ_CALCTARGET
%
%   Copyright (C) 2008-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

bin=floor((size-1)/binsize)+1;
