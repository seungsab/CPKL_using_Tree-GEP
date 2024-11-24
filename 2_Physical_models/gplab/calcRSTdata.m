function state=calcRSTdata(data,state,rst_rss,rst_rsr,rst_ri)
%CALCRSTDATA    Determines the RST data subset to be used in this GPLAB generation.
%   CALCRSTDATA(DATA,STATE,RST_RSS,RST_RSR,RST_RI) returns a
%   subset of the training data according to the RST parameters.
%
%   Input arguments:
%      DATA - datasets being used in this run (struct)
%      STATE - the current state of the algorithm (struct)
%      RST_RSS - random subset size: number of samples or proportion of the entire training set (double)
%      RST_RSR - random subset reset: how many generations before changing subset (integer)
%      RST_RI - random interleaved : probability of using entire set in each generation (double) 
%   Output arguments:
%      STATE - updated state of the algorithm with data for RST (struct)
%
%   References:
%   - Goncalves I, Silva S, Melo JB, Carreiras JMB. Random Sampling Technique
%   for Overfitting Control in Genetic Programming. EuroGP-2012.
%   - Goncalves I, Silva S. Balancing Learning and Overfitting in Genetic
%   Programming with Interleaved Sampling of Training Data. EuroGP-2013.
%
%   Copyright (C) 2014 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

state.rstdata.classes=data.classes;

if state.generation==0 % first time, calculate ndatarows (size of subset)
    if rst_rss>=1
        state.rstdata.ndatarows=min(round(rst_rss),length(data.result));
    else
        state.rstdata.ndatarows=min(round(rst_rss*length(data.result)),length(data.result));
    end
end

r=rand;
if or(r>rst_ri,mod(state.generation,rst_rsr)==0)
    % r>rst_ri : probability of reshuffling is 1-rst_ri (1 - probability of interleaving)
    % mod(gen,rst_rsr)==0 : time to reshuffle, also catches gen=0
	myshuffleddata=shuffle([data.example data.result],1);
	state.rstdata.example=myshuffleddata(1:state.rstdata.ndatarows,1:end-1);
	state.rstdata.result=myshuffleddata(1:state.rstdata.ndatarows,end);
end
