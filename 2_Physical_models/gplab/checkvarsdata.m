function [data,params]=checkvarsdata(start,continuing,data,params);
%CHECKVARSDATA    Fills the dataset variable for the GPLAB algorithm.
%   CHECKVARSDATA(START,CONTINUE,DATA,PARAMS) returns the dataset
%   on which the GPLAB algorithm will run, after prompting the user
%   for the names of the files containing this data.
%
%   [DATA,PARAMS]=CHECKVARSDATA(START,CONTINUE,PARAMS) also
%   returns the updated algorithm parameters.
%
%   Input arguments:
%      START - true if no generations have been run yet (boolean)
%      CONTINUE - true if some generations have been run (boolean)
%      DATA - the current dataset for the algorithm to run (array)
%      PARAMS - the algorithm running parameters (struct)
%   Output arguments:
%      DATA - the dataset on which the algorithm will run (struct)
%      PARAMS - the algorithm running parameters (struct)
%
%   See also CHECKVARSSTATE, CHECKVARSPARAMS, GPLAB
%
%   Copyright (C) 2018 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% check training data variables:

if start
   
   % ask for the name of the files:
   parentdir=strrep(which(mfilename),[mfilename,'.m'],'');
   
   % training data:
   if isempty(params.datafilex) % if not in params, ask user
   	filenamex=input('Please name the file (with extension) containing the input data: ','s');
   else
      filenamex=params.datafilex;
   end
   if isempty(findstr(filenamex,filesep))
	   % if file was not given with path, use parentdir (where this file is):
      filenamex=[parentdir,filenamex];
   end
   params.datafilex=filenamex;
   x=load(params.datafilex); % load the file
   
	if isempty(params.datafiley) % if not in params, ask user
   	filenamey=input('Please name the file (with extension) containing the desired output: ','s');
   else
      filenamey=params.datafiley;
   end
   if isempty(findstr(filenamey,filesep))
      % if file was not given with path, use parentdir (where this file is):
      filenamey=[parentdir,filenamey];
   end
   params.datafiley=filenamey;
   y=load(params.datafiley); % load the file
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % MODIFIED
   testfilenamex=params.testdatafilex;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % test data:
   if params.usetestdata
      
		if isempty(params.testdatafilex)
   		testfilenamex=input('Please name the file (with extension) containing the test input data: ','s');
      else
      	testfilenamex=params.testdatafilex;
      end
      if isempty(findstr(testfilenamex,filesep))
         testfilenamex=[parentdir,testfilenamex];
      end
      params.testdatafilex=testfilenamex;
      testx=load(params.testdatafilex); % load the file
      
		if isempty(params.testdatafiley)
   		testfilenamey=input('Please name the file (with extension) containing the test expected output: ','s');
   	else
      	testfilenamey=params.testdatafiley;
      end
      if isempty(findstr(testfilenamey,filesep))
	      testfilenamey=[parentdir,testfilenamey];
      end
      params.testdatafiley=testfilenamey;
      testy=load(params.testdatafiley); % load the file
      
   end %if params.usetestdata
   
	% create variable for the algorithm:
   data=feval(params.files2data,x,y);
   if params.usetestdata   	
	   data.test=feval(params.files2data,testx,testy);
   end
   
   % examples for training are now vectors containing the input data:
   %      data.example(1,:)=[0.75,0.25,0.25,0.25];
   %      data.result(1)=[0,0];
   %      data.example(2,:)=[0.875,0.625,0.625,0.625];
   %      data.result(2)=[1,2];
   % test data is something like:
   %      data.test.example(1,:)=[0.85,0.5,0.5,0.75];
   %      data.test.result(1)=[0,1];
   %      data.test.example(2,:)=[0.5,0.25,0.25,0.0];
   %      data.test.result(2)=[1,3];
   
   % save the terminals, for efficiency:
   % ***
   
   if params.M3GP
       % assume data.result is unidimensional;
       % also assume that data.result contains all possible classes
       % (no new classes will appear in data.test.result)
       data.classes=unique(data.result,'stable'); % no sorting!
       %[~,~,data.classfrequencies]=countfind(data.result,data.classes);
        if params.usetestdata   	
            %data.test.classes=unique(data.test.result,'stable');
            % if something in test was not in training, issue an error
            % *******TO DO!!!!!!
            data.test.classes=data.classes;
        %    [~,~,data.test.classfrequencies]=countfind(data.test.result,data.test.classes);
        end
   end
   
end
