function complexityvalues=complexity(xdata,ydata)

% error checking:
nrowsx=size(xdata,1);
nrowsy=size(ydata,1);
if nrowsx~=nrowsy
    error('COMPLEXITY: The number of rows must be the same in both matrices.')
elseif nrowsx<2
    error('COMPLEXITY: The number of rows in each matrix must be at least 2.')
end
    
% number of dimensions of input data:
ndims=size(xdata,2);
% number of points in each dimension:
npoints=nrowsx; % = size(xdata,1);

% initialize vector of complexity values (one value per dimension):
complexityvalues=zeros(1,ndims);

% sort the x values on each dimension (by ascending order):
[sortedxdata,i]=sort(xdata,1);
% sort the y values (using the ndims different orders):
ydata=repmat(ydata,1,ndims); % we need one y vector for each dimension
sortedydata=ydata(i);

% calculate complexity values for each dimension separately:
for d=1:ndims

    % because the same x value may have different y values (when ndims>1)
    % we use the median of the y values. to do this:
    
    % find the indexes with repeated x values:
    x=sortedxdata(:,d);
    [ux,i]=unique(x,'legacy'); % i returns the indexes of the unique values
    ri=setdiff(1:npoints,i); % ri returns the indexes of the repetitions
    % (when there is a repeated value, i contains the index of the *last*
    %  occurrence of the value, while ri contains the indexes of the n
    %  first occurrences of the value (n>=1)
    
    if ~isempty(ri) % if there are no repetitions skip this entire part
        
        % count occurrences, to help identify the "groups" of repetitions:
        [ans,c]=countfind(x',ux'); % how many times each value appears
        c=c(c~=1); % we are only interested if it appears more than once
        % (c now contains the number of occurrences of each repeated value)
        % (c and ri can now be used to transform the y vector to replace
        %  repetitions with median values)
    
        % build final y vector with medians where necessary:
        y=sortedydata(:,d);
        yfinal=[];
        endlastrep=0; % index on the y vector, end of last repetition
        for i=1:length(c)
            beginthisrep=ri(sum(c(1:i-1))-(i-1)+1);
            yfinal=[yfinal;y(endlastrep+1:beginthisrep-1)];
            % (appends to yfinal everything until this repetition)
            endthisrep=beginthisrep+c(i)-1;
            % (take note where this repitition ends)
            yfinal(end+1,1)=median(y(beginthisrep:endthisrep));
            % (appends to yfinal the median value)
            % (do NOT remove the ,1 or you may get a row instead of column)
            endlastrep=endthisrep;
        end
        yfinal=[yfinal;y(endthisrep+1:end)];
        % (appends to yfinal everything after the final repetition)
    
        % also build final x vector:
        xfinal=ux; % no secret here, just take the unique values
    
    else
        yfinal=sortedydata(:,d);
        xfinal=x;
    end  % if ~isempty(ri)
    
    % finally, calculate slopes using xfinal and yfinal:
    xfinal1=xfinal(1:end-1);
    xfinal2=xfinal(2:end);
    yfinal1=yfinal(1:end-1);
    yfinal2=yfinal(2:end);
    xfinal0=xfinal2-xfinal1; % differences between consecutive x points
    yfinal0=yfinal2-yfinal1; % differences between consecutive y points
    slopes=yfinal0./xfinal0;
    
    % absolute differences between consecutive slopes:
    slopes1=slopes(1:end-1);
    slopes2=slopes(2:end);
    slopes0=abs(slopes2-slopes1);
    
    % complexity of this dimension is the sum of the differences:
    cvalue=sum(slopes0);
    complexityvalues(d)=cvalue;
end

% - Now we use whatever we want from this vector
% - I suggest using (and comparing results with) average and maximum values


