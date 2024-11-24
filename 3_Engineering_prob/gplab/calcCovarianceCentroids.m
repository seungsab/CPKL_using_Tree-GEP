function [covM,invcovM,centroids]=calcCovarianceCentroids_1class_nclusters(mapping,data);

nclasses=length(data.classes);
[~,nfeatures]=size(mapping);
centroids=zeros(nfeatures,nclasses);
for c=1:nclasses
    covM{c}=zeros(nfeatures,nfeatures);
    invcovM{c}=zeros(nfeatures,nfeatures);
end


for c=1:nclasses
    % get from mapping only the observations of the c-th class:
    cmapping=mapping(data.result==data.classes(c),:);
    % calculate covariance and centroids:
    covM{c}=cov(cmapping);
    centroids(:,c)=mean(cmapping);
    % calculate inverse of covM, if possible (if not, remains zeros):
    [~,pd] = chol(covM{c});
    if and(rcond(covM{c})>1e-10,not(pd))
        % if the covariance matrix is not "nearly singular"
        % AND is "positive semidefinite" (pd==0), invert it
        invcovM{c}=inv(covM{c});
    else
        % otherwise make the inverse equal to the identity matrix
        invcovM{c}=eye(nfeatures);
    end
end