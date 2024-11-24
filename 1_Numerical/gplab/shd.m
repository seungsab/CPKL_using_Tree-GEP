function d=shd(tree1,tree2)
%SHD    Structural Hamming Distance between two GPLAB trees.
%   DISTANCE=SHD(TREE1,TREE2) calculates and returns the Structural Hamming Distance
%   between two trees TREE1 and TREE2.
%
%   Input arguments:
%      TREE1 - the first tree (struct)
%      TREE2 - the second tree (struct)
%   Output arguments:
%      DISTANCE - the calculated Structural Hamming Distance (double)
%
%   References:
%   Moraglio A, Silva S (2010). Geometric Differential Evolution on
%   the Space of Genetic Programs, Proceedings of EuroGP-2010, 171-183.
%
%   Copyright (C) 2010-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

k1=length(tree1.kids);
k2=length(tree2.kids);

if k1*k2==0 % both trees are terminals
    if strcmp(tree1.op,tree2.op)
        d=0;
    else
        d=1;
    end
elseif k1~=k2 % the root nodes have different arities
    d=1;
else % the trees are not terminals and the root nodes have the same arity
    if strcmp(tree1.op,tree2.op)
        hd=0; % hamming distance
    else
        hd=1;
    end
    sd_kids=0; % sum of distances measured for all the kids
    for i=1:k1
        d_kids=shd(tree1.kids{i},tree2.kids{i});
        sd_kids=sd_kids+d_kids;
    end
    d=1/(k1+1)*(hd+sd_kids);
end

    