function ret=evaluate_tree(tree,X)
%EVALUATE_TREE    Alternative to using 'eval' in 'regfitness' or 'RMSE' in GPLAB.
%   EVALUATE_TREE(TREE,X) evaluates a GPLAB tree when the "nesting 32"
%   error is issued by the 'eval' function in the 'regfitness.m' function.
%   It returns the same as 'eval' would. This is not needed in older 
%   MATLAB releases.
%
%   See also REGFITNESS, RMSE
%
%   Created (2007) by Marco Medori (marco.medori@poste.it) and Bruno Morelli
%   Modified (2008) by Sara Silva (sara@dei.uc.pt)
%   Acknowledgements: Zheng Yin (zheng.yin@ucdconnect.ie)
%   This file is part of the GPLAB Toolbox

    [tmp ll]=size(X);
    op=tree.op;
    if isempty(tree.kids)
        for i=ll:-1:1
            op=strrep(op,['X',num2str(i)],['X(:,',num2str(i),')']);
        end       
        ret=op;
        return    
    end
    op=[tree.op,'('];
    ll=length(tree.kids);
    for i=1:ll
        argument=evaluate_tree(tree.kids{1,i},X);
        op=[op,argument];
        if i<ll
            op=[op,','];
        end
    end
    op=[op,')'];
    ret=mat2str(eval(op));
return