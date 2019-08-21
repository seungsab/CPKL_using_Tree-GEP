
function x1 = normalize_input(x,lb,ub,target_space)
%%
switch target_space
    case 'x' % [0 1]^p => original space
        x1 = x.*(ub - lb) + lb;
        
    case 'u' % original space => [0 1]^p
        x1 = (x - lb)./(ub - lb);
        
end
end