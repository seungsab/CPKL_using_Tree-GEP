
function Best = GP_SINGLE_KERNEL_STRUCTURE(x,y,k);
%% INPUT
% x : Input
% y : output
% k : kernel type

D = size(x,2);

%% INITIALIZE HYPER-PARAMETERS
% temp = [std(y)/sqrt(2)*[1e-4 4]]; % sigma
% temp = [temp; mean(std(x)) * [1e-1 3]]; % length scale
% temp = [temp; 0.01 3]; % Multiplier
% minn = temp(:,1)'; maxn = temp(:,2)';
% x0 = LHS(minn,maxn,30);

x0 = [std(y)/sqrt(2) mean(std(x)) 1];

%% Define Likelihood function
inf = @infGaussLik; lik = {@likGauss};
n_cnt = 1;
for i0 = 1:size(x0,1)
    hyp0=[]; covfunc1 = [];
    
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    switch k
        case 1 % covSEiso
            covfunc1 = {@covSEiso}; hyp0.cov = [ log(ell); log(sf)  ]; % Isotropic squared exponential (SE)
            
        case 2 % covSEard
            covfunc1 = {@covSEard}; hyp0.cov = [ log(ell)*ones(D,1); log(sf)  ]; % Isotropic squared exponential (SE)
            
        case 3 % covRQiso
            covfunc1 = {@covRQiso}; hyp0.cov = [ log(ell); log(sf); log(p)]; % Isotropic squared exponential (SE)
            
        case 4 % covRQard
            covfunc1 = {@covRQard}; hyp0.cov = [ log(ell)*ones(D,1); log(sf); log(p)]; % Isotropic squared exponential (SE)
            
        case 5 % covMaterniso, 3
            covfunc1 = {@covMaterniso, 3}; hyp0.cov = [ log(ell); log(sf)]; % Isotropic squared exponential (SE)
            
        case 6 % covMaternard, 3
            covfunc1 = {@covMaternard,3}; hyp0.cov = [ log(ell)*ones(D,1); log(sf)]; % Isotropic squared exponential (SE)
            
        case 7 % covMaterniso, 5
            covfunc1 = {@covMaterniso,5}; hyp0.cov = [ log(ell); log(sf)]; % Isotropic squared exponential (SE)
            
        case 8 % covMaternard, 5
            covfunc1 = {@covMaternard,5}; hyp0.cov = [ log(ell)*ones(D,1); log(sf)]; % Isotropic squared exponential (SE)
            
    end
    
    p = feval(covfunc1{:}); p = eval(p);
    covfunc{1,n_cnt} = covfunc1;
    
    try
    % FIT GPR
    hyp{1,n_cnt} = minimize(hyp0,@gp,-100,inf,[], covfunc1, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp{1,n_cnt}, inf, [], covfunc1, lik, x, y); % negative log marginal likelihood
    
    BIC(1,n_cnt) = 2*nlml + p*log(size(x,1));
    
    if isnan(BIC(1,n_cnt))
        BIC(1,n_cnt) = 1e8;
    end
    catch
        BIC(1,n_cnt) = 1e8;
    end
    n_cnt = n_cnt + 1;
end

[B IX] = min(BIC);
Best.BIC = B; Best.hyp = hyp{1,IX};
Best.covfunc = covfunc{1,IX};
end