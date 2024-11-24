function [Optimal, Final, BIC_hist] = GP_Structure_Discovery_CKL(x,y,n_depth)
%% Initialize hyper-parameters
% temp = [std(y)/sqrt(2)*[1e-4 4]]; % sigma
% temp = [temp; mean(std(x)) * [1e-1 3]]; % length scale
% temp = [temp; 0.01 3]; % Multiplier
% minn = temp(:,1)'; maxn = temp(:,2)';
% x0 = LHS(minn,maxn,30);
x0 = [std(y)/sqrt(2) mean(std(x)) 1];

%% SEARCH STRUCTURAL DISCOVERY
BIC_hist = [];
ind_depth = 1; stop_ind = 1;
while (ind_depth<=n_depth && stop_ind)
    if ind_depth == 1
        % DEFINE initial values of KERNEL FUNCTION
        Best = GP_struct_1_depth(x0,x,y);
        BIC1 = Best.BIC; hyp1 = Best.hyp; covfunc1 = Best.covfunc;
        
    else
        BIC0 = BIC1; hyp0 = hyp1; covfunc0 = covfunc1;
        Best = GP_struct_n_depth(x0,covfunc0,hyp0,x,y);
        BIC1 = Best.BIC; hyp1 = Best.hyp; covfunc1 = Best.covfunc;
        
    end
    
    if isempty(BIC_hist)
        BIC_hist = [BIC_hist BIC1];
        Optimal.covfunc = covfunc1; Optimal.hyp = hyp1;
        Optimal.BIC = BIC1;
        Final.covfunc = covfunc1; Final.hyp = hyp1;
        Final.BIC = BIC1;
        
    else
        BIC_hist = [BIC_hist BIC1];
        if BIC_hist(end) >= BIC_hist(end-1)
            stop_ind = 0;
            Optimal.covfunc = covfunc0; Optimal.hyp = hyp0;
            Final.covfunc = covfunc1; Final.hyp = hyp1;
            Optimal.BIC = BIC_hist(end-1);  Final.BIC = BIC_hist(end);
        end
    end
    ind_depth = ind_depth + 1;
end

end

function Best = GP_struct_1_depth(x0,x,y)
%% Define Likelihood function
inf = @infGaussLik; lik = {@likGauss};

nlml = []; hyp = {}; Best = []; n_cnt = 1;
D = size(x,2);
%% BASE KERNEL #1: White Noise
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    covfunc{1,n_cnt} = 'covNoise'; hyp0.cov = [log(sf)]; % Noise
    hyp{1,n_cnt} = minimize(hyp0,@gp,-100,inf,[], covfunc{1,n_cnt}, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp{1,n_cnt}, inf, [], covfunc{1,n_cnt}, lik, x, y); % negative log marginal likelihood
    try
        p = eval(feval(covfunc{1,n_cnt}{:}));
    catch
        p = eval(feval(covfunc{1,n_cnt}));
    end
    BIC(1,n_cnt) = 2*nlml + p*log(size(x,1));
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #2: Constant
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    covfunc{1,n_cnt} = 'covConst'; hyp0.cov = [log(sf)];
    hyp0.cov = [max([log(ell) log(sf)])];
    hyp{1,n_cnt} = minimize(hyp0,@gp,-100,inf,[], covfunc{1,n_cnt}, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp{1,n_cnt}, inf, [], covfunc{1,n_cnt}, lik, x, y); % negative log marginal likelihood
    try
        p = eval(feval(covfunc{1,n_cnt}{:}));
    catch
        p = eval(feval(covfunc{1,n_cnt}));
    end
    BIC(1,n_cnt) = 2*nlml + p*log(size(x,1));
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #3: LIN
covLIN = {'covSum', {'covLINard','covConst'}};
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    covfunc{1,n_cnt} = covLIN;
    hyp0.cov = [min([log(ell) log(sf)])*ones(D,1);  max([log(ell) log(sf)])];
    hyp{1,n_cnt} = minimize(hyp0,@gp,-100,inf,[], covfunc{1,n_cnt}, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp{1,n_cnt}, inf, [], covfunc{1,n_cnt}, lik, x, y); % negative log marginal likelihood
    try
        p = eval(feval(covfunc{1,n_cnt}{:}));
    catch
        p = eval(feval(covfunc{1,n_cnt}));
    end
    BIC(1,n_cnt) = 2*nlml + p*log(size(x,1));
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #4: Periodic
% periodic Matern with ARD
% covPER = {'covProd', {'covPERard','covSEard'}};
% covPER = {'covProd', {'covPeriodic', 'covSEard'}};
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    
    % Now construct the kernel.
    kernel_components = cell(1, D);
    %     kernel_hypers = cell(1, D);
    kernel_hypers = [];
    for i = 1:D
        kernel_components{i} = { 'covMask', { i, 'covPeriodic'}};
        kernel_hypers = [kernel_hypers; log(ell); log(p); log(sf)];
        %         kernel_hypers{i} = [log(ell); log(p); log(sf)];
    end
    % concatenate all additive components
    covPER = { 'covSum', kernel_components };
    
    covfunc{1,n_cnt} = covPER;
    hyp0.cov = kernel_hypers;
    %     try
    %         hyp0.cov = unwrap(kernel_hypers);
    %     catch
    %         hyp0.cov = unwrap(kernel_hypers{:});
    %     end
    hyp{1,n_cnt} = minimize(hyp0,@gp,-100,inf,[], covfunc{1,n_cnt}, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp{1,n_cnt}, inf, [], covfunc{1,n_cnt}, lik, x, y); % negative log marginal likelihood
    try
        p = eval(feval(covfunc{1,n_cnt}{:}));
    catch
        p = eval(feval(covfunc{1,n_cnt}));
    end
    BIC(1,n_cnt) = 2*nlml + p*log(size(x,1));
    n_cnt = n_cnt + 1;
end


%% BASE KERNEL #5: covSEard
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    covfunc{1,n_cnt} = 'covSEard';
    hyp0.cov = [ log(ell)*ones(D,1); log(sf)  ]; % Isotropic squared exponential (SE)
    hyp{1,n_cnt} = minimize(hyp0,@gp,-100,inf,[], covfunc{1,n_cnt}, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp{1,n_cnt}, inf, [], covfunc{1,n_cnt}, lik, x, y); % negative log marginal likelihood
    try
        p = eval(feval(covfunc{1,n_cnt}{:}));
    catch
        p = eval(feval(covfunc{1,n_cnt}));
    end
    BIC(1,n_cnt) = 2*nlml + p*log(size(x,1));
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #6: covRQiard
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    covfunc{1,n_cnt} = 'covRQard';
    hyp0.cov = [ log(ell)*ones(D,1); log(sf); log(p) ]; % Isotropic rational quadratic (RQ)
    hyp{1,n_cnt} = minimize(hyp0,@gp,-100,inf,[], covfunc{1,n_cnt}, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp{1,n_cnt}, inf, [], covfunc{1,n_cnt}, lik, x, y); % negative log marginal likelihood
    try
        p = eval(feval(covfunc{1,n_cnt}{:}));
    catch
        p = eval(feval(covfunc{1,n_cnt}));
    end
    BIC(1,n_cnt) = 2*nlml + p*log(size(x,1));
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #7: covMaternard3
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    covfunc{1,n_cnt} = {@covMaternard,3};
    hyp0.cov = [ log(ell)*ones(D,1); log(sf) ]; % Isotropic rational quadratic (RQ)
    hyp{1,n_cnt} = minimize(hyp0,@gp,-100,inf,[], covfunc{1,n_cnt}, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp{1,n_cnt}, inf, [], covfunc{1,n_cnt}, lik, x, y); % negative log marginal likelihood
    try
        p = eval(feval(covfunc{1,n_cnt}{:}));
    catch
        p = eval(feval(covfunc{1,n_cnt}));
    end
    BIC(1,n_cnt) = 2*nlml + p*log(size(x,1));
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #8: covMaternard5
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    covfunc{1,n_cnt} = {@covMaternard,5};
    hyp0.cov = [ log(ell)*ones(D,1); log(sf) ]; % Isotropic rational quadratic (RQ)
    hyp{1,n_cnt} = minimize(hyp0,@gp,-100,inf,[], covfunc{1,n_cnt}, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp{1,n_cnt}, inf, [], covfunc{1,n_cnt}, lik, x, y); % negative log marginal likelihood
    try
        p = eval(feval(covfunc{1,n_cnt}{:}));
    catch
        p = eval(feval(covfunc{1,n_cnt}));
    end
    BIC(1,n_cnt) = 2*nlml + p*log(size(x,1));
    n_cnt = n_cnt + 1;
end



[B IX] = min(BIC);
Best.BIC = B; Best.hyp = hyp{1,IX};
Best.covfunc = covfunc{1,IX};

end

function Best = GP_struct_n_depth(x0,covfunc1,hyp1,x,y)
%% Define Likelihood function
inf = @infGaussLik; lik = {@likGauss};
nlml = []; hyp = {}; Best = []; n_cnt = 1;
D = size(x,2);
%% Summation part
%% BASE KERNEL #1: + White Noise
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(sf)]; % Noise
    covfunc{1,n_cnt} = {'covSum', {covfunc1,'covNoise'}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #2: + CONST
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; max([log(ell) log(sf)])];
    covfunc{1,n_cnt} = {'covSum', {covfunc1,'covConst'}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #3: + LIN
covLIN = {'covSum', {'covLINard','covConst'}};
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; min([log(ell) log(sf)])*ones(D,1);  max([log(ell) log(sf)])];
    covfunc{1,n_cnt} = {'covSum', {covfunc1,covLIN}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end
%% BASE KERNEL #4: + Periodic
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    
    % Now construct the kernel.
    kernel_components = cell(1, D);
    kernel_hypers = [];
    for i = 1:D
        kernel_components{i} = { 'covMask', { i, 'covPeriodic'}};
        kernel_hypers = [kernel_hypers; log(ell); log(p); log(sf)];
    end
    % concatenate all additive components
    covPER = { 'covSum', kernel_components };
    covfunc{1,n_cnt} = {'covSum', {covfunc1,covPER}};
    
    %     try
    %         temp_cov = unwrap(kernel_hypers);
    %     catch
    %         temp_cov = unwrap(kernel_hypers{:});
    %     end
    hyp0.cov = [hyp1.cov; kernel_hypers];
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end
%% BASE KERNEL #5: + SE
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(ell)*ones(D,1); log(sf)]; % Linear with bias
    covfunc{1,n_cnt} = {'covSum', {covfunc1,'covSEard'}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end
%% BASE KERNEL #6: + RQ
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(ell)*ones(D,1); log(sf); log(p) ]; % Isotropic rational quadratic (RQ)
    covfunc{1,n_cnt} = {'covSum', {covfunc1,'covRQard'}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end
%% BASE KERNEL #7: + covMaternard3
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(ell)*ones(D,1); log(sf) ]; % Isotropic rational quadratic (RQ)
    temp_cov1 = {@covMaternard,3};
    covfunc{1,n_cnt} = {'covSum', {covfunc1,temp_cov1}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end
%% BASE KERNEL #8: + covMaternard5
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(ell)*ones(D,1); log(sf) ]; % Isotropic rational quadratic (RQ)
    temp_cov1 = {@covMaternard,5};
    covfunc{1,n_cnt} = {'covSum', {covfunc1,temp_cov1}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end


%% Product part
%% BASE KERNEL #1: x White Noise
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(sf)]; % Noise
    covfunc{1,n_cnt} = {'covProd', {covfunc1,'covNoise'}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #2: X CONST
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; max([log(ell) log(sf)])];
    covfunc{1,n_cnt} = {'covProd', {covfunc1,'covConst'}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #3: x LINEAR
covLIN = {'covSum', {'covLINard','covConst'}};
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(ell)*ones(D,1);  log(sf) ];
    covfunc{1,n_cnt} = {'covProd', {covfunc1,covLIN}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end
%% BASE KERNEL #4: x Periodic
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    
    % Now construct the kernel.
    kernel_components = cell(1, D);
    kernel_hypers = [];
    for i = 1:D
        kernel_components{i} = { 'covMask', { i, 'covPeriodic'}};
        kernel_hypers = [kernel_hypers; log(ell); log(p); log(sf)];
    end
    % concatenate all additive components
    covPER = { 'covSum', kernel_components };
    covfunc{1,n_cnt} = {'covProd', {covfunc1,covPER}};
    %     try
    %         temp_cov = unwrap(kernel_hypers);
    %     catch
    %         temp_cov = unwrap(kernel_hypers{:});
    %     end
    hyp0.cov = [hyp1.cov; kernel_hypers];
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end

%% BASE KERNEL #5: x SE
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(ell)*ones(D,1); log(sf)];
    covfunc{1,n_cnt} = {'covProd', {covfunc1,'covSEard'}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end
%% BASE KERNEL #6: X RQ
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(ell)*ones(D,1); log(sf); log(p) ]; % Isotropic rational quadratic (RQ)
    covfunc{1,n_cnt} = {'covProd', {covfunc1,'covRQard'}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end
%% BASE KERNEL #7: X covMaternard3
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(ell)*ones(D,1); log(sf) ]; % Isotropic rational quadratic (RQ)
    temp_cov1 = {@covMaternard,3};
    covfunc{1,n_cnt} = {'covProd', {covfunc1,temp_cov1}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end
%% BASE KERNEL #8: X covMaternard5
for i0 = 1:size(x0,1)
    % DEFINE COV FUNCTION
    hyp0.lik = x0(i0,1); sf = x0(i0,1); ell = x0(i0,2); p = x0(i0,3);
    hyp0.cov = [hyp1.cov; log(ell)*ones(D,1); log(sf) ]; % Isotropic rational quadratic (RQ)
    temp_cov1 = {@covMaternard,5};
    covfunc{1,n_cnt} = {'covProd', {covfunc1,temp_cov1}};
    [BIC(1,n_cnt), hyp{1,n_cnt}] = Composition_kernel(covfunc{1,n_cnt},hyp0,x,y);
    n_cnt = n_cnt + 1;
end

[B IX] = min(BIC);
Best.BIC = B; Best.hyp = hyp{IX};
Best.covfunc = covfunc{IX};
end

function [BIC, hyp] = Composition_kernel(covfunc,hyp0,x,y)
D = size(x,2);
%% Define Likelihood function
inf = @infGaussLik; lik = {@likGauss};
hyp = minimize(hyp0,@gp,-100,inf,[], covfunc, lik,x,y); % optimise hyperparameters
nlml = gp(hyp, inf, [], covfunc, lik, x, y); % negative log marginal likelihood
try
    p = eval(feval(covfunc{:}));
catch
    p = eval(feval(covfunc));
end
BIC = 2*nlml + p*log(size(x,1));
end