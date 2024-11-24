
function ind = GP_Structure_Discovery_GEP(ind,params,data,~,~)
global x y lb ub fun

% x1: WN
% x2: CON
% x3: LIN
% x4: PER
% x5: SEard
% x6: RQard
% x7: Maternard3
% x8: Maternard5

%% Define Likelihood function
inf = @infGaussLik; lik = {@likGauss};

%% Initialize hyper-parameters
hyp0 = [std(y)/sqrt(2) mean(std(x)) 1];

%% TREE STRUCTURE
tree = ind.tree;
% drawtree(tree);

tree_trace = get_treeinfor(tree,[],[]);
tree_trace(4,:) =  tree_trace(3,:)<1;
tree_trace(5,:) =  0;
try
    if size(tree_trace,2) == 1
        kernel_type = tree_trace(3,1);
        [covfunc_final, hyp_final] = Single_kernel(hyp0,x,y,kernel_type);
        
    else
        % CONSTUCT TABLE FOR hirachical tree
        label_compo = find(tree_trace(4,:)==1);
        COV_compo = cell(1,size(label_compo,2));
        hyp_compo = cell(1,size(label_compo,2));
        nCount = 1;
        
        for i000=1:size(tree_trace,2)
            if tree_trace(4,i000) == 1
                tree_trace(4,i000) = nCount; nCount = nCount + 1;
            end
        end
        tree_trace0 = tree_trace;
        
        %% CONSTRUCT COMPOSITIONAL KERNEL
        for i_depth = 1:size(tree_trace,2)-1
            if tree_trace(4,i_depth) == tree_trace(4,i_depth+1)
                kernel_type = [tree_trace(3,i_depth-1) tree_trace(3,i_depth) tree_trace(3,i_depth+1)];
                
                if kernel_type(1)<0
                    % COMPONENTs for COMPOSITIONAL KERNEL
                    if kernel_type(2)>0
                        [covfuncL, hypL] = Single_kernel(hyp0,x,y,kernel_type(2));
                    else % Previously compositional kernel
                        [covfuncL, hypL] = Single_kernel(hyp0,x,y,kernel_type(2));a
                    end
                    
                    if kernel_type(3)>0
                        [covfuncR, hypR] = Single_kernel(hyp0,x,y,kernel_type(3));
                    else % Previously compositional kernel
                        [covfuncR, hypR] = Single_kernel(hyp0,x,y,kernel_type(3));
                    end
                    
                    [COV1, hyp1] = Composition_kernel(covfuncL,hypL,covfuncR,hypR,x,y,kernel_type(1));
                    temp_ind = tree_trace(2,i_depth);
                    temp_ind=find(temp_ind == tree_trace(1,:));
                    temp_ind = tree_trace(4,temp_ind);
                    COV_compo{1,temp_ind} = COV1; hyp_compo{1,temp_ind} = hyp1;
                    tree_trace(5,i_depth:i_depth+1)=1;
                end
            end
        end
        
        while isempty(COV_compo{1,1})
            temp_ind = find(tree_trace(5,:)==0);
            tree_trace=tree_trace(:,temp_ind);
            
            for i_depth = 1:size(tree_trace,2)-1
                if tree_trace(2,i_depth) == tree_trace(2,i_depth+1)
                    kernel_type = [tree_trace(3,i_depth-1) tree_trace(3,i_depth) tree_trace(3,i_depth+1)];
                    run_comp = 1;
                    
                    % COMPONENTs for COMPOSITIONAL KERNEL
                    if kernel_type(2)>0
                        [covfuncL, hypL] = Single_kernel(hyp0,x,y,kernel_type(2));
                    else % Previously compositional kernel
                        ind1 = tree_trace(4,i_depth);
                        if isempty(COV_compo{1,ind1})
                            run_comp = 0;
                        else
                            covfuncL = COV_compo{1,ind1};
                            hypL = hyp_compo{1,ind1};
                        end
                    end
                    
                    if kernel_type(3)>0
                        [covfuncR, hypR] = Single_kernel(hyp0,x,y,kernel_type(3));
                    else % Previously compositional kernel
                        ind1 = tree_trace(4,i_depth+1);
                        if isempty(COV_compo{1,ind1})
                            run_comp = 0;
                        else
                            covfuncR = COV_compo{1,ind1};
                            hypR = hyp_compo{1,ind1};
                        end
                    end
                    if run_comp
                        [COV1, hyp1] = Composition_kernel(covfuncL,hypL,covfuncR,hypR,x,y,kernel_type(1));
                        temp_ind = tree_trace(2,i_depth); 
                        temp_ind=find(temp_ind == tree_trace(1,:));
                        temp_ind = tree_trace(4,temp_ind);
                        COV_compo{1,temp_ind} = COV1; hyp_compo{1,temp_ind} = hyp1;
                        tree_trace(5,i_depth:i_depth+1)=1;
                    end
                    
                end
            end
        end
        
        covfunc_final = COV_compo{1};
        hyp_final = hyp_compo{1}; hyp_final.lik = hyp0(1);
        
    end
catch
    drawtree(tree);
    tree
end

np = size(hyp_final.cov,1);

try
    %% INFERENCE of GPR
    hyp_final = minimize(hyp_final,@gp,-100,inf,[], covfunc_final, lik,x,y); % optimise hyperparameters
    nlml = gp(hyp_final, inf, [], covfunc_final, lik, x, y); % negative log marginal likelihood
    ind.result = nlml;
    
    %% FITNESS FUNCTION (BIC)
    
    ind.fitness = 2*nlml + np*log(size(x,1)); % BIC
    
    % now limit fitness precision, to eliminate rounding error problem:
    ind.fitness=fixdec(ind.fitness,params.precision);
    
    if isnan(ind.fitness)
        ind.fitness = 1e8;
    end
catch
    ind.fitness = 1e8;
end

end

function [covfunc, hyp0] = Single_kernel(x0,x,y,kernel_type)
D = size(x,2);
hyp0.lik = x0(1,1); sf = x0(1,1); ell = x0(1,2); p = x0(1,3);

switch kernel_type
    case 1
        %% BASE KERNEL #1: White Noise
        covfunc = 'covNoise';
        hyp0.cov = [log(sf)]; % Noise
        
    case 2
        %% BASE KERNEL #2: Constant
        covfunc = 'covConst';
        hyp0.cov = [max([log(ell) log(sf)])];
        
    case 3
        %% BASE KERNEL #3: LIN
        covLIN = {'covSum', {'covLINard','covConst'}};
        covfunc = covLIN;
        hyp0.cov = [min([log(ell) log(sf)])*ones(D,1);  max([log(ell) log(sf)])];
        
    case 4
        %% BASE KERNEL #4: Periodic
        % periodic Matern with ARD
        % covPER = {'covProd', {'covPERard','covSEard'}};
        % covPER = {'covProd', {'covPeriodic', 'covSEard'}};
        kernel_components = cell(1, D);
        kernel_hypers = [];
        for i = 1:D
            kernel_components{i} = { 'covMask', { i, 'covPeriodic'}};
            kernel_hypers = [kernel_hypers; log(ell); log(p); log(sf)];
        end
        % concatenate all additive components
        covPER = { 'covSum', kernel_components };
        covfunc = covPER;
        hyp0.cov = kernel_hypers;
        
    case 5
        %% BASE KERNEL #5: covSEard
        covfunc = 'covSEard';
        hyp0.cov = [ log(ell)*ones(D,1); log(sf)  ]; % Isotropic squared exponential (SE)
        
    case 6
        %% BASE KERNEL #6: covRQard
        covfunc = 'covRQard';
        hyp0.cov = [ log(ell)*ones(D,1); log(sf); log(p) ]; % Isotropic rational quadratic (RQ)
        
    case 7
        %% BASE KERNEL #7: covMaternard3
        covfunc = {@covMaternard,3};
        hyp0.cov = [ log(ell)*ones(D,1); log(sf) ]; % Isotropic rational quadratic (Maternard3)
        
    case 8
        %% BASE KERNEL #8: covMaternard5
        covfunc = {@covMaternard,5};
        hyp0.cov = [ log(ell)*ones(D,1); log(sf) ]; % Isotropic rational quadratic (Maternard5)

end

end

function [covfunc, hyp] = Composition_kernel(covfuncL,hypL,covfuncR,hypR,x,y,kernel_type);
%% Define Likelihood function
inf = @infGaussLik; lik = {@likGauss};
hyp0 = [std(y)/sqrt(2) mean(std(x)) 1];

D = size(x,2);
switch kernel_type
    case -1 % plus
        covfunc = {'covSum', {covfuncL,covfuncR}};
    case -2 % times
        covfunc = {'covProd', {covfuncL,covfuncR}};
end

hyp.cov = [hypL.cov; hypR.cov]; hyp.lik = hyp0(1);
hyp = minimize(hyp,@gp,-100,inf,[], covfunc, lik,x,y); % optimise hyperparameters

end

function s =get_treeinfor(tree, s, parentID)

nkids=size(tree.kids,2);


temp(1,1) = tree.nodeid;
if isempty(parentID)
    parentID = 0;
end
temp(2,1) = parentID;

switch tree.op
    case 'plus'
        temp(3,1) = -1;
    case 'times'
        temp(3,1) = -2;
    otherwise
        temp(3,1) = str2num(tree.op(2:end));
        
end

s = [s temp];

if nkids>0
    for i=1:nkids
        [s]=get_treeinfor(tree.kids{i},s, tree.nodeid);
    end
end
end