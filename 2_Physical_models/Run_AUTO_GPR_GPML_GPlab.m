%% Tutorial for optimzing the kernel function in Gaussian Process
function main
clc, clear all, close all, warning('off','all')

%% ADD GPML TOOLBOX PATH
global x y lb ub fun
addpath(genpath([pwd '\Training_data']))
addpath(genpath([pwd '\test_fcnt']))
addpath(genpath([pwd '\gplab']))
addpath(genpath([pwd '\gpml']))
addpath(genpath([pwd '\utils']))

%% Define Test data
% noise_amp = 0; %
% noise_amp = 0.05; %
% noise_amp = 0.1; %

for noise_amp = [0.003 0.005]
    for fcnt_type = [1 2 4 6]
        switch fcnt_type
            case 1 % 3D
                range = [3.5 800 200
                    6.5 3200 800];
                lb = range(1,:); ub = range(2,:);
                fun = @(x)shortcol(x);
                load Training_3D
                
            case 2 % 6D
                range = [50 25 0.5 1.2 0.25 50
                    150 70 3 2.5 1.2 300];
                lb = range(1,:); ub = range(2,:);
                fun = @(x)otlcircuit(x);
                load Training_6D
                
            case 3 % 7D
                range = [30 0.005 0.002 1000 90000  290 340
                    60 0.020 0.010 5000 110000 296 360];
                lb = range(1,:); ub = range(2,:);
                fun = @(x)piston_7D(x);
                load Training_7D
                
            case 4 % 8D borehole
                range = [0.05 100   63070  990  63.1 700 1120 1500;
                    0.15 50000 115600 1110 116  820 1680 15000];
                lb = range(1,:); ub = range(2,:);
                fun = @(x)borehole(x);
                load Training_8D
                
            case 5 % 8D robot
                range = [0 0 0 0 0 0 0 0
                    2*pi*ones(1,4) ones(1,4)];
                lb = range(1,:); ub = range(2,:);
                fun = @(x)Robotarm(x);
                load Training_8D
                
            case 6 % 10D
                range = [150 220 6 -10 16 0.5 0.08 2.5 1700 0.025
                    200 300 10 10 45 1 0.18 6 2500 0.08];
                lb = range(1,:); ub = range(2,:);
                fun = @(x)wingweight(x);
                load Training_10D
        end
        
        x1 = normalize_input(x_val,lb,ub,'x'); y_val = [];
        for i = 1:size(x1,1)
            y_val(i,1) = fun(x1(i,:));
        end
        
        %% Convention: Train GP with optimizing hyper-parameters of fixed kernel type
        Optimal2 = {}; RMSE2 = []; RMAE2 = [];
        for i = 1:size(x0,2)
            for j=1:10
                x = []; y = [];
                % ASSIGN SUBSET OF TRAINING DATA
                x = x0{j,i};
                temp_x = normalize_input(x,lb,ub,'x');
                for ijk = 1:size(temp_x,1)
                    y(ijk,:) = fun(temp_x(ijk,:));
                end
                
                %% NOISE INJECTION
                seedNum = j*123456789;                  % 2. Seed for random number generator; put [] for autmatic randomization
                %% Randomization
                if isempty( seedNum ) == false
                    rand('state',seedNum);
                end
                
                RMS=sqrt(sum((y.^2))/size(y,1));
                noise = RMS*noise_amp*randn(size(y,1),1);
                RMS_n=sqrt(sum((noise.^2))/size(noise,1));
                
                % CONSTRUCT NOISY DATA
                y = y + noise;
                
                %% GE PROGRAMMING
                % Set parameters of the GE programming
                Set_GPlab_SS
                
                % RUN GE Programming
                pop = 200; gen = 10;
                [v,b]=gplab(gen,pop,p);
                
                Best = composite_best_kernel(b.tree,x,y);
                
                % PREDICT
                [ypred s2] = gp(Best.hyp_final, @infGaussLik, [], Best.covfunc_final, @likGauss, x, y, x_val);
                
                % i => N  : # of Training Samples
                % j => r1 : # of MC Runs
                % k => N  : # of Kernel
                RMSE2(j,i) = sqrt(sum((repmat(y_val,1,size(ypred,2))-ypred).^2)/size(y_val,1));
                RMAE2(j,i)= max(abs(repmat(y_val,1,size(ypred,2))-ypred));
                Optimal2{j,i} = Best;
                disp([ num2str(i) '/' num2str(size(x0,2))]);
                disp([ num2str(j) '/' num2str(size(x0,1))]);
                
                switch noise_amp
                    case 0
                        save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type)]);
                    case 0.003
                        save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type) '_Noise_0003']);
                    case 0.005
                        save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type) '_Noise_0005']);
                    case 0.01
                        save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type) '_Noise_001']);
                    case 0.03
                        save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type) '_Noise_003']);
                    case 0.05
                        save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type) '_Noise_005']);
                    case 0.1
                        save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type) '_Noise_01']);
                end
            end
        end
        
    end
end
end

function Best = composite_best_kernel(tree,x,y)
%% Define Likelihood function
inf = @infGaussLik; lik = {@likGauss};

%% Initialize hyper-parameters
hyp0 = [std(y)/sqrt(2) mean(std(x)) 1];

%% TREE STRUCTURE
tree_trace = get_treeinfor(tree,[],[]);
tree_trace(4,:) =  tree_trace(3,:)<1;
tree_trace(5,:) =  0;

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
                temp_ind = tree_trace(2,i_depth); temp_ind = tree_trace(4,temp_ind);
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
                    temp_ind = tree_trace(2,i_depth); temp_ind = tree_trace0(4,temp_ind);
                    COV_compo{1,temp_ind} = COV1; hyp_compo{1,temp_ind} = hyp1;
                    tree_trace(5,i_depth:i_depth+1)=1;
                end
                
            end
        end
    end
    
    covfunc_final = COV_compo{1};
    hyp_final = hyp_compo{1}; hyp_final.lik = hyp0(1);
    
end

%% INFERENCE of GPR
hyp_final = minimize(hyp_final,@gp,-100,inf,[], covfunc_final, lik,x,y); % optimise hyperparameters
nlml = gp(hyp_final, inf, [], covfunc_final, lik, x, y); % negative log marginal likelihood

np = size(hyp_final.cov,1);

BIC = 2*nlml + np*log(size(x,1)); % BIC

Best.BIC = BIC; Best.hyp_final = hyp_final; Best.covfunc_final = covfunc_final;
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