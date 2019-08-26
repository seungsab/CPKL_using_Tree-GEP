%% Compositional Kernel Learning using Tree-based Genetic Programming (CPKL_Tree-GEP) for Gaussian Proecess Regression
%% Method #2: Gaussian Process using Hierarchical CPK

% Updated for Replication of results

% This code is only working for the following test function
% Test function : Branin function (2-D problem)

% If you need source code for CPK using tree-GEP, feel free to contact me
% via e-mail (seungsab@gmail.com)
% Coded by SS (2019.08.21)

function main
clc, clear all, close all, warning('off','all')

%% ADD GPML TOOLBOX PATH
addpath(genpath([pwd '\Training_data']))
addpath(genpath([pwd '\test_fcnt']))
addpath(genpath([pwd '\gpml']))
addpath(genpath([pwd '\utils']))

%% Define Test data
noise_amp = 0;
% noise_amp = 0.005;
% noise_amp = 0.05;

N_rep = 3;
N_rep = min(N_rep,10); % Training data is generated 10 times

for fcnt_type = 1
    switch fcnt_type
        case 1 % 1D
            lb = [-5 0]; ub=[10 15];
            fun = @(x)branin_fcnt(x);
            load Training_2D
    end
    
    x1 = normalize_input(x_val,lb,ub,'x');
    for i = 1:size(x1,1)
        y_val(i,1) = fun(x1(i,:));
    end
    
    %% Convention: Train GP with optimizing hyper-parameters of fixed kernel type
    Optimal2 = {}; RMSE2 = []; RMAE2 = [];
    for i = 1:size(x0,2)
        for j=1:N_rep
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
            
            % FIT GPR ACCORDING TO DEFINED KERNEL
            try
                n_depth = 10;
                [Best, Final, BIC_hist] = GP_Structure_Discovery_CKL(x,y,n_depth,fun);
            catch
                n_depth = 10;
                [Best, Final, BIC_hist] = GP_Structure_Discovery_CKL(x,y,n_depth,fun);
                
            end
            BIC1 = Best.BIC; hyp1 = Best.hyp; covfunc1 = Best.covfunc;
            try
                % PREDICT
                [ypred s2] = gp(hyp1, @infGaussLik, [], covfunc1, @likGauss, x, y, x_val);
            catch
                % PREDICT
                [ypred s2] = gp(hyp1, @infGaussLik, [], covfunc1, @likGauss, x, y, x_val);
            end
            % i => N  : # of Training Samples
            % j => r1 : # of MC Runs
            % k => N  : # of Kernel
            RMSE2(j,i) = sqrt(sum((repmat(y_val,1,size(ypred,2))-ypred).^2)/size(y_val,1));
            RMAE2(j,i)= max(abs(repmat(y_val,1,size(ypred,2))-ypred));
            Optimal2{j,i} = Best;
            disp([ num2str(i) '/' num2str(size(x0,2))]);
            disp([ num2str(j) '/' num2str(N_rep)]);
            
            switch noise_amp
                case 0
                    save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type)]);
                case 0.005
                    save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type) '_Noise_0005']);
                case 0.05
                    save(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type) '_Noise_005']);
            end
        end
    end
end
end
