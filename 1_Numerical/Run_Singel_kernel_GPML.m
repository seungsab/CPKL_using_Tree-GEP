%% Tutorial for optimzing the kernel function in Gaussian Process
function main
clc, clear all, close all, warning('off','all')

%% ADD GPML TOOLBOX PATH
addpath(genpath([pwd '\Training_data']))
addpath(genpath([pwd '\test_fcnt']))
addpath(genpath([pwd '\gpml']))
addpath(genpath([pwd '\utils']))

%% Define Noise level
% noise_amp = 0; %
% noise_amp = 0.003; %
% noise_amp = 0.005; %
% noise_amp = 0.01; %
% noise_amp = 0.03; %
% noise_amp = 0.05; %
% noise_amp = 0.1; %
for noise_amp = [0.003 0.005]
    for fcnt_type = [2 4 5 6 7]
        % fcnt_type = 6;
        switch fcnt_type
            case 1 % 1D
                lb = 0.5; ub= 2.5;
                fun = @(x)grlee12_fcnt(x);
                load Training_1D
                
            case 2 % 2D
                lb = [-5 0]; ub=[10 15];
                fun = @(x)branin_fcnt(x);
                load Training_2D
                
            case 3 % 2D
                lb = [0 0]; ub=[1 1];
                fun = @(x)curretal88exp(x);
                load Training_2D
                
            case 4 % 3D
                lb = [0 0 0]; ub=[1 1 1];
                fun = @(x)detpep10curv(x);
                load Training_3D
                
            case 5 % 5D
                lb = [0 0 0 0 0]; ub=[1 1 1 1 1];
                fun = @(x)fried_fcnt(x);
                load Training_5D
                
            case 6 % 8D
                lb = zeros(1,8); ub=ones(1,8);
                fun = @(x)detpep108d(x);
                load Training_8D
                
            case 7 % 20D
                lb = -0.5*ones(1,20); ub=0.5*ones(1,20);
                fun = @(x)welchetal92(x);
                load Training_20D
                
            case 8 % 2D
                lb = -10*ones(1,2); ub=10*ones(1,2);
                fun = @(x)levy2D(x);
                load Training_2D
                
            otherwise % 2D Adaptive sampling test fcnt
                % TO BE CORRECTED!!
                lb = 0; ub= 1;
                fun = @(x)grlee08_fcnt(x);
        end
        
        x1 = normalize_input(x_val,lb,ub,'x');
        for i = 1:size(x1,1)
            y_val(i,1) = fun(x1(i,:));
        end
        
        %% DEFINE KERNEL FUNCTION AND HYPER-PARAMETERS
        % sf = sigma;
        % ell (or L) = length scale
        % al = Multiplier (Power)
        % 1. @covSEiso; hypgi = log([ell;sf]);         % isotropic Gaussian
        % 2. @covSEard; hypga = log([L;sf]);           % Gaussian with ARD
        % 3. @covRQiso; hypri = log([ell;sf;al]);      % isotropic ration. quad.
        % 4. @covRQard; hypra = log([L;sf; al]);        % ration. quad. with ARD
        % 5. {@covMaterniso,3}; hypma = log([ell;sf]); % isotropic Matern class d=3
        % 6. {@covMaternard,3}; hypma = log([ell;sf]); % Matern class with d=3 and ARD
        % 7. {@covMaterniso,5}; hypma = log([ell;sf]); % isotropic Matern class d=5
        % 8. {@covMaternard,5}; hypma = log([ell;sf]); % Matern class with d=5 and ARD
        
        %% Convention: Train GP with optimizing hyper-parameters of fixed kernel type
        Optimal1 = {}; RMSE1 = []; RMAE1 = [];
        for i = 1:size(x0,2)
            for j=1:10
                % ASSIGN SUBSET OF TRAINING DATA
                x = x0{j,i};
                temp_x = normalize_input(x,lb,ub,'x'); y=[];
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
                
                n_ker = 1;
                for k=[2 4 6 8]
                    % FIT GPR ACCORDING TO DEFINED KERNEL
                    Best = GP_SINGLE_KERNEL_STRUCTURE(x,y,k);
                    BIC1 = Best.BIC; hyp1 = Best.hyp; covfunc1 = Best.covfunc;
                    
                    % PREDICT
                    try
                        [ypred s2] = gp(hyp1, @infGaussLik, [], covfunc1, @likGauss, x, y, x_val);
                    catch
                        [ypred s2] = gp(hyp1, @infGaussLik, [], covfunc1, @likGauss, x, y, x_val);
                    end
                    
                    % i => N  : # of Training Samples
                    % j => r1 : # of MC Runs
                    % k => N  : # of Kernel
                    RMSE1{i}(j,n_ker) = sqrt(sum((repmat(y_val,1,size(ypred,2))-ypred).^2)/size(y_val,1));
                    RMAE1{i}(j,n_ker)= max(abs(repmat(y_val,1,size(ypred,2))-ypred));
                    Optimal1{j,i} = Best;
                    n_ker = n_ker+1;
                end
                
                switch noise_amp
                    case 0
                        save(['Single_Kernel_GP' num2str(fcnt_type)]);
                    case 0.001
                        save(['Single_Kernel_GP' num2str(fcnt_type) '_Noise_0001']);
                    case 0.003
                        save(['Single_Kernel_GP' num2str(fcnt_type) '_Noise_0003']);
                    case 0.005
                        save(['Single_Kernel_GP' num2str(fcnt_type) '_Noise_0005']);
                    case 0.01
                        save(['Single_Kernel_GP' num2str(fcnt_type) '_Noise_001']);
                    case 0.03
                        save(['Single_Kernel_GP' num2str(fcnt_type) '_Noise_003']);
                    case 0.05
                        save(['Single_Kernel_GP' num2str(fcnt_type) '_Noise_005']);
                    case 0.1
                        save(['Single_Kernel_GP' num2str(fcnt_type) '_Noise_01']);
                end
            end
        end
    end
end
end