%% Tutorial for optimzing the kernel function in Gaussian Process
function main
clc, clear all, close all, warning('off','all')

%% ADD GPML TOOLBOX PATH
addpath(genpath([pwd '\Training_data']))
addpath(genpath([pwd '\test_fcnt']))
addpath(genpath([pwd '\gpml']))
addpath(genpath([pwd '\utils']))

%% Define Test data
% noise_amp = 0; %
% noise_amp = 0.05; %
% noise_amp = 0.1; %

for noise_amp = [0.003 0.005]
    for fcnt_type = [2 4 6]
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
                
            case 4 % 8D
                range = [0.05 100   63070  990  63.1 700 1120 1500;
                    0.15 50000 115600 1110 116  820 1680 15000];
                lb = range(1,:); ub = range(2,:);
                fun = @(x)borehole(x);
                load Training_8D
                
            case 5 % 8D
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
                
                % FIT GPR ACCORDING TO DEFINED KERNEL
                try
                    n_depth = 10;
                    [Best, Final, BIC_hist] = GP_Structure_Discovery_CKL(x,y,n_depth);
                catch
                    n_depth = 10;
                    [Best, Final, BIC_hist] = GP_Structure_Discovery_CKL(x,y,n_depth);
                    
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
                disp([ num2str(j) '/' num2str(size(x0,1))]);
                
                switch noise_amp
                    case 0
                        save(['AUTO_COMPOSITIONAL_GPR' num2str(fcnt_type)]);
                    case 0.003
                        save(['AUTO_COMPOSITIONAL_GPR' num2str(fcnt_type) '_Noise_0003']);
                    case 0.005
                        save(['AUTO_COMPOSITIONAL_GPR' num2str(fcnt_type) '_Noise_0005']);
                    case 0.01
                        save(['AUTO_COMPOSITIONAL_GPR' num2str(fcnt_type) '_Noise_001']);
                    case 0.03
                        save(['AUTO_COMPOSITIONAL_GPR' num2str(fcnt_type) '_Noise_003']);
                    case 0.05
                        save(['AUTO_COMPOSITIONAL_GPR' num2str(fcnt_type) '_Noise_005']);
                    case 0.1
                        save(['AUTO_COMPOSITIONAL_GPR' num2str(fcnt_type) '_Noise_01']);
                end
            end
        end
        
    end
end
end
