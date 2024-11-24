
%% INITIALIZE
clear all, close all, clc
addpath(genpath([pwd '\gpml']))
addpath(genpath([pwd '\utils']))

%% GENERATE R.V.
if isdir('MC_1E6.mat')
    load MC_1E6.mat X_MC
else    
    mu = [0.07433 0.1 13 4751 -648];
    sigma = [0.005 0.01 60 48 11];
    
    N_MC = 1e6;
    X_MC = [];
    for i=1:size(mu,2)
        m = mu(1,i);
        v = sigma(1,i);
        X_MC(:,i)=normrnd(m,v,N_MC,1);
    end
    
    save MC_1E6.mat X_MC
end

%% MCS
lb = [0.04,0.03,-375,4500,-720];
ub = [0.15,0.20,400,5100,-580];
X_MC_normal = normalize_input(X_MC,lb,ub,'u');

y = uq_liquidHydrogenTank(X_MC);

%% PLOT PDF
hist(y)

%% PLOT structural failure probability
[f,x] = ecdf(y);
[B IX] = find(x<0);

p_f0 = size(B,1)/size(X_MC,1);

%% Single KERNEL
load Single_Kernel_GP11

p_f1 = [];

for i=1:size(Optimal1,2)
    for j = 1:size(Optimal1,1)
        x = x0{j,i}; y=[];
        temp_x = normalize_input(x,lb,ub,'x');
        for ijk = 1:size(temp_x,1)
            y(ijk,:) = fun(temp_x(ijk,:));
        end
                
        hyp1 = Optimal1{j,i}.hyp; covfunc1 = Optimal1{j,i}.covfunc;
        [ypred s2] = gp(hyp1, @infGaussLik, [], covfunc1, @likGauss, x, y, X_MC_normal);
        
        [f,x] = ecdf(ypred);
        [B IX] = find(x<0);
        
        p_f1(j,i) = size(B,1)/size(X_MC_normal,1);
    end
end

%% CPK
load AUTO_COMPOSITIONAL_GPR11
RMSE = RMSE2; RMAE = RMAE2;

p_f2 = [];

for i=1:size(Optimal2,2)
    for j = 1:size(Optimal2,1)
        x = x0{j,i}; y=[];
        temp_x = normalize_input(x,lb,ub,'x');
        for ijk = 1:size(temp_x,1)
            y(ijk,:) = fun(temp_x(ijk,:));
        end
                
        hyp1 = Optimal2{j,i}.hyp; covfunc1 = Optimal2{j,i}.covfunc;
        [ypred s2] = gp(hyp1, @infGaussLik, [], covfunc1, @likGauss, x, y, X_MC_normal);
        
        [f,x] = ecdf(ypred);
        [B IX] = find(x<0);
        
        p_f2(j,i) = size(B,1)/size(X_MC_normal,1);
    end
end

%% CPK
clear Optimal2

load AUTO_COMPOSITIONAL_GEP_GPR11
RMSE3 = RMSE2; RMAE3 = RMAE2;
RMSE2 = RMSE; RMAE2 = RMAE;

p_f3 = [];

for i=1:size(Optimal2,2)
    for j = 1:size(Optimal2,1)
        x = x0{j,i}; y=[];
        temp_x = normalize_input(x,lb,ub,'x');
        for ijk = 1:size(temp_x,1)
            y(ijk,:) = fun(temp_x(ijk,:));
        end
                
        hyp1 = Optimal2{j,i}.hyp; covfunc1 = Optimal2{j,i}.covfunc;
        [ypred3 s2] = gp(hyp1, @infGaussLik, [], covfunc1, @likGauss, x, y, X_MC_normal);
        
        [f,x] = ecdf(ypred);
        [B IX] = find(x<0);
        
        p_f3(j,i) = size(B,1)/size(X_MC_normal,1);
    end
end

save MCR_result_liquid_tank x0 p_f0 p_f1 p_f2 p_f3

% figure;
% plot(x,p_f,'b-','linewidth',2); hold on
% plot(x(B(end)),p_f(B(end)),'ro','linewidth',2);
% plot(x(IX)*[1 1],[0 p_f(IX)],'r:','linewidth',2);
% xlabel('Performance function');
% ylabel('Probability of Exceedance');
% set(gca,'yscale','log','fontsize',15,'fontweight','bold');
% hold on; grid on
% xlim([0 1.5]);
%%
figure
plot(p_f1,'ro'); hold on
plot(p_f2,'go');
plot(p_f3,'bo');

