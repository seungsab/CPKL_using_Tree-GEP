

% t   = xx(1); % t: Thickness of tube (mm) //  Normal(5.0,0.5) => [2 8]
% d   = xx(2); % d: Outside diameter (mm) //  Normal(42.0,4.2) => [20 65]
% L1  = xx(3); % L1:  Distance of F1 (mm) //  Normal(120.0,12.0) => [60 180]
% L2  = xx(4); % L2:  Distance of F2 (mm) //  Normal(50.0,5.0) => [20 80]
% F1  = xx(5); % F1: External Force #1 with the angle of 15 deg (kN) // Log-normal (3.0,0.1)
% F2  = xx(6); % F2: External Force #2 with the angle of 25 deg (kN) // Log-normal (3.0,0.1)
% P   = xx(7); % P: External Force #3 (axial force) (kN) // Log-normal (12.0,0.1)
% T   = xx(8); % T: Tortional Force (kN)  // Log-normal (100.0,0.1)
% Sy  = xx(9); % Sy: Material Yield Strength (MPa) // Log-normal (220.0,0.1)
clear all, close all, clc

X = [];
mu = [0.07433 0.1 13 4751 -648];
sigma = [0.005 0.01 60 48 11];
N_MC = 1e6;

n_cnt = 1;
for i0=1:size(mu,2)
    for i=1:1e4
        if i0<6
            m = mu(1,n_cnt);
            v = sigma(1,n_cnt);
            X1=normrnd(m,v,N_MC,1);
        else
            m = mu(1,n_cnt);
            v = m*sigma(1,n_cnt);
            m1 = log((m^2)/sqrt(v^2+m^2));
            v1 = sqrt(log(v^2/(m^2)+1));
            X1 = lognrnd(m1,v1,N_MC,1);
        end
        
        if i == 1
            y_min(1,n_cnt) = min(X1); y_max(1,n_cnt) = max(X1);
        end
        
        
        if min(X1)<y_min
            y_min(1,n_cnt) = min(X1);
        end
        
        if max(X1)>y_max
            y_max(1,n_cnt) = max(X1);
        end
    end
    X = [X X1];
    n_cnt = n_cnt + 1;
end

%% MCS
% [y] = uq_fourBranch(X,6);
y = uq_liquidHydrogenTank(X);

%% PLOT PDF
hist(y)

%% PLOT structural failure probability
[f,x] = ecdf(y); p_f = 1-f;
[B IX] = find(x<0)

figure;
plot(x,p_f,'b-','linewidth',2); hold on
plot(1,p_f(IX),'ro','linewidth',2);
plot(x(IX)*[1 1],[0 p_f(IX)],'r:','linewidth',2);
xlabel('Performance function');
ylabel('Probability of Exceedance');
set(gca,'yscale','log','fontsize',15,'fontweight','bold');
hold on; grid on
xlim([0 1.5]);

X_MC = X;
save MC_1E5.mat X