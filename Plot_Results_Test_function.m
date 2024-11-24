
function main
%% Compositional Kernel Learning using Tree-based Genetic Programming (CPKL_Tree-GEP) for Gaussian Proecess Regression
% Updated for Replication of results

% This code is for plotting the results in manuscript
% Mathematical test fucntion

% If you need source code for CPK using tree-GEP, feel free to contact me
% via e-mail (seungsab@gmail.com)
% Coded by SS (2019.08.21)

%% INITIALIZATION
clc, clear all, close all, warning('off','all')
addpath([pwd '\1_Test_functions']);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select the types for results
% Select function type
IND_fcnt_type = 1; % Test function #1
% IND_fcnt_type = 2; % Test function #2
% IND_fcnt_type = 3; % Test function #3
% IND_fcnt_type = 4; % Test function #4

% Select noise level
noise_amp = 0;
% noise_amp = 0.005;
% noise_amp = 0.05;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD RESULTS FROM GPs using single kernel, HCPK, and CPK-GEP
disp(['Function #' num2str(IND_fcnt_type) ' // Noise level:' num2str(noise_amp)]);

switch noise_amp
    case 0
        load(['Result3_CPK_GEP_GP' num2str(IND_fcnt_type)]);
        load(['Result2_HCPK_GP' num2str(IND_fcnt_type)]);
        load(['Result1_GP' num2str(IND_fcnt_type)]);
    otherwise
        temp = num2str(noise_amp);
        load(['Result3_CPK_GEP_GP' num2str(IND_fcnt_type) '_Noise_' [temp(1) temp(3:end)]]);
        load(['Result2_HCPK_GP' num2str(IND_fcnt_type) '_Noise_' [temp(1) temp(3:end)]]);
        load(['Result1_GP' num2str(IND_fcnt_type) '_Noise_' [temp(1) temp(3:end)]]);
end

%% MERGE ALL RESULTS (RMSE AND MAE)
n_count = 1; RMSE =[]; MAE=[]; group=[]; IND_colm = 1:10;
mu_RMSE = []; mu_MAE = []; std_RMSE = []; std_MAE = [];
max_RMSE = []; max_MAE = []; min_RMSE = []; min_MAE = [];
legend_str={};
for i=1:size(x0,2);
    legend_str{1,n_count} = num2str(size(x0{1,i},1));
    RMSE=[RMSE RMSE1{i}(IND_colm,1)' RMSE1{i}(IND_colm,2)'...
        RMSE1{i}(IND_colm,3)' RMSE1{i}(IND_colm,4)'...
        RMSE2(IND_colm,i)' RMSE3(IND_colm,i)'];
    MAE=[MAE MAE1{i}(IND_colm,1)' MAE1{i}(IND_colm,2)'...
        MAE1{i}(IND_colm,3)' MAE1{i}(IND_colm,4)'...
        MAE2(IND_colm,i)' MAE3(IND_colm,i)'];
    % LABELING FOR GROUPED BOXPLOT
    group=[group (6*(n_count-1)+1)*ones(1,10) (6*(n_count-1)+2)*ones(1,10) (6*(n_count-1)+3)*ones(1,10) ...
        (6*(n_count-1)+4)*ones(1,10) (6*(n_count-1)+5)*ones(1,10) (6*(n_count-1)+6)*ones(1,10)];
    xline = 6;
    
    temp1 = [RMSE1{i}(IND_colm,1) RMSE1{i}(IND_colm,2) RMSE1{i}(IND_colm,3) RMSE1{i}(IND_colm,4) RMSE2(IND_colm,i) RMSE3(IND_colm,i)];
    temp2 = [MAE1{i}(IND_colm,1) MAE1{i}(IND_colm,2) MAE1{i}(IND_colm,3) MAE1{i}(IND_colm,4) MAE2(IND_colm,i) MAE3(IND_colm,i)];
    
    mu_RMSE = [mu_RMSE mean(temp1)']; mu_MAE = [mu_MAE mean(temp2)'];
    std_RMSE = [std_RMSE std(temp1)']; std_MAE = [std_MAE std(temp2)'];
    max_RMSE = [max_RMSE max(temp1)']; max_MAE = [max_MAE max(temp2)'];
    min_RMSE = [min_RMSE min(temp1)']; min_MAE = [min_MAE min(temp2)'];
    n_count = n_count + 1;
end

%% PLOT RESULTS
% Normalize performance measures by smaller value
RMSE = RMSE/min(RMSE); MAE = MAE/min(MAE);

% PLOT RMSE
figure(1); set(gcf,'name','RMAE')
GROUPED_BOXPLOT(RMSE,group,legend_str); %grid on
set(gca,'yscale','log'); hold on;
plot(xline+1*[1 1],[min(RMSE) max(RMSE)],'k-');
plot(2*(xline+1)*[1 1],[min(RMSE) max(RMSE)],'k-');
plot(3*(xline+1)*[1 1],[min(RMSE) max(RMSE)],'k-');
plot(4*(xline+1)*[1 1],[min(RMSE) max(RMSE)],'k-');
c = get(gca, 'Children');
hleg1 = legend(c(5:10), 'SE', 'RQ', 'ME3','ME5','CKL-HKL','CKL-GEP','Location','best');
title('Normalized Root Mean Square Error (RMSE)','fontsize',15,'fontweight','bold')
xlabel('# of Training Samples');
ylabel(['Test function #' num2str(IND_fcnt_type)]);
set(gca,'fontsize',15,'fontweight','bold'); grid on

% PLOT MAE
figure(2); set(gcf,'name','RMAE')
GROUPED_BOXPLOT(MAE,group,legend_str); %grid on
set(gca,'yscale','log'); hold on;
plot(xline+1*[1 1],[min(MAE) max(MAE)],'k-');
plot(2*(xline+1)*[1 1],[min(MAE) max(MAE)],'k-');
plot(3*(xline+1)*[1 1],[min(MAE) max(MAE)],'k-');
plot(4*(xline+1)*[1 1],[min(MAE) max(MAE)],'k-');
c = get(gca, 'Children');
title('Normalized Mean Absolute Error (MAE)','fontsize',15,'fontweight','bold')
xlabel('# of Training Samples');
ylabel(['Test function #' num2str(IND_fcnt_type)]);
set(gca,'fontsize',15,'fontweight','bold'); grid on
axis tight

rmpath([pwd '\1_Test_functions']);
end

function GROUPED_BOXPLOT(x,group,legend_str)
n_OBJ=group(end)/size(legend_str,2);

positions = []; n_shift=0;
for i=1:max(group)/n_OBJ
    temp = [n_OBJ*(i-1)+1:n_OBJ*i] + n_shift;
    positions = [positions temp];
    n_shift = n_shift + 1;
end
boxplot(x,group, 'positions', positions);

mu_pos=[];  color=[];
for i=1:size(legend_str,2)
    mu_pos(1,i)=mean(positions(n_OBJ*(i-1)+1:n_OBJ*i));
    if n_OBJ == 5
        color = [color 'b', 'r', 'y', 'c', 'm'];
    else
        color = [color 'b', 'r', 'y', 'c', 'm', 'g'];
    end
    
end

set(gca,'xtick',mu_pos)
set(gca,'xticklabel',legend_str)
color=[color color];

h = findobj(gca,'Tag','Box');
for j=1:length(h)
%    patch(get(h(j),'XData'),get(h(j),'YData'),color(j));
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.55);
end

end