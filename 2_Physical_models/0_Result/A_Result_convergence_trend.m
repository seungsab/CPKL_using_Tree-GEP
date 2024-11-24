
%% INITIALIZATION
clc, clear all, close all

fn_size = 15;

%% LOAD FN
n_FCNT = 1; 
noise_amp = 0;
IND_fcnt_type = [2 4 6];

color_type = [
    0.00 0.45 0.74
    0.85 0.33 0.10
    0.93 0.69 0.13
    0.49 0.18 0.56
    0.47 0.67 0.19
    0.30,0.75,0.93
    ];

for i_fn = IND_fcnt_type
    switch noise_amp
        case 0
            try
                load(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(i_fn)]); RMSE3=RMSE2; RMAE3=RMAE2;
            end
            load(['AUTO_COMPOSITIONAL_GPR' num2str(i_fn)]);
            load(['Single_Kernel_GP_ENG' num2str(i_fn)]);
        otherwise
            temp = num2str(noise_amp);
            try
                load(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(i_fn) '_Noise_' [temp(1) temp(3:end)]]);
                RMSE3=RMSE2; RMAE3=RMAE2;
            end
            load(['AUTO_COMPOSITIONAL_GPR' num2str(i_fn) '_Noise_' [temp(1) temp(3:end)]]);
            load(['Single_Kernel_GP_ENG' num2str(i_fn) '_Noise_' [temp(1) temp(3:end)]]);
    end
        
    %% PLOT RESULT
    n_count = 1; RMSE =[]; RMAE=[]; group=[]; Nx = [];
    
    %% RMSE history of Validation data
    for n_type = 1:4
        for i = 1:size(x0,2)
            max_RMSE(n_type,i) = max(RMSE1{i}(:,n_type));
            min_RMSE(n_type,i) = min(RMSE1{i}(:,n_type));
            mean_RMSE(n_type,i) = mean(RMSE1{i}(:,n_type));
            
            max_RMAE(n_type,i) = max(RMAE1{i}(:,n_type));
            min_RMAE(n_type,i) = min(RMAE1{i}(:,n_type));
            mean_RMAE(n_type,i) = mean(RMAE1{i}(:,n_type));
            
            Nx(1,i) = size(x0{1,i},1);
        end
    end
end

max_RMSE(5,:) = max(RMSE2); min_RMSE(5,:) = min(RMSE2); mean_RMSE(5,:) = mean(RMSE2);
max_RMAE(5,:) = max(RMAE2); min_RMAE(5,:) = min(RMAE2); mean_RMAE(5,:) = mean(RMAE2);

max_RMSE(6,:) = max(RMSE3); min_RMSE(6,:) = min(RMSE3); mean_RMSE(6,:) = mean(RMSE3);
max_RMAE(6,:) = max(RMAE3); min_RMAE(6,:) = min(RMAE3); mean_RMAE(6,:) = mean(RMAE3);

for i_fn = 1:size(max_RMSE,1)
    figure(1); hold on
    fill([Nx flip(Nx)],[max_RMSE(n_type,:) flip(min_RMSE(n_type,:))],color_type(n_type,:),'FaceAlpha',0.2,'EdgeColor','w');
    a(i_fn)=plot(Nx,mean_RMSE(i_fn,:),'-','color',color_type(i_fn,:),'linewidth',5);
    
    figure(2); hold on
    a(i_fn)=plot(Nx,mean_RMAE(i_fn,:),'-','color',color_type(i_fn,:),'linewidth',5);
    fill([Nx flip(Nx)],[max_RMAE(n_type,:) flip(min_RMAE(n_type,:))],color_type(n_type,:),'FaceAlpha',0.2,'EdgeColor','w');
end

figure(1);
xlabel('# Sample','fontsize',fn_size,'fontweight','bold');
title('Normalized Root Mean Square Error (RMSE)','fontsize',15,'fontweight','bold')
set(gca,'yscale','log','fontsize',fn_size,'fontweight','bold');
set(gcf,'position',[576 427 791 472]);
legend(a, 'SEard', 'RQard', 'Maternard3','Maternard5','CKL','Location','sw');
axis tight; grid on

figure(2);
xlabel('# Sample','fontsize',fn_size,'fontweight','bold');
set(gca,'yscale','log','fontsize',fn_size,'fontweight','bold');
set(gcf,'position',[576 427 791 472]);
title('Normalized Mean Absolute Error (MAE)','fontsize',15,'fontweight','bold')
legend(a, 'SEard', 'RQard', 'Maternard3','Maternard5','CKL','Location','sw');
axis tight; grid on
