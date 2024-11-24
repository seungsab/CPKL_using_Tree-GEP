
%% INITIALIZATION
clc, clear all
%  close all

%% LOAD DATA
n_FCNT = 1;
% noise_amp = 0;
% noise_amp = 0.003;
% noise_amp = 0.005;
% noise_amp = 0.01;
% noise_amp = 0.03;
% noise_amp = 0.05;
% noise_amp = 0.1;

% noise_amp = 0.005;
noise_amp = 0.05;

IND_fcnt_type = [2 4 6];

for fcnt_type0 = IND_fcnt_type
    switch noise_amp
        case 0
            try
                load(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type0)]); RMSE3=RMSE2; RMAE3=RMAE2;
            end
            load(['AUTO_COMPOSITIONAL_GPR' num2str(fcnt_type0)]);
            load(['Single_Kernel_GP_ENG' num2str(fcnt_type0)]);
        otherwise
            temp = num2str(noise_amp);
            try
                load(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type0) '_Noise_' [temp(1) temp(3:end)]]);
                RMSE3=RMSE2; RMAE3=RMAE2;
            end
            load(['AUTO_COMPOSITIONAL_GPR' num2str(fcnt_type0) '_Noise_' [temp(1) temp(3:end)]]);
            load(['Single_Kernel_GP_ENG' num2str(fcnt_type0) '_Noise_' [temp(1) temp(3:end)]]);
    end
    
    
    IND = 1:size(x0,2);
    
    %% PLOT RESULT
    n_count = 1; RMSE =[]; RMAE=[]; group=[]; IND_colm = 1:10;
    mu_RMSE = []; mu_RMAE = []; std_RMSE = []; std_RMAE = [];
    max_RMSE = []; max_RMAE = []; min_RMSE = []; min_RMAE = [];
    legend_str={};
    for i=IND
        legend_str{1,n_count} = num2str(size(x0{1,i},1));
        try
            RMSE=[RMSE RMSE1{i}(IND_colm,1)' RMSE1{i}(IND_colm,2)'...
                RMSE1{i}(IND_colm,3)' RMSE1{i}(IND_colm,4)'...
                RMSE2(IND_colm,i)' RMSE3(IND_colm,i)'];
            RMAE=[RMAE RMAE1{i}(IND_colm,1)' RMAE1{i}(IND_colm,2)'...
                RMAE1{i}(IND_colm,3)' RMAE1{i}(IND_colm,4)'...
                RMAE2(IND_colm,i)' RMAE3(IND_colm,i)'];
            % LABELING FOR GROUPED BOXPLOT
            group=[group (6*(n_count-1)+1)*ones(1,10) (6*(n_count-1)+2)*ones(1,10) (6*(n_count-1)+3)*ones(1,10) ...
                (6*(n_count-1)+4)*ones(1,10) (6*(n_count-1)+5)*ones(1,10) (6*(n_count-1)+6)*ones(1,10)];
            xline = 6;
        catch
            RMSE=[RMSE RMSE1{i}(IND_colm,1)' RMSE1{i}(IND_colm,2)'...
                RMSE1{i}(IND_colm,3)' RMSE1{i}(IND_colm,4)' RMSE2(IND_colm,i)'];
            RMAE=[RMAE RMAE1{i}(IND_colm,1)' RMAE1{i}(IND_colm,2)'...
                RMAE1{i}(IND_colm,3)' RMAE1{i}(IND_colm,4)' RMAE2(IND_colm,i)'];
            % LABELING FOR GROUPED BOXPLOT
            group=[group (5*(n_count-1)+1)*ones(1,10) (5*(n_count-1)+2)*ones(1,10) (5*(n_count-1)+3)*ones(1,10) ...
                (5*(n_count-1)+4)*ones(1,10) (5*(n_count-1)+5)*ones(1,10)];
            xline = 5;
        end
        
        temp1 = [RMSE1{i}(IND_colm,1) RMSE1{i}(IND_colm,2) RMSE1{i}(IND_colm,3) RMSE1{i}(IND_colm,4) RMSE2(IND_colm,i) RMSE3(IND_colm,i)];
        temp2 = [RMAE1{i}(IND_colm,1) RMAE1{i}(IND_colm,2) RMAE1{i}(IND_colm,3) RMAE1{i}(IND_colm,4) RMAE2(IND_colm,i) RMAE3(IND_colm,i)];
        
        mu_RMSE = [mu_RMSE mean(temp1)']; mu_RMAE = [mu_RMAE mean(temp2)'];
        std_RMSE = [std_RMSE std(temp1)']; std_RMAE = [std_RMAE std(temp2)'];
        max_RMSE = [max_RMSE max(temp1)']; max_RMAE = [max_RMAE max(temp2)'];
        min_RMSE = [min_RMSE min(temp1)']; min_RMAE = [min_RMAE min(temp2)'];
        n_count = n_count + 1;
        
    end
    
    %%
    RMSE = RMSE/min(RMSE); RMAE = RMAE/min(RMAE);
    
    %     figure(1); set(gcf,'name','RMSE','position',[508.3 63.7 1973.3 1290.7]);
    %     subplot(5,2,2*n_FCNT-1)
    
    figure(1); set(gcf,'name','RMAE','position',[631.7 174.3 1048.7 1118.7]);
    subplot(size(IND_fcnt_type,2),1,n_FCNT)
    GROUPED_BOXPLOT(RMSE,group,legend_str); %grid on
    set(gca,'yscale','log'); hold on;
    plot(xline+1*[1 1],[min(RMSE) max(RMSE)],'k-');
    plot(2*(xline+1)*[1 1],[min(RMSE) max(RMSE)],'k-');
    plot(3*(xline+1)*[1 1],[min(RMSE) max(RMSE)],'k-');
    plot(4*(xline+1)*[1 1],[min(RMSE) max(RMSE)],'k-');
    c = get(gca, 'Children');
    if n_FCNT == 1
%         hleg1 = legend(c(5:10), 'SE', 'RQ', 'ME3','ME5','CKL-HKL','CKL-GEP','Location','best');
        title('Normalized Root Mean Square Error (RMSE)','fontsize',15,'fontweight','bold')
    elseif n_FCNT == 5
        xlabel('# of Training Samples');
    end
    switch fcnt_type
        case 2
            ylabel('OTL Circuit model');
        case 4
            ylabel('Borehole model');
        case 6
            ylabel('Wing weight model');
    end
    
    set(gca,'fontsize',15,'fontweight','bold'); grid on
    axis tight
    
    %     subplot(5,2,2*n_FCNT)
    figure(2); set(gcf,'name','RMAE','position',[631.7 174.3 1048.7 1118.7]);
    subplot(size(IND_fcnt_type,2),1,n_FCNT)
    GROUPED_BOXPLOT(RMAE,group,legend_str); %grid on
    set(gca,'yscale','log'); hold on;
    plot(xline+1*[1 1],[min(RMAE) max(RMAE)],'k-');
    plot(2*(xline+1)*[1 1],[min(RMAE) max(RMAE)],'k-');
    plot(3*(xline+1)*[1 1],[min(RMAE) max(RMAE)],'k-');
    plot(4*(xline+1)*[1 1],[min(RMAE) max(RMAE)],'k-');
    c = get(gca, 'Children');
    if n_FCNT == 1
        title('Normalized Mean Absolute Error (MAE)','fontsize',15,'fontweight','bold')
    elseif n_FCNT == 5
        xlabel('# of Training Samples');
    end
    switch fcnt_type
        case 2
            ylabel('OTL Circuit model');
        case 4
            ylabel('Borehole model');
        case 6
            ylabel('Wing weight model');
    end
    set(gca,'fontsize',15,'fontweight','bold'); grid on
    axis tight
    
    n_FCNT = n_FCNT + 1;
end