
%% INITIALIZATION
clc, clear all, close all

color_type = [
    0.00 0.45 0.74
    0.85 0.33 0.10
    0.93 0.69 0.13
    0.49 0.18 0.56
    0.47 0.67 0.19
    0.30,0.75,0.93
    ];

%% LOAD DATA
n_FCNT = 1;
% noise_amp = 0;
% noise_amp = 0.003;
% noise_amp = 0.005;
% noise_amp = 0.01;
% noise_amp = 0.03;
% noise_amp = 0.05;
% noise_amp = 0.1;

fn_size = 20;

noise_amp = [0 0.005 0.05];

fcnt_type = 2;
for ii0 = noise_amp
    temp = num2str(ii0);
    switch ii0
        case 0
            load(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type)]); noise = 0;
        otherwise
            load(['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fcnt_type) '_Noise_' [temp(1) temp(3:end)]]);
    end
    
    figure(n_FCNT)
    a = plot(y,y+noise,'o','markersize',4,'markerfacecolor',color_type(n_FCNT,:),'color',color_type(n_FCNT,:));
    hold on; axis tight
    plot([2 10],[2 10],'k:','linewidth',1);
    axis([2 10 2 10]);
    set(gca,'fontsize',fn_size,'fontweight','bold');
    xlabel('y','fontsize',fn_size,'fontweight','bold');
    ylabel('y+\beta*RMS(y)*N(0,1)','fontsize',fn_size,'fontweight','bold');
    title(['Noise amplitude (\beta): ' temp]);
    grid on
    
    n_FCNT = n_FCNT + 1;
end
