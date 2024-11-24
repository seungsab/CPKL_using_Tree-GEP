
clc, clear all, close all

noise_amp = 0;
% noise_amp = 0.005;
% noise_amp = 0.05;

fn_sv_num = [2 4 6];
for IND_fcnt_type = 1:3
    
    IND = 1:5;
%     switch IND_fcnt_type
%         case 1 % 2D // branin_fcnt
%             IND = 1:size(x0,2);            
%         case 2 % 5D // fried_fcnt
%             IND = 1:size(x0,2)-1;            
%         case 3 % 8D // detpep108d
%             IND = 1:size(x0,2);            
%         case 4 % 20D // detpep108d
%             IND = 1:size(x0,2);
%     end
    
    
    switch noise_amp
        case 0
            fn3_str = ['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fn_sv_num(IND_fcnt_type))];
            svfn3_str = ['Result3_ENG_CPK_GEP_GP' num2str(IND_fcnt_type)];
            load(fn3_str);
                        
            Optimal3 = Optimal2(1:10,IND);
            RMSE3 = RMSE2(1:10,IND);
            MAE3 = RMAE2(1:10,IND);
            x0 = x0(1:10,IND);
            save(svfn3_str,'fun','noise_amp','y_val','x_val','Optimal3','RMSE3','MAE3','x0');
            
            
            
            fn1_str = ['Single_Kernel_GP_ENG' num2str(fn_sv_num(IND_fcnt_type))];
            svfn1_str = ['Result1_ENG_GP' num2str(IND_fcnt_type)];
            load(fn1_str);
                        
            Optimal1 = Optimal1(1:10,IND);
            RMSE1{1} = RMSE1{1}(1:10,:); RMSE1{2} = RMSE1{2}(1:10,:);
            RMSE1{3} = RMSE1{3}(1:10,:); RMSE1{4} = RMSE1{4}(1:10,:);
            RMSE1{5} = RMSE1{5}(1:10,:);
            
            MAE1{1} = RMAE1{1}(1:10,:); MAE1{2} = RMAE1{2}(1:10,:);
            MAE1{3} = RMAE1{3}(1:10,:); MAE1{4} = RMAE1{4}(1:10,:);
            MAE1{5} = RMAE1{5}(1:10,:);
            x0 = x0(1:10,IND);
            save(svfn1_str,'fun','noise_amp','y_val','x_val','Optimal1','RMSE1','MAE1','x0');
            
            
            fn2_str = ['AUTO_COMPOSITIONAL_GPR' num2str(fn_sv_num(IND_fcnt_type))];
            svfn2_str = ['Result2_ENG_HCPK_GP' num2str(IND_fcnt_type)];
            load(fn2_str);
            
            Optimal2 = Optimal2(1:10,IND);
            RMSE2 = RMSE2(1:10,IND);
            MAE2 = RMAE2(1:10,IND);
            x0 = x0(1:10,IND);
            save(svfn2_str,'fun','noise_amp','y_val','x_val','Optimal2','RMSE2','MAE2','x0');
            
            
        otherwise
            temp = num2str(noise_amp);
            
            fn3_str = ['AUTO_COMPOSITIONAL_GEP_GPR' num2str(fn_sv_num(IND_fcnt_type)) '_Noise_' [temp(1) temp(3:end)]];
            svfn3_str = ['Result3_ENG_CPK_GEP_GP' num2str(IND_fcnt_type) '_Noise_' [temp(1) temp(3:end)]];
            load(fn3_str);
            
            Optimal3 = Optimal2(1:10,IND);
            RMSE3 = RMSE2(1:10,IND);
            MAE3 = RMAE2(1:10,IND);
            x0 = x0(1:10,IND);
            save(svfn3_str,'fun','noise_amp','y_val','x_val','Optimal3','RMSE3','MAE3','x0');
            
            
            fn1_str = ['Single_Kernel_GP_ENG' num2str(fn_sv_num(IND_fcnt_type)) '_Noise_' [temp(1) temp(3:end)]];
            svfn1_str = ['Result1_ENG_GP' num2str(IND_fcnt_type) '_Noise_' [temp(1) temp(3:end)]];
            load(fn1_str);
            
            Optimal1 = Optimal1(1:10,IND);
            RMSE1{1} = RMSE1{1}(1:10,:); RMSE1{2} = RMSE1{2}(1:10,:);
            RMSE1{3} = RMSE1{3}(1:10,:); RMSE1{4} = RMSE1{4}(1:10,:);
            RMSE1{5} = RMSE1{5}(1:10,:);
            
            MAE1{1} = RMAE1{1}(1:10,:); MAE1{2} = RMAE1{2}(1:10,:);
            MAE1{3} = RMAE1{3}(1:10,:); MAE1{4} = RMAE1{4}(1:10,:);
            MAE1{5} = RMAE1{5}(1:10,:);
            x0 = x0(1:10,IND);
            save(svfn1_str,'fun','noise_amp','y_val','x_val','Optimal1','RMSE1','MAE1','x0');
            
            fn2_str = ['AUTO_COMPOSITIONAL_GPR' num2str(fn_sv_num(IND_fcnt_type)) '_Noise_' [temp(1) temp(3:end)]];
            svfn2_str = ['Result2_ENG_HCPK_GP' num2str(IND_fcnt_type) '_Noise_' [temp(1) temp(3:end)]];
            load(fn2_str);
            
            Optimal2 = Optimal2(1:10,IND);
            RMSE2 = RMSE2(1:10,IND);
            MAE2 = RMAE2(1:10,IND);
            x0 = x0(1:10,IND);
            save(svfn2_str,'fun','noise_amp','y_val','x_val','Optimal2','RMSE2','MAE2','x0');
            
    end
    
end