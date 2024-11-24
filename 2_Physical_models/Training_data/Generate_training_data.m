

n_MC = 30;

%% Define Test data
for dim = [7 8 10]  
    x0 = {}; x_val = [];
        switch dim
            case 1 % 1D
                N = 5:5:30;
                n_val = 300;
            case 2 % 2D
                N = 10:10:50;
                n_val = 500;
                
            case 3 % 3D
                N = 20:10:60;
                n_val = 600;
                
            case 5 % 5D
                N = 50:10:100;
                n_val = 1000;
                
            case 7% 7D
                N = 120:20:200;
                n_val = 3000;
                
            case 8 % 8D
                N = 120:20:200;
                n_val = 3000;
                
            case 10 % 10D
                N = 170:20:250;
                n_val = 5000;
        end
        
        for i=1:size(N,2)
            for j=1:n_MC
                x0{j,i} = lhsdesign(N(i),dim,'iterations',100);
            end
        end
        
        x_val = lhsdesign(n_val,dim,'iterations',100);
        
        save(['Training_' num2str(dim) 'D.mat'],'x0','x_val');
end