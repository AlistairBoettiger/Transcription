%%                              simulate_cellvar.m
% Alistair Boettiger                                   Date Begun: 11/23/10
%                                                   Last Modified: 11/23/10
%

% Simulate independent enhancer function in a whole embryo
% Take an embryo nuclear map, randomly assign a fraction of cells as on,
% randomly assign all cells a second time with same probability of being
% on.  Combine the two by finding all cells on in either the first or the
% second instance. 
%
% Plot all the results and color shift to yellow-cyan.

clear all; 

     
     fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Markov Modeling/Results/';

 folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed'; 
  base_emb = 'MP05xYW_30C_sna_y-full03_data.mat';
  % base_emb = 'SogP_C81_22C_LacZ_sog11_data.mat';
load([folder,'/',base_emb]); 
Ntot = max(H(:));
[h,w] = size(H); 





    Cyto = H.*Nuc_cntr;
    
load  compfull_pdist_simdata;   eqT = 1;
%load  compfull_pdist_simdata_X;   eqT = 5.25;
%load  compfull_pdist_simdata_X2;   eqT = 3.25;

%figure(1); clf; plot(FE);


% Construct density function for ER model results
    us = FE*1E7;
    FE_dist = [];
    step = 5;
    for k=1:length(FE)/step
        v =k*step;
        FE_dist = [FE_dist, v*ones(1,us(v))];
    end
    figure(2); clf; subplot(1,2,1); hist(FE_dist,0:step:length(FE));
title(['ER dist  mean_T: ',num2str(mean_E,3),'  std_T: ',num2str(std_E,3)]);
    
% Construct density function for IR model results
    us = FI*1E7;
    FI_dist = [];
    step = 5;
    for k=1:length(FI)/step
        v =k*step;
        FI_dist = [FI_dist, v*ones(1,us(v))];
    end
    figure(2); subplot(1,2,2); hist(FI_dist,0:step:length(FI));
title(['IR dist  mean_T: ',num2str(mean_I,3),'  std_T: ',num2str(std_I,3)]);

%%
time = lam/60; 
[jnk,Tmax] = find(time>120,1);



    % simulate draws under the ER Model
    NE = zeros(1,Ntot); 
    for k = 1:Ntot
        T = 0; N=0;
        while T<Tmax
        j = 1+round(rand*(length(FE_dist)-1)); 
        t_mRNA = FE_dist(j);
        N = N+1;
        T = t_mRNA + T;
        end
    NE(k) = N;
    end
    mean(NE)
 


    ER_dots = uint8(zeros(h,w)); 
        reg_data = regionprops(Cyto,'PixelIdxList');
            for k=1:Ntot;
                pixes = reg_data(k).PixelIdxList;      
                for m = 1:NE(k)
                    i = round(length(pixes)*rand);
                    if i>0
                    dot = pixes(i); 
                    end 
                    ER_dots(dot)=1;
                end
            end
            
                    
I = uint8(zeros(h,w,3));
I(:,:,3) = handles.In;%
I(:,:,2) = uint8(1-Nuc_line)*55;
I(:,:,1) = ER_dots*255;% 
     
Izoom_ER = I(400:500,1200:1400,:);
figure(1); clf; subplot(1,2,1); imshow(Izoom_ER);
title(['ER scheme:  ' ,num2str(mean(NE),2),' +/- ',num2str(std(NE),2),...
    '  transcripts'] );
    


% simulate draws under the IR Model
    NI = zeros(1,Ntot);
    for k = 1:Ntot
        T = 0; N=0;
        while T<Tmax*eqT
        j = 1+round(rand*(length(FI_dist)-1)); 
        t_mRNA = FI_dist(j);
        N = N+1;
        T = t_mRNA + T;
        end
       NI(k) = N;
    end      
mean(NI)

    xmin = 0; xmax = 40; bins = 20; 
    hs= figure(4); clf; subplot(1,2,1); 
    hist(NE,linspace(xmin,xmax,bins)); xlim([xmin,xmax]); 
    xlabel('number of transcripts'); ylabel('frequency'); 
    title('ER scheme'); 
    subplot(1,2,2); hist(NI,linspace(xmin,xmax,bins)); xlim([xmin,xmax]); 
     xlabel('number of transcripts'); ylabel('frequency'); 
     title('IR scheme');
    set(gcf,'color','w'); 
    saveas(hs,[fout,'/hist_Dpars_N.ai']);
  
    IR_dots = uint8(zeros(h,w)); 
        reg_data = regionprops(Cyto,'PixelIdxList');
            for k=1:Ntot;
                pixes = reg_data(k).PixelIdxList;      
                for m = 1:NI(k)
                    i = round(length(pixes)*rand);
                    if i>0
                    dot = pixes(i); 
                    end 
                    IR_dots(dot)=1;
                end
            end
            
                    
I = uint8(zeros(h,w,3));
I(:,:,3) = handles.In;%
I(:,:,2) = uint8(1-Nuc_line)*55;
I(:,:,1) = IR_dots*255;% 
     
Izoom_IR = I(400:500,1200:1400,:);
figure(1); subplot(1,2,2); imshow(Izoom_IR);
title(['IR scheme:  ' ,num2str(mean(NI),2),' +/- ',num2str(std(NI),2),...
    '  transcripts'] );
  





% title(['IR scheme, \mu= ',num2str(mean(NI),2),'  \sigma= ',num2str(std(NI),2),...
%     '  \sigma/\mu: ',num2str(std(NI)/mean(NI),2) ]);

%%



 
 
 % % IR model mu_t = 16, sigma_t = 12,   ER model mu_t = 5, sigma_t = 4; 
% 
% mu_t = 5; % mean time to transcription (minutes). 
% sigma_t = 4;
% mu_N = 5;  
% T = mu_t*mu_N; 
% mu_N = T/mu_t; % (number of transcripts)
% eta = sigma_t^2/mu_t;
% sigma_N = sqrt(sigma_t^2*T/mu_t^3);
% cov = sigma_N^2/mu_N^2;
% 
% mRNA = mu_N + sigma_N*randn(1E6,1); % poissrnd(mu,Ntot,1);
% mRNA(mRNA<0) = []; 
% mRNA = mRNA(1:Ntot); 
% disp(['\mu_N = ',num2str(mu_N,2), '   \sigma_N = ',num2str(sigma_N,2) ,'   \eta = ', num2str(eta,1)]);
% 
% h = figure(2); clf; hist(mRNA,0:30);
% saveas(h,[fout,'muT',num2str(mu_t,2),'sigmaT',num2str(sigma_t,2),'.ai']); 
% 
% mean(mRNA)
% std(mRNA)
% 
% 
% Cyto = H.*Nuc_cntr;
% 
%  cell_int = uint8(zeros(h,w)); 
%         reg_data = regionprops(Cyto,'PixelIdxList');
%             for k=1:Ntot;
%                 pixes = reg_data(k).PixelIdxList;      
%                 for m = 1:round(mRNA(k))
%                     i = round(length(pixes)*rand);
%                     if i>0
%                     dot = pixes(i); 
%                     end 
%                     cell_int(dot)=1;
%                 end
%             end
%             
%                     
% I = uint8(zeros(h,w,3));
% I(:,:,3) = handles.In;%
% I(:,:,2) = uint8(1-Nuc_line)*55;
% I(:,:,1) = cell_int*255;% 
% figure(1); clf; imshow(I);
%      
% Izoom = I(400:600,1200:1600,:);
% figure(1); clf; imshow(Izoom);
% 
% 
% 
%      imwrite(Izoom,[fout,'sim_','mu',num2str(mu_N,2),'_eta',num2str(eta,2),'.tif'],'tif');
%      
% C = [2,2,0;
%     1,1,0;
%     0,.61,.61];
% T = [0,1;
%     0,1;
%     0,1];
% f = [1,0];
% Iz = im_recolor(Izoom,C,T,f) ;
%  figure(1); clf; imshow(Iz);
%  
