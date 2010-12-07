
%% Intuit eta
clear all; 



% 
fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Markov Modeling/Results/';
 folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed'; 
  base_emb = 'MP05xYW_30C_sna_y-full03_data.mat';
  % base_emb = 'SogP_C81_22C_LacZ_sog11_data.mat';
load([folder,'/',base_emb]); 


save Cell_data H Nuc_cntr Nuc_line all_nucs; 
clear all;
load Cell_data; 

    Ntot = max(H(:));
    [h,w] = size(H); 
   Cyto = H.*Nuc_cntr;

% mu_NE = 20;
% sigma_NE = 2.5;
% mu_NI = 18.5;
% sigma_NI = 11.7;


mu_NE = 8;
sigma_NE = 3;
mu_NI = 0;
sigma_NI = 10;

% mu_NE = 30;
% sigma_NE = 5;
% mu_NI = 14;
% sigma_NI = 30;


mRNA_E = mu_NE + sigma_NE*randn(1E6,1); % poissrnd(mu,Ntot,1);
mRNA_E(mRNA_E<0) = []; 
mRNA_E = mRNA_E(1:Ntot); 




mRNA_I = mu_NI + sigma_NI*randn(1E6,1); % poissrnd(mu,Ntot,1);
mRNA_I(mRNA_I<0) = []; 
mRNA_I = mRNA_I(1:Ntot); 



     % Histograms    (trivial here, just the truncated gaussians from above)   
                xmin = 0; xmax = 50; bins = 20; 
    hs= figure(4); clf; subplot(1,2,1); 
    hist(mRNA_E,linspace(xmin,xmax,bins)); xlim([xmin,xmax]); 
    xlabel('number of transcripts'); ylabel('frequency'); 
    title(['ER scheme,  mean=',num2str(mean(mRNA_E),2), '   std=' num2str(std(mRNA_E),2)]); 
    subplot(1,2,2); hist(mRNA_I,linspace(xmin,xmax,bins)); xlim([xmin,xmax]); 
     xlabel('number of transcripts'); ylabel('frequency'); 
     title(['IR scheme,  mean=',num2str(mean(mRNA_I),2), '   std=' num2str(std(mRNA_I),2)]); 
    set(gcf,'color','w'); 
    
    

%% Simulate transcripts placement under model
        ER_dots = uint8(zeros(h,w)); 
        reg_data = regionprops(Cyto,'PixelIdxList');
            for k=1:Ntot;
                pixes = reg_data(k).PixelIdxList;      
                for m = 1:mRNA_E(k)
                    i = round(length(pixes)*rand);
                    if i>0
                    dot = pixes(i); 
                    end 
                    ER_dots(dot)=1;
                end
            end
            

        IR_dots = uint8(zeros(h,w)); 
        reg_data = regionprops(Cyto,'PixelIdxList');
            for k=1:Ntot;
                pixes = reg_data(k).PixelIdxList;      
                for m = 1:mRNA_I(k)
                    i = round(length(pixes)*rand);
                    if i>0
                    dot = pixes(i); 
                    end 
                    IR_dots(dot)=1;
                end
            end
   

    
    %%  Simulated Cells
         
    ER_dots2 = imdilate(ER_dots,strel('disk',1));        
  
    IR_excess = all_nucs(mRNA_I>2*mean(mRNA_I));
    IRX = ismember(H,IR_excess); 
    IR_lack = all_nucs(2*mRNA_I<mean(mRNA_I));
    IRL =  ismember(H,IR_lack); 
    
      ER_excess = all_nucs(mRNA_E>2*mean(mRNA_E));
    ERX = ismember(H,ER_excess); 
    ER_lack = all_nucs(2*mRNA_E<mean(mRNA_E));
    ERL =  ismember(H,ER_lack); 
    
   % figure(5); clf; imshow(IRX)
    
I = uint8(255*ones(h,w,3));
I(:,:,3) = I(:,:,3) - ER_dots2*255 - uint8(1-Nuc_line)*155-uint8(100*ERX) ;%
I(:,:,2) = I(:,:,2)  - ER_dots2*255 -uint8(1-Nuc_line)*155-uint8(100*ERX) -uint8(100*ERL);
I(:,:,1) = I(:,:,1) -uint8(1-Nuc_line)*155-uint8(100*ERL);%  
     
Izoom_ER = I(400:500,1200:1400,:);
figure(1); clf; subplot(1,2,1); imshow(Izoom_ER);
title(['ER scheme:  ' ,num2str(mean(mRNA_E),2),' +/- ',num2str(std(mRNA_E),2),...
    '  transcripts'] );

IR_dots2 = imdilate(IR_dots,strel('disk',1));        
I = uint8(255*ones(h,w,3));
I(:,:,3) = I(:,:,3) - IR_dots2*255 - uint8(1-Nuc_line)*155 -uint8(100*IRX) ;%
I(:,:,2) = I(:,:,2)  - IR_dots2*255 -uint8(1-Nuc_line)*155 -uint8(100*IRX) -uint8(100*IRL);
I(:,:,1) = I(:,:,1) -uint8(1-Nuc_line)*155 -uint8(100*IRL);% 
     
Izoom_IR = I(400:500,1200:1400,:);
figure(1); subplot(1,2,2); imshow(Izoom_IR);
title(['IR scheme:  ' ,num2str(mean(mRNA_I),2),' +/- ',num2str(std(mRNA_I),2),...
    '  transcripts'] );
%   
%  imwrite(Izoom_ER,[fout,'mu',num2str(mu_NE,2),'sigma',num2str(sigma_NE,2),'.tif'],'tif'); 
%  imwrite(Izoom_IR,[fout,'mu',num2str(mu_NI,2),'sigma',num2str(sigma_NI,2),'.tif'],'tif'); 
