%%                          multiloop_pdist.m
%
%
% Alistair Boettiger                                   Date Begun: 11/24/10
% Levine Lab                                Functionally Complete: 11/24/10
%                                                   Last Modified: 11/24/10
%

%% Description
%   Numerically Invert Laplace Transforms for a chosen set of model
%   parameters for the Full ER and IR models 

%% Modified


clear all;
load Model2b_pdist;

% Specific Model Constants
kab = sym('kab','real');
kba = sym('kba','real');
k12 = sym('k12','real'); 
k21 = sym('k21','real'); 
k23 = sym('k23','real'); 
k24 = sym('k24','real'); 
k32 = sym('k32','real'); 
k35 = sym('k35','real');
k53 = sym('k53','real');
k54 = sym('k54','real');
k42 = sym('k42','real'); 
k45 = sym('k45','real'); 
k56 = sym('k56','real'); 
k65 = sym('k65','real'); 
k67 = sym('k67','real'); 
k78 = sym('k78','real'); 


%              1        2       3       4        5       6       7       8        9         10       11       12       13        14     15     16
pars =   [  k12,      k21,   k23,      k24,     k32,    k35,    k53,    k54,     k42,       k45,     k56,     k65,    k67,       k78,   kab,   kba];         
% vals = 10*[6*.0216, 6*.145, 6*.0216, 6*.0216, 6*.145, 6*.0216 6*.145,
% 6*.145,  6*.0216,  6*.0216,  3*.00159, .001,  3*.00159,  3*.00159,.001, .1, ];
   vals = 5*[.0216, .145, 2, 2, 2, 2, 2,  2,  2,  2,  .00159, .1,  2,  2,.01, 1 ];
%  vals = [.00216, .145, 2, 2, 2, 2, 2,  2,  2,  2,  .159, .1,  2,  2,.1, 1];  % extreme case.    
% vals = rand(16,1);

lam = linspace(10,1E4,1000);  % time interval (in seconds);
fI_18 = subs(vI(1),pars,vals); 
fE_18 = subs(vE(1),pars,vals); 

fI_58 = subs(vI(3),pars,vals); 
fE_58 = subs(vE(3),pars,vals); 
      
        
% Here we actually numerically invert the Laplace space solutions of the
% model.  
   FI_18 = invlap2(fI_18, lam'); 
   FE_18 = invlap2(fE_18, lam'); 

   FI_58 = invlap2(fI_58, lam'); 
   FE_58 = invlap2(fE_58, lam'); 
 
%  xlim([0,1E4/60]);

  dt = (max(lam)-min(lam)) / length(lam) ; 


  save comp_pdist_multloop;

 figure(1); clf; 
 plot(lam/60,FI_18,'b','LineWidth',2); hold on; 
 plot(lam/60,FE_18,'r','LineWidth',2);   
 
 plot(lam/60,FI_58,'b--','LineWidth',2); 
 plot(lam/60,FE_58,'r--','LineWidth',2);   
 
 clear vE vI fI_18 fI_58 fE_18 fE_58 Gen GI GE;
 
%% Build Density functions

load comp_pdist_multloop;
 clear vE vI fI_18 fI_58 fE_18 fE_58 Gen GI GE;

% Construct density function for ER 1->8 model results
    us = FE_18*1E7;
    FE18_dist = zeros(1,100000);
    FE18_dist = [];
    step = 5;
    for k=1:floor( length(FE_18)/step )
        v =k*step;
        FE18_dist = [FE18_dist, v*ones(1,round(us(v)))];
    end
    figure(2); clf; subplot(2,2,1); hist(FE18_dist,0:step:length(FE_18));
title(['ER dist']);
    
% Construct density function for IR 1->8 model results
    us = FI_18*1E7;
    FI18_dist = [];
    step = 5;
    for k=1:floor( length(FI_18)/step )
        v =k*step;
        FI18_dist = [FI18_dist, v*ones(1,round(us(v)))];
    end
 subplot(2,2,2); hist(FI18_dist,0:step:length(FI_18));
title(['IR dist  ']);

% Construct density function for ER 5->8 model results
    us = FE_58*1E7;
    FE58_dist = [];
    step = 5;
    for k=1:floor( length(FE_58)/step )
        v =k*step;
        FE58_dist = [FE58_dist, v*ones(1,round(us(v)))];
    end
subplot(2,2,3); hist(FE58_dist,0:step:length(FE_58));
title(['ER dist']);
    
% Construct density function for IR 5->8 model results
    us = FI_58*1E7;
    FI58_dist = [];
    step = 5;
    for k=1:floor( length(FI_58)/step )
        v =k*step;
        FI58_dist = [FI58_dist, v*ones(1,round(us(v)))];
    end
subplot(2,2,4); hist(FI58_dist,0:step:length(FI_58));
title(['IR dist  ']);

%%
     fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Markov Modeling/Results/';
 folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed'; 
  base_emb = 'MP05xYW_30C_sna_y-full03_data.mat';
  % base_emb = 'SogP_C81_22C_LacZ_sog11_data.mat';
load([folder,'/',base_emb]); 
Ntot = max(H(:));
[h,w] = size(H); 
    Cyto = H.*Nuc_cntr;

%% Draw from density functions 

Tmax = 600; b = .02; % prob return to state 1

    % simulate draws under the ER Model
    NE = zeros(1,Ntot); 
    for k = 1:Ntot % k cycles of Ntot different cells
        
        T = 0; N=0;
        % first draw is always under the 1->8 distribution
            j = 1+round(rand*(length(FE18_dist)-1)); 
            t_mRNA = FE18_dist(j);
            N = N+1;
            T = t_mRNA + T;
        % ~~~~~~~~
        
        while T<Tmax  % during expression interval
            % choose 5->8 or 1->8 distribution based on bernuolli draw b. 
            if rand < b
                bdist = FE18_dist;
            else
                bdist = FE58_dist;
            end
            
            j = 1+round(rand*(length(bdist)-1)); 
            t_mRNA = bdist(j);
            N = N+1;
            T = t_mRNA + T;
        end    
    NE(k) = N;  % number of transcripts made by cell k.  
    end
    mean(NE)
 
    
    
    % simulate draws under the IR Model
    NI = zeros(1,Ntot); 
    for k = 1:Ntot % k cycles of Ntot different cells
        
        T = 0; N=0;
        % first draw is always under the 1->8 distribution
            j = 1+round(rand*(length(FI18_dist)-1)); 
            t_mRNA = FI18_dist(j);
            N = N+1;
            T = t_mRNA + T;
        % ~~~~~~~~
        
        while T<Tmax  % during expression interval
            % choose 5->8 or 1->8 distribution based on bernuolli draw b. 
            if rand < b
                bdist = FI18_dist;
            else
                bdist = FI58_dist;
            end
            
            j = 1+round(rand*(length(bdist)-1)); 
            t_mRNA = bdist(j);
            N = N+1;
            T = t_mRNA + T;
        end    
    NI(k) = N;  % number of transcripts made by cell k.  
    end
    mean(NI)
 
    
    %% plotting 
    xmin = 0; xmax = 40; bins = 20; 
    hs= figure(4); clf; subplot(2,1,1); 
    hist(NE,linspace(xmin,xmax,bins)); xlim([xmin,xmax]); ylim([0,600]);
    xlabel('number of transcripts'); 
title(['ER dist  mean_N=',num2str(mean(NE),3),...
    '  std_N=',num2str(std(NE),3), '  b=',num2str(b),'  T=',num2str(Tmax)]);
    subplot(2,1,2); hist(NI,linspace(xmin,xmax,bins)); xlim([xmin,xmax]); 
     xlabel('number of transcripts');  ylim([0,600]);
title(['IR dist  mean_N=',num2str(mean(NI),3),'  std_N=',...
    num2str(std(NI),3), '  b=',num2str(b),'  T=',num2str(Tmax)]);
    set(gcf,'color','w');  
    saveas(hs,[fout,'/hist_Dpars_N.ai']);
    
    %%
    
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

% imwrite(Izoom_ER,[fout,'simER_n_','b',num2str(b),'T',num2str(T),'.tif'],'tif'); 
%   imwrite(Izoom_IR,[fout,'simIR_n_','b',num2str(b),'T',num2str(T),'.tif'],'tif'); 
 

  %%  Simulated Cells
  err_thresh = 1.5;       
  
    ER_dots2 = imdilate(ER_dots,strel('disk',1));        
  
    IR_excess = all_nucs(NI> err_thresh *mean(NI));
    IRX = ismember(H,IR_excess); 
    IR_lack = all_nucs( err_thresh *NI<mean(NI));
    IRL =  ismember(H,IR_lack); 
    
      ER_excess = all_nucs(NE> err_thresh *mean(NE));
    ERX = ismember(H,ER_excess); 
    ER_lack = all_nucs( err_thresh *NE<mean(NE));
    ERL =  ismember(H,ER_lack); 
    
   % figure(5); clf; imshow(IRX)
    
I = uint8(255*ones(h,w,3));
I(:,:,3) = I(:,:,3) - ER_dots2*255 - uint8(1-Nuc_line)*155-uint8(100*ERX) ;%
I(:,:,2) = I(:,:,2)  - ER_dots2*255 -uint8(1-Nuc_line)*155-uint8(100*ERX) -uint8(100*ERL);
I(:,:,1) = I(:,:,1) -uint8(1-Nuc_line)*155-uint8(100*ERL);%  
     
Izoom_ER = I(400:600,1200:1600,:);
figure(1); clf; subplot(1,2,1); imshow(Izoom_ER);
title(['ER scheme:  ' ,num2str(mean(NE),2),' +/- ',num2str(std(NE),2),...
    '  transcripts'] );

IR_dots2 = imdilate(IR_dots,strel('disk',1));        
I = uint8(255*ones(h,w,3));
I(:,:,3) = I(:,:,3) - IR_dots2*255 - uint8(1-Nuc_line)*155 -uint8(100*IRX) ;%
I(:,:,2) = I(:,:,2)  - IR_dots2*255 -uint8(1-Nuc_line)*155 -uint8(100*IRX) -uint8(100*IRL);
I(:,:,1) = I(:,:,1) -uint8(1-Nuc_line)*155 -uint8(100*IRL);% 
     
Izoom_IR = I(400:600,1200:1600,:);
figure(1); subplot(1,2,2); imshow(Izoom_IR);
title(['IR scheme:  ' ,num2str(mean(NI),2),' +/- ',num2str(std(NI),2),...
    '  transcripts'] );
%   
%   imwrite(Izoom_ER,[fout,'DX_ER_pred.tif'],'tif'); 
%   imwrite(Izoom_IR,[fout,'DX_IR_pred.tif'],'tif'); 

  