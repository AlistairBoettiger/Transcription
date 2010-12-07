
%% Plotting Commands 

fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Markov Modeling/Results/'

%% Original Version
load Markov_simp_dist; 


% Plotting parameters 
xmin = -5;  % min log difference
xmax = 20;  % max log difference
offset = .1;  % little buffer on the axis
xs = linspace(xmin,xmax,50); % xaxis and number of bins
F = 10; % Fontesize

 
% convert ratios of parameters into 
 m = log2(dM);  
 s = log2(dV);
 n = log2(dN); 
 

% Figure comparing ratio of values over all parameter space. 
f1 = figure(1); clf; 
subplot(1,3,1);
mp = m(m>0); 
mn = m(m<0); 
hist(mn,xs);
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r'); hold on;
 hist(mp,xs); 
  xlim([xmin-offset,xmax+offset]);  
 xlabel('log_{2}(\mu_{IR}/\mu_{ER})'); set(gca,'YTickLabel',' ');
 ylabel('frequency');  set(gca,'FontSize',F);
title('Mean Expression Speed, \mu','FontSize',F);

subplot(1,3,2);
sp = s(s>0);
sn = s(s<0);
 hist(sn,xs); 
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r'); hold on;
 hist(sp,xs); 
  xlim([xmin-offset,xmax+offset]); 
xlabel('log_{2}(\sigma^2_{IR}/\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing, \sigma^2','FontSize',F); set(gca,'FontSize',F);


subplot(1,3,3);
np = n(n>0);
nn = n(n<0);
 hist(nn,xs); hold on;
     h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r');
 hist(np,xs); 
  xlim([xmin-offset,xmax+offset]); 
xlabel('log_{2}(\eta_{IR}/\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number, \eta','FontSize',F); 
set(gcf,'color','w'); set(gca,'FontSize',F);

 
 

 % Plot parameter skews when IR is smaller.  
 var_names = {'K12','K21','K23','K34','Kab','Kba'};
 xx = linspace(0,1,6);
 figure(3); clf; set(gcf,'color','w'); 
 kk = 0;
  IR_better = log2(dM)<-.5;
 for k=1:6
     kk = kk+1;
    subplot(3,6,kk); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
    title(['\mu ', var_names{k}]); xlim([-.05,1.05]);
 end
 IR_better = log2(dV)<-.5;
  for k=1:6
     kk = kk+1;
    subplot(3,6,kk); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
    title(['\sigma^2 ', var_names{k}]); xlim([-.05,1.05]);
  end
  IR_better = log2(dN)<-.3;
   for k=1:6
     kk = kk+1;
    subplot(3,6,kk); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
    title(['\eta ', var_names{k}]); xlim([-.05,1.05]);
   end
 
%% Regulate k32
%
%  
% % Model Definition:
% % Simple Initiation Regulated Model
% GI{1} = zeros(3);
% GI{2} =[[-kab,kab,0,0,0];
%         [kba,-k12-kba,k12,0,0];
%         [0,k21,-k21-k23,k23,0];
%         [0,0,0,-k34,k34];
%         [0,0,0,0,0 ]];  
% 
% % Revised Elongation model: regulat k23 instead of k34
% % Implies state 3A is inaccesable  
% GE{1} = zeros(3);
% %            1A 2A 3A  1B  2B  3B   4 
% GE{2} = [[-k12-kab, k12,    kab,    0,  0,  0];       % 1A
%         [k21,   -k21-kab,   0,    kab,    0,  0];  % 2A
%         [kba,0,-kba-k12,k12,0,0];        % 1B
%         [0,kba,k21,-kba-k21-k23,k23,0];  % 2B
%         [0,0,0,0,-k34,k34] ;       % 3B
%         [0,0,0,0,0,0]];                  % 4   


% load data for plotting:
load plotdata_regk23



% Plotting parameters 
xmin = -3;  % min log difference
xmax = 3;  % max log difference
offset = .1;  % little buffer on the axis
xs = linspace(xmin,xmax,50); % xaxis and number of bins
F = 10; % Fontesize

 
% convert ratios of parameters into 
 m = log2(dM);  
 s = log2(dV);
 n = log2(dN); 
 

% Figure comparing ratio of values over all parameter space. 
f1 = figure(1); clf; 
subplot(1,3,1);
mp = m(m>0); 
mn = m(m<0); 
hist(mn,xs);
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r'); hold on;
 hist(mp,xs); 
  xlim([xmin-offset,xmax+offset]);  
 xlabel('log_{2}(\mu_{IR}/\mu_{ER})'); set(gca,'YTickLabel',' ');
 ylabel('frequency');  set(gca,'FontSize',F);
title('Mean Expression Speed, \mu','FontSize',F);

subplot(1,3,2);
sp = s(s>0);
sn = s(s<0);
 hist(sn,xs); 
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r'); hold on;
 hist(sp,xs); 
  xlim([xmin-offset,xmax+offset]); 
xlabel('log_{2}(\sigma^2_{IR}/\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing, \sigma^2','FontSize',F); set(gca,'FontSize',F);


subplot(1,3,3);
np = n(n>0);
nn = n(n<0);
 hist(nn,xs); hold on;
     h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r');
 hist(np,xs); 
  xlim([xmin-offset,xmax+offset]); 
xlabel('log_{2}(\eta_{IR}/\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number, \eta','FontSize',F); 
set(gcf,'color','w'); set(gca,'FontSize',F);

 
 

 % Plot parameter skews when IR is smaller.  
 var_names = {'K12','K21','K23','K34','Kab','Kba'};
 xx = linspace(0,1,6);
 figure(3); clf; set(gcf,'color','w'); 
 kk = 0;
  IR_better = log2(dM)<-.5;
 for k=1:6
     kk = kk+1;
    subplot(3,6,kk); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
    title(['\mu ', var_names{k}]); xlim([-.05,1.05]);
 end
 IR_better = log2(dV)<-.5;
  for k=1:6
     kk = kk+1;
    subplot(3,6,kk); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
    title(['\sigma^2 ', var_names{k}]); xlim([-.05,1.05]);
  end
  IR_better = log2(dN)<-.3;
   for k=1:6
     kk = kk+1;
    subplot(3,6,kk); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
    title(['\eta ', var_names{k}]); xlim([-.05,1.05]);
   end
 
   
   
   %%  Add k32
%    % Simple Initiation Regulated Model
% GI{1} = zeros(3);
% GI{2} =[[-kab,kab,0,0,0];
%         [kba,-k12-kba,k12,0,0];
%         [0,k21,-k21-k23,k23,0];
%         [0,0,k32,-k34-k32,k34];
%         [0,0,0,0,0 ]];
%    
%    [m1I, m2I] = SeriesDecomp(GI); 
% 
% % Simple Elongation model
% GE{1} = zeros(3);
% %            1A 2A 3A  1B  2B  3B   4 
% GE{2} = [[-k12-kab,k12,0,kab,0,0,0];       % 1A
%         [k21,-k21-kab-k23,k23,0,kab,0,0];  % 2A
%         [0,k32,-kab-k32,0,0,kab,0 ];             % 3A
%         [kba,0,0,-kba-k12,k12,0,0];        % 1B
%         [0,kba,0,k21,-kba-k21-k23,k23,0];  % 2B
%         [0,0,kba,0,k32,-k32-k34-kba,k34] ;       % 3B
%         [0,0,0,0,0,0,0]];   
   
  clear all; 
   load plotdata_addk32



% Plotting parameters 
xmin = -3;  % min log difference
xmax = 3;  % max log difference
offset = .1;  % little buffer on the axis
xs = linspace(xmin,xmax,50); % xaxis and number of bins
F = 10; % Fontesize

 
% convert ratios of parameters into 
 m = log2(dM);  
 s = log2(dV);
 n = log2(dN); 
 

% Figure comparing ratio of values over all parameter space. 
figure(1); clf; 
subplot(1,3,1);
mp = m(m>0); 
mn = m(m<0); 
hist(mn,xs);
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r'); hold on;
 hist(mp,xs); 
  xlim([xmin-offset,xmax+offset]);  
 xlabel('log_{2}(\mu_{IR}/\mu_{ER})'); set(gca,'YTickLabel',' ');
 ylabel('frequency');  set(gca,'FontSize',F);
title('Mean Expression Speed, \mu','FontSize',F);

subplot(1,3,2);
sp = s(s>0);
sn = s(s<0);
 hist(sn,xs); 
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r'); hold on;
 hist(sp,xs); 
  xlim([xmin-offset,xmax+offset]); 
xlabel('log_{2}(\sigma^2_{IR}/\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing, \sigma^2','FontSize',F); set(gca,'FontSize',F);


subplot(1,3,3);
np = n(n>0);
nn = n(n<0);
 hist(nn,xs); hold on;
     h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r');
 hist(np,xs); 
  xlim([xmin-offset,xmax+offset]); 
xlabel('log_{2}(\eta_{IR}/\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number, \eta','FontSize',F); 
set(gcf,'color','w'); set(gca,'FontSize',F);

 
 

 % Plot parameter skews when IR is smaller.  
 var_names = {'K12','K21','K23','K34','Kab','Kba'};
 xx = linspace(0,1,6);
 figure(3); clf; set(gcf,'color','w'); 
 kk = 0;
  IR_better = log2(dM)<-.5;
 for k=1:6
     kk = kk+1;
    subplot(3,6,kk); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
    title(['\mu ', var_names{k}]); xlim([-.05,1.05]);
 end
 IR_better = log2(dV)<-.5;
  for k=1:6
     kk = kk+1;
    subplot(3,6,kk); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
    title(['\sigma^2 ', var_names{k}]); xlim([-.05,1.05]);
  end
  IR_better = log2(dN)<-.3;
   for k=1:6
     kk = kk+1;
    subplot(3,6,kk); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
    title(['\eta ', var_names{k}]); xlim([-.05,1.05]);
   end
 
   
   