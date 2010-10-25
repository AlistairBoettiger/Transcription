 
%  Model Analysis 
% 
% Alistair Boettiger                                   Date Begun: 07/02/09
% Levine Lab                                        Last Modified: 08/05/09


clear all; 
 load modelFull 
%load modelFull3; % 10,000 x 16
% load modelFull2; % 5000 x 16, independent of modelFull.mat

% modelFull.mat contains
%  * m1En, m1In, for 1000 random k values, the mean transition time for
%  elongation regulated and initiation regulated models 
%   *  m2En, m2In, for 1000 random k values, 2nd moment transition time for
%  elongation regulated and initiation regulated models   
%  * 1000 16-vectors that are the gradients for the mean transition time for
%  elongation, initation, difference I-E, standard deviation of E, I,
%  difference, noise for E, I and difference,

%-------------------------------------------------------------%
%                      Sensitivity Analysis                   %
%-------------------------------------------------------------%
% Compute percent of variance explained by parameter, v(k)
% definition: v(k) = (partial q/partial k)^2 / (length grad q)^2
% where q is speed, variance, or noise, and k are the transition rates.  

% elongation mean
lmE = sum(grad_mEn.^2,2); % length
lmE = meshgrid(lmE,1:16); % duplicate in 16 columns to divide
vmE = grad_mEn.^2 ./ lmE'; % percent of variance explained

% iniation mean
lmI = sum(grad_mIn.^2,2); % length
lmI = meshgrid(lmI,1:16); % duplicate in 16 columns to divide
vmI = grad_mIn.^2 ./ lmI'; % percent of variance explained

% difference in means
ldm = sum(grad_dmn.^2,2); % length
ldm = meshgrid(ldm,1:16); % duplicate in 16 columns to divide
vdm = grad_dmn.^2 ./ ldm'; % percent of variance explained

% elongation variance
lsE = sum(grad_sEn.^2,2); % length
lsE = meshgrid(lsE,1:16); % duplicate in 16 columns to divide
vsE = grad_sEn.^2 ./ lsE'; % percent of variance explained

% initation variance
lsI = sum(grad_sIn.^2,2); % length
lsI = meshgrid(lsI,1:16); % duplicate in 16 columns to divide
vsI = grad_sIn.^2 ./ lsI'; % percent of variance explained

% difference in variance
lds = sum(grad_dsn.^2,2); % length
lds = meshgrid(lds,1:16); % duplicate in 16 columns to divide
vds = grad_dsn.^2 ./ lds'; % percent of variance explained

% elongation noise
lnE = sum(grad_nEn.^2,2); % length
lnE = meshgrid(lnE,1:16); % duplicate in 16 columns to divide
vnE = grad_nEn.^2 ./ lnE'; % percent of variance explained

% initation noise
lnI = sum(grad_nIn.^2,2); % length
lnI = meshgrid(lnI,1:16); % duplicate in 16 columns to divide
vnI = grad_nIn.^2 ./ lnI'; % percent of variance explained

% difference in noise
ldn = sum(grad_dnn.^2,2); % length
ldn = meshgrid(ldn,1:16); % duplicate in 16 columns to divide
vdn = grad_dnn.^2 ./ ldn'; % percent of variance explained


mean(vmE)
mean(vmI)
mean(vdm)

F=10;

% Plot comparisons of means
figure(2); clf; colormap(jet);
   subplot(3,1,1);
     histin = vmE; histin(histin==0)=NaN; 
     hist(log10(abs(histin)));      
     ylim([0,length(grad_mEn)]); 
     set(gca,'YtickLabel',str2mat('0','0.2','0.4','0.6','0.8','1'),'fontsize',F);
     ylabel('normalized frequency'); 
     xlabel('log_{10} (variance explained)') 
     title('ER model mean time to expression','FontSize',12);
   subplot(3,1,2);
      histin = vmI; histin(histin==0)=NaN; 
      hist(log10(abs(histin)));
      ylim([0,length(grad_mEn)]); 
      set(gca,'YtickLabel',str2mat('0','0.2','0.4','0.6','0.8','1'),'fontsize',F);
      ylabel('normalized frequency'); 
      xlabel('log_{10} (variance explained)');
      title('IR model mean time to expression','FontSize',12); 
      legend('K01','K10','K12','K21','K23','K24','K32','K35','K42','K45',...
     'K56','K65','K67','K78','Kab','Kba','Location','West'); 
 subplot(3,1,3);
    histin = vdm;  histin(histin==0)=NaN; 
    hist(log10(abs(histin)));      
    ylim([0,length(grad_mEn)]); 
    set(gca,'YtickLabel',str2mat('0','0.2','0.4','0.6','0.8','1'),'fontsize',F);
    ylabel('normalized frequency'); 
    xlabel('log_{10} (variance explained)') 
    title('difference in mean time to expression','FontSize',12); 

set(gcf,'color','white');

    
    
 % Plot comparisons of variances   
    figure(3); clf; colormap(jet);
  histin = vsE; histin(histin==0)=NaN; 
   subplot(3,1,1);
   hist(log10(abs(histin)));      
   ylim([0,length(grad_mEn)]); 
   set(gca,'YtickLabel',str2mat('0','0.2','0.4','0.6','0.8','1'),'fontsize',F);
    ylabel('normalized frequency'); 
    xlabel('log_{10} (variance explained)') 
       title('ER model variance in time to expression','FontSize',12);

    subplot(3,1,2);
     histin = vsI; histin(histin==0)=NaN; 
     hist(log10(abs(histin)));
      ylim([0,length(grad_mEn)]); 
      set(gca,'YtickLabel',str2mat('0','0.2','0.4','0.6','0.8','1'),'fontsize',F);
      ylabel('normalized frequency'); 
      xlabel('log_{10} (variance explained)');
      title('IR model variance in time to expression','FontSize',12); 
      legend('K01','K10','K12','K21','K23','K24','K32','K35','K42','K45',...
     'K56','K65','K67','K78','Kab','Kba','Location','West');
  
 subplot(3,1,3);
     histin = vds;  histin(histin==0)=NaN; 
   hist(log10(abs(histin)));      
   ylim([0,length(grad_mEn)]); 
   set(gca,'YtickLabel',str2mat('0','0.2','0.4','0.6','0.8','1'),'fontsize',F);
    ylabel('normalized frequency'); 
    xlabel('log_{10} (variance explained)') 
    title('difference in varaince in time to expression','FontSize',12); 

set(gcf,'color','white');

    
    
% Plot comparisons of noise    
    figure(4); clf; colormap(jet);
  histin = vnE; histin(histin==0)=NaN; 
   subplot(3,1,1);
   hist(log10(abs(histin)));      
   ylim([0,length(grad_mEn)]); 
   set(gca,'YtickLabel',str2mat('0','0.2','0.4','0.6','0.8','1'),'fontsize',F);
    ylabel('normalized frequency'); 
    xlabel('log_{10} (variance explained)') 
    title('ER model noise time to expression','FontSize',12);   
   
 subplot(3,1,2);
    histin = vnI; histin(histin==0)=NaN; 
    hist(log10(abs(histin)));
    ylim([0,length(grad_mEn)]); 
   set(gca,'YtickLabel',str2mat('0','0.2','0.4','0.6','0.8','1'),'fontsize',F);
    ylabel('normalized frequency'); 
    xlabel('log_{10} (variance explained)');
    title('IR model noise time to expression','FontSize',12); 
    legend('K01','K10','K12','K21','K23','K24','K32','K35','K42','K45',...
     'K56','K65','K67','K78','Kab','Kba','Location','West');
  
 subplot(3,1,3);
     histin = vdn;  histin(histin==0)=NaN; 
   hist(log10(abs(histin)));      
   ylim([0,length(grad_mEn)]); 
   set(gca,'YtickLabel',str2mat('0','0.2','0.4','0.6','0.8','1'),'fontsize',F);
    ylabel('normalized frequency'); 
    xlabel('log_{10} (variance explained)') 
    title('difference in noise in time to expression','FontSize',12); 

  
set(gcf,'color','white');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%  
% Difference in Mean, Variance and Noise   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

load modelFull_diffs;  % 5000 data points

% %  distribution of model differences in mean 
am = sum(dmv>0)/N;
dmp = log10(dmv);
dmn = log10(-dmv(imag(dmp)>0));
dmp(imag(dmp)>0) = NaN;

% %   distribution of model differences in variance
as = sum(dsv>0)/N;
dsp = log10(dsv);
dsn = log10(-dsv(imag(dsp)>0));
dsp(imag(dsp)>0) = NaN;

% %  distribution of model differences in noise
an = sum(dnv>0)/N;
dnp = log10(dnv);
dnn = log10(-dnv(imag(dnp)>0));
dnp(imag(dnp)>0) = NaN;

% % distribution of model differences in noise
% an = sum(dn2>0)/N;
% dnp = log10(dn2);
% dnn = log10(-dn2(imag(dnp)>0));
% dnp(imag(dnp)>0) = NaN;

figure(1); clf;
subplot(3,1,1); 
hist(dmp,20);  
h = findobj(gca,'Type','patch'); hold on;  
set(h,'FaceColor','r'); % alpha(.2); 
hist(dmn,20,'b');
xlabel('log_{10}(\mu_I-\mu_E)','FontSize',16); 
legend(['\mu_I > \mu_E  ',num2str(100*am,3),' %'],...
    ['\mu_I < \mu_E  ',num2str(100-100*am,3),' %'],'Location','Best');
set(gca,'FontSize',16);
title('Mean \tau','FontSize',16);

subplot(3,1,2); 
hist(dsp,20); hold on; 
h = findobj(gca,'Type','patch'); 
set(h,'FaceColor','r');
hist(dsn,20);
xlabel('log_{10}(\sigma^2_I-\sigma^2_E)','FontSize',16);
legend(['\sigma^2_I > \sigma^2_E  ',  num2str(100*as,3), ' %'],...
    ['\sigma^2_I < \sigma^2_E  ', num2str(100-100*as,3), ' %'],'Location','Best');
set(gca,'FontSize',16);
title('Variance \tau','FontSize',16);

subplot(3,1,3); 
hist(dnp,20); hold on; 
h = findobj(gca,'Type','patch'); 
set(h,'FaceColor','r');
hist(dnn,20);
xlabel('log_{10}(n_I-n_E)','FontSize',16);
legend(['n_I > n_E  ',  num2str(100*an,3), ' %'],...
    ['n_I < n_E  ', num2str(100-100*an,3), ' %'],'Location','Best');
set(gca,'FontSize',16);
title('Noise \tau','FontSize',16);

  
set(gcf,'color','white');

% saveas(gcf, 'figure', 'png')









% 
% %------
%  
%  figure(2); clf; colormap(jet);
%   histin = grad_mEn; histin(histin==0)=NaN; 
%    subplot(3,1,1);
%    hist(log10(abs(histin)));     set(gca,'FontSize',6);
%     histin = grad_mIn; histin(histin==0)=NaN; 
%     title('ER model mean time to expression','FontSize',12);
%     subplot(3,1,2);
%    hist(log10(abs(histin)));
%    set(gca,'FontSize',6);
%       title('IR model mean time to expression','FontSize',12); 
%     legend('K01','K10','K12','K21','K23','K24','K32','K35','K42','K45',...
%      'K56','K65','K67','K78','Kab','Kba');
%    subplot(3,1,3);
%      histin = grad_dmn;  histin(histin==0)=NaN; 
%    hist(log10(abs(histin)));     set(gca,'FontSize',6);
%     title('difference in mean time to expression','FontSize',12); 
%    
% 
%  figure(5); clf; 
%   histin = grad_dmn; histin(histin==0)=NaN; 
%    subplot(3,1,1);
%    hist(log10(abs(histin)));     set(gca,'FontSize',6);
%     histin = grad_dsn; histin(histin==0)=NaN; 
%     subplot(3,1,2);
%    hist(log10(abs(histin)));
%    set(gca,'FontSize',6);
%     legend('K01','K10','K12','K21','K23','K24','K32','K35','K42','K45',...
%      'K56','K65','K67','K78','Kab','Kba','Location','Best');
%    subplot(3,1,3);
%      histin = grad_dnn;  histin(histin==0)=NaN; 
%    hist(log10(abs(histin)));     set(gca,'FontSize',6);
% 
% 
%  figure(5); clf; hist(log10(abs(grad_dmn(:,1:2))));
%  
% % Ks = 5; 
% % figure(2); clf; figure(3); clf;
% % C = colormap(cool(Ks));
% % for k=1:Ks
% %     [value,frequency,norm] = my_hist(log(abs(grad_mEn(:,k))),20);
% % figure(2); hold on; plot(value,frequency/norm,'Color',C(k,:),'linewidth',3);
% %  figure(3); hold on; bar(value, frequency./norm, 'FaceColor',C(k,:),'EdgeColor',C(k,:));
% %  end
% 
% % plot normalized histograms
%  
%  mean_mE = mean(log(abs(grad_mEn))); 
%  mean_sE = mean(log(abs(grad_sEn )));
%  mean_mI = mean(log(abs(grad_mIn )));
%  mean_sI = mean(log(abs(grad_sIn))); 
%  mean_nE = mean(log(abs(grad_nEn))) ;
%  mean_nI = mean(log(abs(grad_nIn))) ;
% 
%  mean_dm = mean(log(abs(grad_dmn))) ;
%  mean_ds = mean(log(abs(grad_dsn))) ;
%  mean_dn = mean(log(abs(grad_dnn))) ;
%      
%      figure(1); clf;  colormap(summer(3));
% subplot(3,1,1); 
%     bar(([mean_mI;mean_mE;mean_dm]'));
%      set(gca,'XtickLabel',...
%          str2mat('K01','K10','K12','K21','K23','K24','K32','K35','K42','K45','K56','K65','K67','K78','Kab','Kba'),...
%          'XTick',1:16,'fontsize',12);
%       ylabel('log(abs( gradient ))');
%      title('expression speed, mean sensitivity'); 
% subplot(3,1,2);
%     bar(([mean_sI;mean_sE;mean_ds]'));
%      set(gca,'XtickLabel',...
%      str2mat('K01','K10','K12','K21','K23','K24','K32','K35','K42','K45','K56','K65','K67','K78','Kab','Kba'),...
%         'XTick',1:16,'fontsize',12);
%       ylabel('log(abs( gradient ))');
%      title('expression variance, mean sensitivity');
% subplot(3,1,3); 
%     bar(([mean_nI;mean_nE;mean_dn]'));
%      set(gca,'XtickLabel',...
%      str2mat('K01','K10','K12','K21','K23','K24','K32','K35','K42','K45','K56','K65','K67','K78','Kab','Kba'),...
%          'XTick',1:16,'fontsize',12);
%       ylabel('log(abs( gradient ))');
%      title('expression noise power, mean sensitivity');
%     legend('Initiation Model','Elongation Model','Difference','Location','Best');