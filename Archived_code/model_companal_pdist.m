%   Numerically Invert Laplace Transforms for a chosen set of model
%   parameters for the Full ER and IR models 
%
%
% Alistair Boettiger                                   Date Begun: 10/18/09
% Levine Lab                                Functionally Complete: 10/19/09
%                                                   Last Modified: 10/19/09

% load modelFull_solns % analysis 
%save modelFull_pdists vI vE; 
clear all;
load modelFull_pdists; 

% Specific Model Constants
kab = sym('kab','real');
kba = sym('kba','real');
k01 = sym('k01','real'); 
k10 = sym('k10','real'); 
k12 = sym('k12','real'); 
k21 = sym('k21','real'); 
k23 = sym('k23','real'); 
k24 = sym('k24','real'); 
k32 = sym('k32','real'); 
k35 = sym('k35','real'); 
k42 = sym('k42','real'); 
k45 = sym('k45','real'); 
k56 = sym('k56','real'); 
k65 = sym('k65','real'); 
k67 = sym('k67','real'); 
k78 = sym('k78','real'); 
kF = sym('kF','real');

pars = [ kab,  kba,  k01,   k10,   k12,  k21,  k23,   k24,   k32,   k35,  k42,  k45,   k56,   k65    ,k67,    k78];
%vals = [1E-2, 1E-2, .0216, .145, .0216, .145, .0216, .0216, .145, .0216 .145, .0216, .00159, .01,  .00159,  .00159];

%vals = [1E-2, 1E-2, 6*.0216, 6*.145, 6*.0216, 6*.145, 6*.0216, 6*.0216, 6*.145, 6*.0216 6*.145, 6*.0216, 3*.00159, .001,  3*.00159,  3*.00159];

vals = 3*[1E-2, 7E-2, 6*.0216, 6*.145, 6*.0216, 6*.145, 6*.0216, 6*.0216, 6*.145, 6*.0216 6*.145, 6*.0216, 3*.00159, .001,  3*.00159,  3*.00159];

fI = subs(vI,pars,vals); 
fE = subs(vE,pars,vals); 



%        lam = logspace(-2,8,1000);
%        lam = logspace(-2,6,1000);
     lam = linspace(1E-5,4E3,1000); 
   FI = invlap2(fI(1), lam'); 
FE = invlap2(fE(1), lam'); 

  figure(1); clf; plot(lam/60,FI,'b.'); hold on; 
  plot(lam/60,FE,'r.');    xlim([0,lam(end)/60]); %ylim([0,1.2E-3]);
%  xlim([0,1E4/60]);
  ylabel('probability','FontSize',15); xlabel('time (minutes)','FontSize',15); 

  dt = (max(lam)-min(lam)) / length(lam) ;

  norm_E = sum(FE*dt)  % integrate p(t)*dt 
  norm_I = sum(FI*dt)
  
  mean_E = lam*FE*dt/60  % integrate t*p(t)*dt 
mean_I = lam*FI*dt/60

std_E =  sqrt(lam.^2*FE*dt -  (lam*FE*dt)^2)/60
std_I = sqrt(lam.^2*FI*dt - (lam*FI*dt)^2)/60

figure(1);
title(['\mu_{ER} = ',num2str(mean_E,2), ' \sigma_{ER} = ', num2str(std_E,2),  ...
    '       \mu_{IR} = ',num2str(mean_I,2), ' \sigma_{IR} = ', num2str(std_I,2),...
    '  minutes'],'FontSize',15);
set(gca,'FontSize',15); set(gcf,'color','w');
legend('IR Model','ER Model')


% %    FI = invlap2(sum(fI), lam'); 
% % FE = invlap2(sum(fE), lam'); 
% 
% % optimized for vals 2
% 
%   figure(1); clf; plot(lam,FI,'b.'); hold on; 
%   plot(lam,FE,'r.');    xlim([0,1E4]); ylim([0,1.2E-3]);
%   ylabel('probability'); xlabel('time (seconds)'); 
% 
% % optimized for vals 1
%   figure(1); clf; plot(lam,FI,'b.'); hold on; 
%   plot(lam,FE,'r.');    xlim([0,1E5]); ylim([0,1.2E-4]);
%   ylabel('probability'); xlabel('time (seconds)'); 
%   
%   
%   figure(2); subplot(2,1,1); 
%    plot(lam,FE,'r.');
%   xlim([0,50000]); ylim([0,2E-4]);
%   
%   subplot(2,1,2); 
%    plot(lam,FI,'b.');
%   xlim([0,500000]); ylim([0,2E-5]);
%   

  