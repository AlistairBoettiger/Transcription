%   Numerically Invert Laplace Transforms for a chosen set of model
%   parameters for the Full ER and IR models 
%
%
% Alistair Boettiger                                   Date Begun: 10/18/09
% Levine Lab                                Functionally Complete: 10/19/09
%                                                   Last Modified: 11/01/09
%
% Modified 11/01/09
% Modified 11/19/10 : Use Darzaq's steps as single rate limiting.  


% load modelFull_solns % analysis 
%save modelFull_pdists vI vE; 
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

lambda = sym('lambda','real');
 
%              1        2       3       4        5       6       7       8        9         10       11       12       13        14     15     16
pars =   [  k12,      k21,   k23,      k24,     k32,    k35,    k53,    k54,     k42,       k45,     k56,     k65,    k67,       k78,   kab,   kba];         
% vals = 10*[6*.0216, 6*.145, 6*.0216, 6*.0216, 6*.145, 6*.0216 6*.145,  6*.145,  6*.0216,  6*.0216,  3*.00159, .001,  3*.00159,  3*.00159,.001, .1, ];

  % vals = 5*[.0216, .145, 2, 2, 2, 2, 2,  2,  2,  2,  .00159, .1,  2,  2,.01, 1 ];
  

  vals = [.00216, .145, 2, 2, 2, 2, 2,  2,  2,  2,  .159, .1,  2,  2,.1, 1 ];  % extreme case. 

% vals = rand(16,1);

% pars = [ kab,  kba,    k12,  k21,  k23,   k24,   k32,   k35,  k42,  k45,   k56,   k65    ,k67,  k78];
% vals = 3*[1E-3, 3E-1, 6*.0216, 6*.145, 6*.0216, 6*.0216, 6*.145, 6*.0216 6*.145, 6*.0216, 3*.00159, .001,  3*.00159,  3*.00159];

fI = subs(vI,pars,vals); 

fE = subs(vE,pars,vals); 

fIs = char(fI(1));
fIs = strrep(fIs,'*','.*');
fIs = strrep(fIs,'/','./');
fIs = strrep(fIs,'^','.^');
lambda = linspace(1E-8,4E3,1000); 

fIv = eval(fIs);  
fIl = fIv(fIv>0);

%     figure(2); clf; subplot(2,1,1); ezplot(fI(1),[1,50]); ylim([0,5E-10]);
%   subplot(2,1,2); ezplot(fE(1),[1,50]); ylim([0,5E-8]);
%    title(' ')
%    

% lE = subs(fE(1),lambda,2);


   %   lam = linspace(10,4E4,1000); 
      
        lam = linspace(10,4E5,1000);   % for extreme case
        
% Here we actually numerically invert the Laplace space solutions of the
% model.  
   FI = invlap2(fI(1), lam'); 
   FE = invlap2(fE(1), lam'); 


 
%  xlim([0,1E4/60]);

  dt = (max(lam)-min(lam)) / length(lam) ;
  
norm_E = sum(FE*dt)  % integrate p(t)*dt 
norm_I = sum(FI*dt) 
  
mean_E = lam*FE*dt/60  % integrate t*p(t)*dt 
mean_I = lam*FI*dt/60

std_E =  sqrt(lam.^2*FE*dt -  (lam*FE*dt)^2)/60
std_I = sqrt(lam.^2*FI*dt - (lam*FI*dt)^2)/60


  save compfull_pdist_simdata_X2; 

% load  compfull_pdist3_data;

 figure(1); clf; plot(lam/60,FI,'b','LineWidth',2); hold on; 
  plot(lam/60,FE,'r','LineWidth',2);    xlim([0,lam(end)/60]); %ylim([0,1.2E-3]);
title(['\mu_{ER} = ',num2str(mean_E,2), ' \sigma_{ER} = ', num2str(std_E,2),  ...
    '       \mu_{IR} = ',num2str(mean_I,2), ' \sigma_{IR} = ', num2str(std_I,2),...
    '  minutes'],'FontSize',15);
 xlabel('time (minutes)','FontSize',15); 
set(gca,'FontSize',15); set(gcf,'color','w');
legend('IR Model','ER Model')



  