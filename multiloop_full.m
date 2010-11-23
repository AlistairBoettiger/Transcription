
%%                                   multiloop_full.m

% Alistair Boettiger                                   Date Begun: 10/12/10
% Levine Lab                                        Last Modified: 10/12/10

% Description: 
%
% Modified 11/22/10 to use pdist3 to make computations
%          and to fix error in definition of p that propigated to this code. 


clear all;
 load compfull_pdist3_data; 

 lambda = sym('lambda','real'); 

 %% Evaluate starting at beginning of third pinch-point.  
 
 % IR regulated
 v = vI(3);
 temp = -diff(v,lambda);
m1I_58 = subs(temp,lambda,0);
m2I_58  = -subs(diff(temp,lambda),lambda,0);


% ER regulated
 v = vE(3);
 temp = -diff(v,lambda);
m1Et = subs(temp,lambda,0);
m2Et = -subs(diff(temp,lambda),lambda,0);

% Rapid equilbirium of enhancer model
syms lambda;
f = length(Gen);
Gen_star = Gen; 
Gen_star(f,:) = zeros(1,f); 
iGen = inv(lambda*eye(f) - Gen_star);
mS = -subs(diff(lambda*iGen(1,f),lambda),lambda,0); % mean for enhancer chain 
m2S = subs(diff(diff(lambda*iGen(1,f),lambda),lambda),lambda,0); % 2nd moment for enhancer chain 
p = kab/(kab+kba); % probability of find ing the system in the final state


m1E_58  = m1Et + mS*(1-p);  % First moment of ER model
m2E_58 = m2Et + (1-p)*(2*m1Et*mS+m2S); % Second moment of ER model

 save multiloop_data; 

 
 
 %% 
 load multiloop_data;
 
load complex_data-10-28-10;

 
 b = .5; 
 % b probability of looping back to start. 
 
 m1bE = b*m1E + (1-b)*m1E_58;
 
 m1bI = b*m1I + (1-b)*m1I_58;
 
s2bE =  m2E_58^2 - (b^2*m1E + (1-b)^2*m1E_58^2 + 2*b*(1-b)*m1E*m1E_58^2);
 
s2bI =  m2I_58^2 - (b^2*m1I + (1-b)^2*m1I_58^2 + 2*b*(1-b)*m1I*m1I_58^2);

nbE = s2bE/m1bE;

nbI = s2bI/m1bI; 
 
  save multiloop_data_bp5; 
%%

load multiloop_data_bp5;


m1E =  m1bE;
s2E = s2bE;
m1I =  m1bI;
s2I = s2bI;
 
clear m1E_58 m2E_58 m1I_58 m2I_58 v vI vE

% Define variables
%       1    2    3  4    5  6   7  8    9   10  11  12 13  14  15  16
vars = [k12,k21,k23,k24,k32,k35,k53,k54,k42,k45,k56,k65,k67,k78,kab,kba];
% vals = [ 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 ,.1 , 1 ];

% Compute gradients
% grad_nE = jacobian((m2E-m1E^2)/m1E,vars);
%  grad_nI = jacobian((m2I-m1I^2)/m1I,vars);
 
 
% % Print solutions to Text File (fast method)
 
%  fid = fopen('modelFull_solns.txt', 'wt');
% fprintf(fid, ['m1E = ', char(m1E), '\n', '  \n',...
%     'm2E = ', char(m2E), '\n', '  \n',...
%     'm1I = ', char(m1I), '\n', '  \n',...
%     'm2I = ', char(m2I), '\n', '  \n',...
%     'grad_mE = ', char(grad_mE), '\n', '  \n',...
%     'grad_mI = ', char(grad_mI), '\n', '  \n',...    
%     ]  );
% fclose(fid);

% convert solutions in symbolic expression to matlab commands with
% variables for inputs.  This is much faster than 'subs' command, sadly. 

% Print solutions in memory using strfind (much faster than 'subs' command.
solns = {char(m1E); char(s2E); char(m1I); char(s2I)};% ; ...
    % char(grad_nE); char(grad_nI)};
for k =1:length(solns) 
    solns{k}(strfind(solns{k},'k'))='K';
    % remove text 'matrix' from beginning of character matrices
    try % this will only work if the statement begins with 'matrix(
    solns{k}(strfind(solns{k},'matrix'):6)='      ';
    catch
    end
%     if solns{k}(end)==')' % remove terminal parenthesis 
%          solns{k}(end)=' ';
%      end
end


% % If you don't want to run this cell again, save it's data
 save multiloop_solns_bp5;


%% Explore Parameter Space
% substitutes in N random draws.
%  

load multiloop_solns_bp5;

N=1000; 
P = 16; 
M1E = zeros(N,1); S2E = M1E; M1I = M1E; S2I = M1E; 
vars = zeros(N,P); 



for i=1:N
 K12 = rand; K23 = rand; K24 = rand;  K35 = rand;  K45 = rand; K53=rand; K54=rand;
 K56 = rand; K67 = rand; K78 = rand;     Kab = rand; Kba = rand;
 K21 = rand; K32 = rand; K42 = rand; K65 = rand;
 vars(i,:) = [K12,K21,K23,K24,K32,K35,K53,K54,K42,K45,K56,K65,K67,K78,Kab,Kba];
M1E(i) = eval(solns{1});  % ER model, 1st moment i.e. mean
S2E(i) = eval(solns{2});   % ER model, 2nd moment
M1I(i) = eval(solns{3}); % IR model, 1st moment
S2I(i) = eval(solns{4});  % IR model, 2nd moment
end



% Get rid of unused cells
M1E = nonzeros(M1E); 
S2E = nonzeros(S2E);
M1I = nonzeros(M1I);
S2I = nonzeros(S2I);


nE = S2E./M1E.^2; % 'noise'/CoV of transition time, ER model
nI = S2I./M1I.^2; % 'noise'/CoV of transition time, IR model


dmv = M1I - M1E; % difference in means
dsv = S2I - S2E; % difference in variances
dnv = nI - nE; % difference in CoV



rmv = M1I./M1E; % ratio of means
rsv = S2I./S2E; % ratio of variances
rtv = nI ./ nE; % ratio of CoVs 



 save multiloop_pdata_bp5;

%%
load multiloop_pdata_bp5;
 
 slowR = min(vars(:,1:14)');
 rtv(rtv<0) = 1E-10; 
 
 kab_slow = vars(:,15) < .1*slowR';
 kba_slow = vars(:,16) < .1*slowR';
 noneq = kab_slow | kba_slow;
 
 rmv2 = rmv(~noneq);
 rsv2 = rsv(~noneq); 
 rtv2 = rtv(~noneq); 

 
 
offset = .025;
bins = 70;
xmax =2 + offset;
m = log10(rmv2); s = log10(rsv2); n = log10(rtv2);
mb = linspace(-2,xmax-.025,bins); % max(m)/xmax*bins;
sb =linspace(-2,xmax-.025,bins); % max(s)/xmax*bins;
nb =linspace(-2,xmax-.025,bins); % max(n)/xmax*bins;

C = [1,.4,.4];

F = 10; 

figure(3); clf; subplot(2,2,1); hist(m,mb);
  mh = hist(m,mb); xlim([-2,xmax]);
  hold on; plot([0,0],[0,4/3*max(mh)],'m--','LineWidth',2);
  ylim([0,4/3*max(mh)]);
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
 xlabel('log_{10}(\mu_{IR}) - log_{10}(\mu_{ER})'); set(gca,'YTickLabel',' ');
 ylabel('frequency');  set(gca,'FontSize',F);
title('Mean Expression Speed, \mu','FontSize',F);

subplot(2,2,2);
 hist(s,sb); xlim([-2,xmax]);
   sh = hist(s,sb); 
  hold on; plot([0,0],[0,4/3*max(sh)],'m--','LineWidth',2);
    h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
  ylabel('frequency'); 
  ylim([0,4/3*max(sh)]);
xlabel('log_{10}(\sigma^2_{IR}) - log_{10}(\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing, \sigma^2','FontSize',F); set(gca,'FontSize',F);


subplot(2,2,3);
 hist(n,nb); xlim([-2,xmax]);
    nh = hist(n,nb); 
  hold on; plot([0,0],[0,4/3*max(nh)],'m--','LineWidth',2);
  
    h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
  ylim([0,4/3*max(nh)]); ylabel('frequency'); 
xlabel('log_{10}(\eta_{IR}) -log_{10}(\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number, \eta','FontSize',F); 
set(gcf,'color','w'); set(gca,'FontSize',F);


%  load compfull_pdist3_data; 
% 
%  
% subplot(2,2,4); 
%  plot(lam/60,FI,'b','LineWidth',2); hold on; 
%   plot(lam/60,FE,'r','LineWidth',2);    xlim([0,lam(end)/60]); %ylim([0,1.2E-3]);
% % title(['\mu_{ER} = ',num2str(mean_E,2), ' \sigma_{ER} = ', num2str(std_E,2),  ...
% %     '       \mu_{IR} = ',num2str(mean_I,2), ' \sigma_{IR} = ', num2str(std_I,2),...
% %     '  minutes'],'FontSize',15);
%  xlabel('time (minutes)','FontSize',F); 
% set(gca,'FontSize',F); set(gcf,'color','w');
% legend('IR Model','ER Model')
% title('PDF for \tau'); 
% 
