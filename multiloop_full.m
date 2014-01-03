
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


 % save multiloop_data; 

 
 
 %% 
 load multiloop_data;
 


m1E =  m1E_58;
m2E = m2E_58;
m1I =  m1I_58;
m2I = m2I_58;
 
clear m1E_58 m2E_58 m1I_58 m2I_58 v vI vE

% Define variables
%       1    2    3  4    5  6   7  8    9   10  11  12 13  14  15  16
vars = [k12,k21,k23,k24,k32,k35,k53,k54,k42,k45,k56,k65,k67,k78,kab,kba];


% Print solutions in memory using strfind (much faster than 'subs' command.
solns = {char(m1E); char(m2E); char(m1I); char(m2I)};% ; ...
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
 % save loop58_data;


%% Explore Parameter Space
% substitutes in N random draws.
%  

tic
load loop58_data;
solns58 = solns; % rename to avoid overwriting

load complex_data-10-28-10 vars; % load vars used in tau_18 calculation


N=length(vars); 
P = 16; 
 M1E = zeros(N,1); M2E = M1E; M1I = M1E; M2I = M1E; 

% These are momemnts now not standard deviations. 

% use same vars 

for i=1:N ; 
 K12 = vars(i,1);  K21 = vars(i,2); K23 = vars(i,3); K24 = vars(i,4); K32 = vars(i,5); K35 = vars(i,6);
 K53=vars(i,7); K54=vars(i,8); K42 = vars(i,9);  K45 = vars(i,10);
 K56 = vars(i,11); K65 = vars(i,12);  K67 = vars(i,13); K78 = vars(i,14);  Kab = vars(i,15); Kba = vars(i,16);
% vars(i,:) = [K12,K21,K23,K24,K32,K35,K53,K54,K42,K45,K56,K65,K67,K78,Kab,Kba];
M1E(i) = eval(solns58{1});  % ER model, 1st moment i.e. mean
M2E(i) = eval(solns58{2});   % ER model, 2nd moment
M1I(i) = eval(solns58{3}); % IR model, 1st moment
M2I(i) = eval(solns58{4});  % IR model, 2nd moment
%M2E(i)
end



% % Get rid of unused cells
% M1E = nonzeros(M1E); 
% M2E = nonzeros(S2E);
% M1I = nonzeros(M1I);
% M2I = nonzeros(S2I);

S2E = M2E - M1E.^2;
S2I = M2I - M1I.^2; 

nE = S2E./M1E.^2; % 'noise'/CoV of transition time, ER model
nI = S2I./M1I.^2; % 'noise'/CoV of transition time, IR model


dmv = M1I - M1E; % difference in means
dsv = S2I - S2E; % difference in variances
dnv = nI - nE; % difference in CoV



rmv = M1I./M1E; % ratio of means
rsv = S2I./S2E; % ratio of variances
rtv = nI ./ nE; % ratio of CoVs 
toc

save multiloop_pdata; 
%%
load multiloop_pdata;
 

M1E_58 = M1E;
M2E_58 = M2E;

M1I_58 = M1I;
M2I_58 = M2I;

load complex_data-10-28-10

M2E_18 = s2E + M1E.^2; 
M2I_18 = s2I + M1I.^2;
M1E_18 = M1E;
M1I_18 = M1I; 

save multiloop_pre_b_data; 
%%
load multiloop_pre_b_data; 
xmin = -1; xmax = 2; bins = 60;

% the correct way to do this, with function calls.  

hout = figure(2); clf; 
subplot(1,4,1); b = 1; 
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins)

figure(2); subplot(1,4,2); b = .9;
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins)

figure(2); subplot(1,4,3);  b = .3;
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins)

figure(2); subplot(1,4,4);  b = 0;
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins)

set(gcf,'color','w');
% saveas(hout,[fout,'/scaffold_effects.ai']); 


%%

S2I_58 = M2I_58 - M1I_58.^2;
S2E_58 = M2E_58 - M1E_58.^2;

NI_58 = S2I_58./M2I_58;
NE_58 = S2E_58./M2E_58;

rmv2 = M1I_58./M1E_58;
rsv2 = S2I_58./S2E_58;
rtv2 = NI_58./NE_58 ;

bins = 200; 
m = log10(rmv2); s = log10(rsv2); n = log10(rtv2);
mb = linspace(xmin,xmax-.025,bins); % max(m)/xmax*bins;
sb =linspace(xmin,xmax-.025,bins); % max(s)/xmax*bins;
nb =linspace(xmin,xmax-.025,bins); % max(n)/xmax*bins;
fano = log10(r_fano); nf = linspace(-1,xmax-.025,bins); 

C = [1,.4,.4];

F = 10; 

figure(3); clf; subplot(3,1,1); hist(m,mb);
  mh = hist(m,mb); xlim([xmin,xmax]);
  hold on; plot([0,0],[0,4/3*max(mh)],'m--','LineWidth',2);
  ylim([0,4/3*max(mh)]);
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
 xlabel('log_{10}(\mu_{IR}) - log_{10}(\mu_{ER})'); set(gca,'YTickLabel',' ');
 ylabel('frequency');  set(gca,'FontSize',F);
title('Mean Expression Speed, \mu','FontSize',F);

subplot(3,1,2);
 hist(s,sb); xlim([xmin,xmax]);
   sh = hist(s,sb); 
  hold on; plot([0,0],[0,4/3*max(sh)],'m--','LineWidth',2);
    h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
  ylabel('frequency'); 
  ylim([0,4/3*max(sh)]);
xlabel('log_{10}(\sigma^2_{IR}) - log_{10}(\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing, \sigma^2','FontSize',F); set(gca,'FontSize',F);


subplot(3,1,3);
 hist(n,nb); xlim([xmin,xmax]);
    nh = hist(n,nb); 
  hold on; plot([0,0],[0,4/3*max(nh)],'m--','LineWidth',2);
  
    h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
  ylim([0,4/3*max(nh)]); ylabel('frequency'); 
xlabel('log_{10}(\eta_{IR}) -log_{10}(\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number, \eta','FontSize',F); 
set(gcf,'color','w'); set(gca,'FontSize',F);

