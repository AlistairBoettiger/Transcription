
%                               model_CompFull2b.m
%
% Alistair Boettiger                                   Date Begun: 04/06/09
% Levine Lab                                     Functional Since: 05/05/09
% Functional computational code                     Last Modified: 11/07/09

% Notes:
% compute moment generating function of first passage times for arbitrary
% pinch point decomposed (series) Markov Chains using Matlab's symbolic
% Algebra routines

% Modified 10/18/09 to remove state 0. 
% Modified 10/20/09 to clean up code
% Modified 11/07/09 to fix bug k53 and k54 missing (had used k35 and k45)


clear all;
% 
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


Gen =  [[-kab,kab];  % enhancer chain
       [kba,-kba]];


% Submatrices for initiation regulated system
GI{1} = zeros(3); 

GI{2}=[[-kab,kab,0];
    [kba, -kba-k12, k12];
    [0, k21, -k21]];

GI{3}=[[-k23-k24, k23, k24, 0];
    [k32, -k32-k35, 0, k35];
    [k42, 0, -k42-k45, k45];
    [0, k53, k54, -k53-k54]];

GI{4} =[[-k56, k56, 0, 0];
    [k65, -k65-k67, k67, 0];
    [0, 0, -k78, k78];
    [0, 0, 0, 0]];

 [m1I,m2I] = SeriesDecomp(GI); 

% Submatrices for elongation regulated system 
GE{1} = zeros(3); 


GE{2}=[[-k12, k12];
       [k21, -k21]];


GE{3}=[[-k23-k24, k23, k24, 0];
    [k32, -k32-k35, 0, k35];
    [k42, 0, -k42-k45, k45];
    [0, k53, k54, -k53-k54]];

GE{4}=[[-k56, k56, 0, 0];
    [k65, -k65-k67, k67, 0];
    [0, 0, -k78, k78];
    [0,0,0,0]]; % kab kba integrated below

 [m1Et,m2Et] = SeriesDecomp(GE); % means without enhnacer yet
%    

syms lambda;
f = length(Gen);
Gen_star = Gen; 
Gen_star(f,:) = zeros(1,f); 
iGen = inv(lambda*eye(f) - Gen_star);
mS = -subs(diff(lambda*iGen(1,f),lambda),lambda,0); % mean for enhancer chain 
m2S = subs(diff(diff(lambda*iGen(1,f),lambda),lambda),lambda,0); % 2nd moment for enhancer chain 

p = kba/(kab+kba); % probability of find ing the system in the final state
% should be computed directly from the final element of the kernel of the
% matrix Gen.  Matlab refuses to compute kernel through left divide.  
m1E = m1Et + mS*(1-p);
m2E = m2Et + (1-p)*(2*m1Et*mS+m2S); 


%       1    2    3  4    5  6   7  8    9   10  11  12 13  14  15  16
vars = [k12,k21,k23,k24,k32,k35,k53,k54,k42,k45,k56,k65,k67,k78,kab,kba];
% vals = [ 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 ,.1 , 1 ];

 grad_mE = jacobian(m1E,vars);
 grad_sE = jacobian(m2E-m1E^2,vars);
 grad_mI = jacobian(m1I,vars);
 grad_sI = jacobian(m2I-m1I^2,vars);
 grad_nE = jacobian((m2E-m1E^2)/m1E,vars);
 grad_nI = jacobian((m2I-m1I^2)/m1I,vars);
 
 
% % Print solutions 
 
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

solns = {char(m1E); char(m2E); char(m1I); char(m2I); ...
    char(grad_mE); char(grad_sE); char(grad_nE); char(grad_mI); char(grad_sI); char(grad_nI)};
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


 save modelFull2b_solns;




%%
 load modelFull2b_solns;


N=10000; 
P = length(vars); 
M1E = zeros(N,1); M2E = M1E; M1I = M1E; M2I = M1E; 
GmE = zeros(N,P); GsE = GmE; GnE = GmE; GmI = GmE; GsI = GmE; GnI = GmE; 
vars = zeros(N,P); 

for i=1:N
 K12 = rand; K23 = rand; K24 = rand;  K35 = rand;  K45 = rand; K53=rand; K54=rand;
 K56 = rand; K67 = rand; K78 = rand;     Kab = rand; Kba = rand;
 K21 = rand; K32 = rand; K42 = rand; K65 = rand;
 vars(i,:) = [K12,K21,K23,K24,K32,K35,K53,K54,K42,K45,K56,K65,K67,K78,Kab,Kba];
M1E(i) = eval(solns{1});
M2E(i) = eval(solns{2});
M1I(i) = eval(solns{3});
M2I(i) = eval(solns{4});
GmE(i,:) = eval(solns{5});
GsE(i,:) = eval(solns{6});
GnE(i,:) = eval(solns{7});
GmI(i,:) = eval(solns{8});
GsI(i,:) = eval(solns{9});
GnI(i,:) = eval(solns{10});
end



% save modelFull2b_data2; % 11/03/09


M1E = nonzeros(M1E);
M2E = nonzeros(M2E);
M1I = nonzeros(M1I);
M2I = nonzeros(M2I);



sE = M2E - M1E.^2;
sI = M2I - M1I.^2;
nE = sE./M1E.^2;
nI = sI./M1I.^2;
nE3 = sE./M1E.^3;
nI3 = sI./M1I.^3;

dmv = M1I - M1E;
dsv = sI - sE;
dnv = nI - nE;
%dn2 = nI2 - nE2;



rmv = M1I./M1E;
rsv = sI./sE;
%rnv = nI./nE;
rtv = (sI.^2./M1I) ./ (sE.^2./M1E); 


% rtv = sort(rtv);

% Export Data


Data = [GmE;GsE;GnE;GmI;GsI;GnI;vars];

fn = '/Users/alistair/Documents/Berkeley/Levine Lab/Data2b.txt';

dlmwrite(fn,Data); 
 %  modelFull2b_data;   11/10/09
 save modelFull2b_data2;  % 11/23/09




%%
load modelFull2b_data2; 

bins = 70;
xmax =1;
m = log10(rmv); s = log10(rsv); n = log10(rtv);
mb = linspace(-.2,xmax,bins); % max(m)/xmax*bins;
sb =linspace(-.2,xmax,bins); % max(s)/xmax*bins;
nb =linspace(-.2,xmax,bins); % max(n)/xmax*bins;

C = [1,.4,.4];

figure(3); clf; subplot(3,1,1); hist(m,mb);
  mh = hist(m,mb); xlim([-.2,xmax]);
  hold on; plot([0,0],[0,4/3*max(mh)],'m--','LineWidth',2);
  ylim([0,4/3*max(mh)]);
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
 xlabel('log(\mu_{IR}) - log(\mu_{ER})'); set(gca,'YTickLabel',' ');
title('Mean Expression Speed, \mu');

subplot(3,1,2);
 hist(s,sb); xlim([-.2,xmax]);
   sh = hist(s,sb); 
  hold on; plot([0,0],[0,4/3*max(sh)],'m--','LineWidth',2);
    h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
  
  ylim([0,4/3*max(sh)]);
xlabel('log(\sigma^2_{IR}) - log(\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing, \sigma^2'); 


subplot(3,1,3);
 hist(n,nb); xlim([-.2,xmax]);
    nh = hist(n,nb); 
  hold on; plot([0,0],[0,4/3*max(nh)],'m--','LineWidth',2);
  
    h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
  ylim([0,4/3*max(nh)]);
xlabel('log(\eta_{IR}) -log(\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number, \eta'); 
set(gcf,'color','w');


%%
load modelFull2b_data2; 
th = .1;
wm = .005;
m = log10(rmv); mr = m(m<=th); mb = m(m>th);
rb = round(th/wm); 
bb = round((max(m)-th)/wm);
ws = .005;
s = log10(rsv); sr = s(s<=th); sb = s(s>th);
srb = round(th/ws); 
sbb = round((max(s)-th)/ws);
wn=.003;
n = log10(rtv); nr = n(n<=0); nb = n(n>0);
nrb = round(-min(n)/wn); 
nbb = round(max(n)/wn);

C = [1,.4,.4];

figure(2); clf; subplot(3,1,1);
 hist(mr,rb); hr = findobj(gca,'Type','patch');
 set(hr,'FaceColor',C,'EdgeColor',C); hold on;
  hist(mb,bb); xlim([-.2,2]);
xlabel('log(\mu_{IR}/\mu_{ER})'); set(gca,'YTickLabel',' ');
title('Mean Expression Speed');
mf = length(mr)/(length(mr)+length(mb)); 
legend([num2str(mf*100,2),'% similar']);

subplot(3,1,2);
 hist(sr,srb); 
 hs = findobj(gca,'Type','patch');
 set(hs,'FaceColor',C,'EdgeColor',C); hold on; 
 hist(sb,sbb); xlim([-.2,2]);
xlabel('log(\sigma^2_{IR}/\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing'); 
sf = length(sr)/(length(sr)+length(sb)); 
legend([num2str(sf*100,2),'% similar']);

subplot(3,1,3);
hist(nr,nrb); hn = findobj(gca,'Type','patch');
set(hn,'FaceColor',C,'EdgeColor',C);
hold on;  hist(nb,nbb); 
% hb = findobj(gca,'Type','patch'); set(hb,'FaceColor','b','EdgeColor','b');
xlim([-.2,2]);
xlabel('log(\eta_{IR}/\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number'); 
nf = length(nr)/(length(nr)+length(nb)); 
legend([num2str(nf*100,2),'% IR < ER']);
set(gcf,'color','w');


%% 

i = 0;
figure(3); clf; figure(4); clf; figure(5); clf;
for j=1:16

k_e_fast = vars(m>th,j);
k_e_slow = vars(m<th,j);
k_e_sync = vars(s>th,j);
k_e_async = vars(s<th,j);
k_e_unif = vars(n>0,j);
k_e_nonunif = vars(n<0,j);

    i=i+1;
     figure(3);   subplot(16,2,i); hist(k_e_slow);
     hr = findobj(gca,'Type','patch'); set(hr,'FaceColor',C,'EdgeColor',C);
     figure(4); subplot(16,2,i); hist(k_e_async);
     hr = findobj(gca,'Type','patch'); set(hr,'FaceColor',C,'EdgeColor',C);
    figure(5);  subplot(16,2,i); hist(k_e_nonunif);
    hr = findobj(gca,'Type','patch'); set(hr,'FaceColor',C,'EdgeColor',C);

    i = i+1;
    figure(3); subplot(16,2,i); hist(k_e_fast); 
     hr = findobj(gca,'Type','patch'); set(hr,'FaceColor','k','EdgeColor','k');
    figure(4); subplot(16,2,i); hist(k_e_sync);
     hr = findobj(gca,'Type','patch'); set(hr,'FaceColor','k','EdgeColor','k');
    figure(5);  subplot(16,2,i); hist(k_e_unif);
     hr = findobj(gca,'Type','patch'); set(hr,'FaceColor','k','EdgeColor','k');
end

figure(3); set(gcf,'color','w'); 

figure(4); set(gcf,'color','w'); 

figure(5); set(gcf,'color','w'); 

% save model_results2b


%% Export Data


