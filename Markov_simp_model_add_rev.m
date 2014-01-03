%%                       Simple model as a control
%

% Alistair Boettiger                                      Date Begun: 04/09
% Levine Lab                                     Functional Since: 06/22/09
%                                                   Last Modified: 11/09/10
%% Notes
% This code uses the directly computable simple model as a control to test
% the implimentation of the Series Decomposition routine in SeriesDecomp.m.
% The model is presented in Figure 1 in the main text and presented again
% with a complete analysis in the Appendix. 
%
%% Modifications
% Modified 08/07/09 to add explicit representation of ``transcribed''
% state.  
% Modified 02/17/10 to change result plotting
% Modified 11/04/10 to make paused state reversibly exitable.  
% Modified 11/09/10 to plot reversible paused state.  
%
%% Required functions
% SeriesDecomp.m 
%
%%

clear all;
tic
% % % parameters
kab = sym('kab','real');
kba = sym('kba','real');
k12 = sym('k12','real'); 
k21 = sym('k21','real'); 
k23 = sym('k23','real'); 
k34 = sym('k34','real'); 
k32 = sym('k32','real'); 

% Simple Initiation Regulated Model
GI{1} = zeros(3);
GI{2} =[[-kab,kab,0,0,0];
        [kba,-k12-kba,k12,0,0];
        [0,k21,-k21-k23,k23,0];
        [0,0,k32,-k34-k32,k34];
        [0,0,0,0,0 ]];
   
   [m1I, m2I] = SeriesDecomp(GI); 

% Simple Elongation model
GE{1} = zeros(3);
%            1A 2A 3A  1B  2B  3B   4 
GE{2} = [[-k12-kab,k12,0,kab,0,0,0];       % 1A
        [k21,-k21-kab-k23,k23,0,kab,0,0];  % 2A
        [0,k32,-kab-k32,0,0,kab,0 ];             % 3A
        [kba,0,0,-kba-k12,k12,0,0];        % 1B
        [0,kba,0,k21,-kba-k21-k23,k23,0];  % 2B
        [0,0,kba,0,k32,-k32-k34-kba,k34] ;       % 3B
        [0,0,0,0,0,0,0]];                  % 4   
    
        
   [m1E, m2E] = SeriesDecomp(GE); 
   

   vars = [k12,k21,k23,k34,k32,kab,kba];
   vals = [1,1,1,1,1,1];

  % save model_simpIvE m1E m2E m1I m2I; 
   % load model_simpIvE

%    % Senstivity Analysis
%  grad_mE = jacobian(m1E,vars);
%  grad_sE = jacobian(m2E-m1E^2,vars);
%  grad_nE = jacobian((m2E-m1E^2)/m1E^2,vars);
%  grad_mI = jacobian(m1I,vars);
%  grad_sI = jacobian(m2I-m1I^2,vars);  
%  grad_nI = jacobian((m2I-m1I^2)/m1I^2,vars);
%  grad_dm = jacobian(m1E-m1I,vars);
%  grad_ds = jacobian(m2E-m1E^2-(m2I-m1I^2),vars);
%  grad_dn = jacobian(((m2E-m1E^2)/m1E^2-(m2I-m1I^2)/m1I^2),vars);
%  % clear m1E m2E m1I m2I; 

 
 vI = m2I - m1I.^2;
 vE = m2E -m1E.^2;
 nI = vI./m1I;
 nE = vE./m1E;
 
%save Markov_simp_data_addrev m1E m2E m1I m2I vI vE nI nE; 
%load Markov_simp_data_addrev 

%% Convert Solutions to text form
% This is a work-around for the 'subs' function, which is functional but
% extremely, phenomenally slow compared to 'eval'.  

%load Markov_simp_data_addrev;

 solns = {char(m1E); char(m2E); char(m1I); char(m2I); ...
    char(vE); char(vI); char(nE); char(nI)};
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
%  save Markov_simp_addrevB

%% Explore parameter space

% load  Markov_simp_addrev

N=10000; 
P = 7; 
M1E = zeros(N,1); M2E = M1E; M1I = M1E; M2I = M1E; 
dM = M1E; dV = M1E; dN = M1E;  
% GmE = zeros(N,P); GsE = GmE; GnE = GmE; GmI = GmE; GsI = GmE; GnI = GmE; 
vars = zeros(N,P); 

for i=1:N
 K12 = rand; K23 = rand; K34 = rand; K32 = rand;
 K21 = rand; Kab = rand; Kba = rand;
 vars(i,:) = [K12,K21,K23,K34,K32,Kab,Kba];

dM(i) = eval(solns{3}) / eval(solns{1});
dV(i) = eval(solns{6}) / eval(solns{5});
dN(i) = eval(solns{8}) / eval(solns{7});

% GmE(i,:) = eval(solns{5});
% GsE(i,:) = eval(solns{6});
% GnE(i,:) = eval(solns{7});
% GmI(i,:) = eval(solns{8});
% GsI(i,:) = eval(solns{9});
% GnI(i,:) = eval(solns{10});
end

 % save Markov_simp_dist_addrev;
%%
% load Markov_simp_dist_rev;  vmax = 4500; 
% load Markov_simp_dist_rev2; 
load Markov_simp_dist_addrev;  % add k32.  
% save plotdata_addk32

x = linspace(-4,18,150); vmax = 9000; 

figure(2); clf; set(gcf,'color','w');

subplot(3,1,1); hist(log(dM),x);
 hold on; plot([0,0],[0,vmax],'r--');
 xlabel('log(\mu_{IR}/\mu_{ER})'); ylim([0,vmax]);
 title('simple model: mean expression delay')

 subplot(3,1,2); hist(log(dV),x);
  hold on; plot([0,0],[0,vmax],'r--');
   xlabel('log(\sigma^2_{IR}/\sigma^2_{ER})');  ylim([0,vmax]);
   title('simple model: variance in delay')

subplot(3,1,3); hist(log(dN),x);
 hold on; plot([0,0],[0,vmax],'r--');
 xlabel('log(\eta_{IR}/\eta_{ER})');  ylim([0,vmax]);
 title('simple model: varibility in transcripts'); 
 
 IR_better = log(dM)<0;
 
 var_names = {'K12','K21','K23','K34','K32','Kab','Kba'};
 xx = linspace(0,1,8);
 figure(3); clf; set(gcf,'color','w'); 
 for k=1:length(var_names); 
    subplot(length(var_names),1,k); hist(vars(IR_better,k),xx); xlim([0,1]);
    title(var_names{k}); xlim([-.05,1.05]);
 end

 
 %%
 rmv = dM; rsv = dV;  rnv = dN;
 
 th = 0;
wm = .03;
m = log10(rmv); mr = m(m<=th); mb = m(m>th);
rb = round(th/wm); 
bb = round((max(m)-th)/wm);
mrb = round(-min(m)/wm); 
mbb = round(max(m)/wm);

ws = .03;
s = log10(rsv); sr = s(s<=th); sb = s(s>th);
srb = round(-min(s)/ws); 
sbb = round(max(s)/ws);

wn=.03;
n = log10(rnv); nr = n(n<=0); nb = n(n>0);
nrb = round(-min(n)/wn); 
nbb = round(max(n)/wn);
%%
C = [1,.4,.4];

figure(4); clf; subplot(3,1,1);
 hist(mr,mrb); hr = findobj(gca,'Type','patch');
 set(hr,'FaceColor',C,'EdgeColor',C); hold on;
  hist(mb,mbb); xlim([-2,2]);
xlabel('log(\mu_{IR}/\mu_{ER})'); set(gca,'YTickLabel',' ');
title('Mean Expression Speed');
mf = length(mr)/(length(mr)+length(mb)); 
legend([num2str(mf*100,2),'%  IR < ER']);

subplot(3,1,2);
 hist(sr,srb); 
 hs = findobj(gca,'Type','patch');
 set(hs,'FaceColor',C,'EdgeColor',C); hold on; 
 hist(sb,sbb); xlim([-2,2]);
xlabel('log(\sigma^2_{IR}/\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing'); 
sf = length(sr)/(length(sr)+length(sb)); 
legend([num2str(sf*100,2),'%  IR < ER']);

subplot(3,1,3);
hist(nr,nrb); hn = findobj(gca,'Type','patch');
set(hn,'FaceColor',C,'EdgeColor',C);
hold on;  hist(nb,nbb); 
% hb = findobj(gca,'Type','patch'); set(hb,'FaceColor','b','EdgeColor','b');
xlim([-2,2]);
xlabel('log(\eta_{IR}/\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number'); 
nf = length(nr)/(length(nr)+length(nb)); 
legend([num2str(nf*100,2),'%  IR < ER']);
set(gcf,'color','w');


toc

