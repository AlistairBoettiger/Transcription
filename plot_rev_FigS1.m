
load save Markov_simp_dist_reg3;

xmin = -4;
xmax = 4; 
x = linspace(xmin,xmax,50); vmax = 1500; 


figure(2); clf; set(gcf,'color','w');

subplot(3,1,1); hist(log(dM),x);
 hold on; plot([0,0],[0,vmax],'r--');
 xlabel('log(\mu_{IR}/\mu_{ER})'); ylim([0,vmax]); xlim([xmin,xmax]);
 title('simple model: mean expression delay')

 subplot(3,1,2); hist(log(dV),x);
  hold on; plot([0,0],[0,vmax],'r--');
   xlabel('log(\sigma^2_{IR}/\sigma^2_{ER})');  ylim([0,vmax]); xlim([xmin,xmax]);
   title('simple model: variance in delay')

subplot(3,1,3); hist(log(dN),x);
 hold on; plot([0,0],[0,vmax],'r--');
 xlabel('log(\eta_{IR}/\eta_{ER})');  ylim([0,vmax]); xlim([xmin,xmax]);
 title('simple model: varibility in transcripts'); 
 
 IR_better = log(dM)<0;
 
 var_names = {'K12','K21','K23','K34','Kab','Kba'};
 xx = linspace(0,1,20);
 figure(3); clf; set(gcf,'color','w'); 
 for k=1:6
    subplot(6,1,k); hist(vars(IR_better,k),xx); xlim([0,1]); xlim([xmin,xmax]);
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
  hist(mb,mbb); xlim([-.2,2]);
xlabel('log(\mu_{IR}/\mu_{ER})'); set(gca,'YTickLabel',' ');
title('Mean Expression Speed');
mf = length(mr)/(length(mr)+length(mb)); 
legend([num2str(mf*100,2),'%  IR < ER']);

subplot(3,1,2);
 hist(sr,srb); 
 hs = findobj(gca,'Type','patch');
 set(hs,'FaceColor',C,'EdgeColor',C); hold on; 
 hist(sb,sbb); xlim([-.2,2]);
xlabel('log(\sigma^2_{IR}/\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing'); 
sf = length(sr)/(length(sr)+length(sb)); 
legend([num2str(sf*100,2),'%  IR < ER']);

subplot(3,1,3);
hist(nr,nrb); hn = findobj(gca,'Type','patch');
set(hn,'FaceColor',C,'EdgeColor',C);
hold on;  hist(nb,nbb); 
% hb = findobj(gca,'Type','patch'); set(hb,'FaceColor','b','EdgeColor','b');
xlim([-.2,2]);
xlabel('log(\eta_{IR}/\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number'); 
nf = length(nr)/(length(nr)+length(nb)); 
legend([num2str(nf*100,2),'%  IR < ER']);
set(gcf,'color','w');



 