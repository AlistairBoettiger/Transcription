



load complex_data-10-28-10;

%% Figure 3 of Main Text
% load modelFull2b_data2; 
 % load model_data-10-06-10; % 10/06/10 

offset = .025;
bins = 70;
xmax =1 + offset;
m = log10(rmv); s = log10(rsv); n = log10(rtv);
mb = linspace(-.2,xmax-.025,bins); % max(m)/xmax*bins;
sb =linspace(-.2,xmax-.025,bins); % max(s)/xmax*bins;
nb =linspace(-.2,xmax-.025,bins); % max(n)/xmax*bins;
fano = log10(r_fano); nf = linspace(-1,xmax-.025,bins); 

C = [1,.4,.4];

F = 10; 

figure(3); clf; subplot(2,2,1); hist(m,mb);
  mh = hist(m,mb); xlim([-.2,xmax]);
  hold on; plot([0,0],[0,4/3*max(mh)],'m--','LineWidth',2);
  ylim([0,4/3*max(mh)]);
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
 xlabel('log_{10}(\mu_{IR}) - log_{10}(\mu_{ER})'); set(gca,'YTickLabel',' ');
 ylabel('frequency');  set(gca,'FontSize',F);
title('Mean Expression Speed, \mu','FontSize',F);

subplot(2,2,2);
 hist(s,sb); xlim([-.2,xmax]);
   sh = hist(s,sb); 
  hold on; plot([0,0],[0,4/3*max(sh)],'m--','LineWidth',2);
    h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
  ylabel('frequency'); 
  ylim([0,4/3*max(sh)]);
xlabel('log_{10}(\sigma^2_{IR}) - log_{10}(\sigma^2_{ER})');  set(gca,'YTickLabel',' ');
title('Variance in Expression Timing, \sigma^2','FontSize',F); set(gca,'FontSize',F);


subplot(2,2,3);
 hist(n,nb); xlim([-.2,xmax]);
    nh = hist(n,nb); 
  hold on; plot([0,0],[0,4/3*max(nh)],'m--','LineWidth',2);
  
    h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
  ylim([0,4/3*max(nh)]); ylabel('frequency'); 
xlabel('log_{10}(\eta_{IR}) -log_{10}(\eta_{ER})');  set(gca,'YTickLabel',' ');
title('Noise in transcript number, \eta','FontSize',F); 
set(gcf,'color','w'); set(gca,'FontSize',F);


 load compfull_pdist2_data; 
 % NEED TO FIX THIS
 
subplot(2,2,4); 
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
% 
% figure(4); clf; 

 hist(fano,nf); xlim([-1,xmax]);
    nh = hist(n,nf); 
  hold on; plot([0,0],[0,4/3*max(nh)],'m--','LineWidth',2);
  
    h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','b');
  ylim([0,max(nh)]); ylabel('frequency'); 
xlabel('log_{10}(F_{IR}) -log_{10}(F_{ER})');  set(gca,'YTickLabel',' ');
title('Fano Factor, F','FontSize',F); 
set(gcf,'color','w'); set(gca,'FontSize',F);






%%
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

sum(rtv>5/4)/length(rtv)

sum(rtv<4/5)/length(rtv)

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


