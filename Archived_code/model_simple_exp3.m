% load simpIvE
% max_dn = max(abs(grad_dn));
% max_En = max(abs(grad_En));
% max_In = max(abs(grad_In));
% max_dm = max(abs(grad_dm));
% max_Em = max(abs(grad_Em));
% max_Im = max(abs(grad_Im));
% max_Se = max(abs(grad_Se));
% max_Si = max(abs(grad_Si));
% max_ds = max(abs(grad_ds));
% 
% std_dn = std(grad_dn);
% std_En = std(grad_En);
% std_In = std(grad_In);
% std_dm = std(grad_dm);
% std_Em = std(grad_Em);
% std_Im = std(grad_Im);
% std_Se = std(grad_Se);
% std_Si = std(grad_Si);
% std_ds = std(grad_ds);
max_dn = max(log(abs(grad_dn)));
max_En = max(log(abs(grad_En)));
max_In = max(log(abs(grad_In)));
max_dm = max(log(abs(grad_dm)));
max_Em = max(log(abs(grad_Em)));
max_Im = max(log(abs(grad_Im)));
max_Se = max(log(abs(grad_Se)));
max_Si = max(log(abs(grad_Si)));
max_ds = max(log(abs(grad_ds)));
% 

mean_dn = mean(log(abs(grad_dn)));
mean_En = mean(log(abs(grad_En)));
mean_In = mean(log(abs(grad_In)));
mean_dm = mean(log(abs(grad_dm)));
mean_Em = mean(log(abs(grad_Em)));
mean_Im = mean(log(abs(grad_Im)));
mean_Se = mean(log(abs(grad_Se)));
mean_Si = mean(log(abs(grad_Si)));
mean_ds = mean(log(abs(grad_ds)));

std_dn = std(log(grad_dn));
std_En = std(log(grad_En));
std_In = std(log(grad_In));
std_dm = std(log(grad_dm));
std_Em =std(log(grad_Em));
std_Im = std(log(grad_Im));
std_Se = std(log(grad_Se));
std_Si =std(log(grad_Si));
std_ds = std(log(grad_ds));


% m=6;
% figure(2); clf;% pie(std_In);
% subplot(m,3,1); pie(max_Em); 
% subplot(m,3,2); pie(max_Im);
% subplot(m,3,3); pie(max_dm);
% subplot(m,3,4); pie(std_Em);
% subplot(m,3,5); pie(std_Im);
% subplot(m,3,6); pie(std_dm);
% subplot(m,3,7); pie(max_En);
% subplot(m,3,8); pie(max_In);
% subplot(m,3,9); pie(max_dn);
% subplot(m,3,10); pie(std_En);
% subplot(m,3,11); pie(std_In);
% subplot(m,3,12); pie(std_dn);
% subplot(m,3,13); pie(max_Se); 
% subplot(m,3,14); pie(max_Si);
% subplot(m,3,15); pie(max_ds);
% subplot(m,3,16); pie(std_Se);
% subplot(m,3,17); pie(std_Si);
% subplot(m,3,18);pie(std_ds);%  legend('kab','kba','k12','k21','k23');% ,'Location','SouthOutside');

figure(3); clf;  colormap(summer(3));
subplot(3,1,1); 
    bar(([std_Im;std_Em;std_dm]'));
    bar(([max_Im;max_Em;max_dm]'));
     set(gca,'XtickLabel',...
         str2mat('kab','kba','k12','k21','k23'),...
         'XTick',1:5,'fontsize',12);
     ylabel('log(abs( gradient ))');
     title('expression speed, max sensitivity'); 
subplot(3,1,2);
    bar(([std_Si;std_Se;std_ds]'));
    bar(([max_Si;max_Se;max_ds]'));
     set(gca,'XtickLabel',...
         str2mat('kab','kba','k12','k21','k23'),...
         'XTick',1:5,'fontsize',12);
       ylabel('log(abs( gradient ))');
     title('expression variance, max sensitivity');
subplot(3,1,3); 
    bar(([std_In;std_En;std_dn]'));
    bar(([max_In;max_En;max_dn]'));
     set(gca,'XtickLabel',...
         str2mat('kab','kba','k12','k21','k23'),...
         'XTick',1:5,'fontsize',12);
        ylabel('log(abs( gradient ))');
     title('expression noise power, max sensitivity');
     legend('Initiation Model','Elongation Model','Difference','Location','Best');

 figure(4); clf;  colormap(summer(3));
subplot(3,1,1); 
    bar(([std_Im;std_Em;std_dm]'));
     set(gca,'XtickLabel',...
         str2mat('kab','kba','k12','k21','k23'),...
         'XTick',1:5,'fontsize',12);
       ylabel('log(abs( gradient ))');
     title('expression speed, stdev sensitivity');
subplot(3,1,2);
    bar(([std_Si;std_Se;std_ds]'));
     set(gca,'XtickLabel',...
        str2mat('kab','kba','k12','k21','k23'),...
         'XTick',1:5,'fontsize',12);
       ylabel('log(abs( gradient ))');
     title('expression variance, stdev sensitivity');
subplot(3,1,3); 
    bar(([std_In;std_En;std_dn]'));
     set(gca,'XtickLabel',...
         str2mat('kab','kba','k12','k21','k23'),...
         'XTick',1:5,'fontsize',12);
      ylabel('log(abs( gradient ))');
     title('expression noise power, stdev sensitivity');
        legend('Initiation Model','Elongation Model','Difference','Location','Best');
     
     figure(5); clf;  colormap(summer(3));
subplot(3,1,1); 
    bar(([mean_Im;mean_Em;mean_dm]'));
     set(gca,'XtickLabel',...
        str2mat('kab','kba','k12','k21','k23'),...
         'XTick',1:5,'fontsize',12);
      ylabel('log(abs( gradient ))');
     title('expression speed, mean sensitivity'); 
subplot(3,1,2);
    bar(([mean_Si;mean_Se;mean_ds]'));
     set(gca,'XtickLabel',...
         str2mat('kab','kba','k12','k21','k23'),...
         'XTick',1:5,'fontsize',12);
      ylabel('log(abs( gradient ))');
     title('expression variance, mean sensitivity');
subplot(3,1,3); 
    bar(([mean_In;mean_En;mean_dn]'));
     set(gca,'XtickLabel',...
        str2mat('kab','kba','k12','k21','k23'),...
         'XTick',1:5,'fontsize',12);
      ylabel('log(abs( gradient ))');
     title('expression noise power, mean sensitivity');
    legend('Initiation Model','Elongation Model','Difference','Location','Best');
 
