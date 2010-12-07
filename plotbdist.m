
function plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins,vout)

F =10;

M1Eb = (1-b)*M1E_18 + (b)*M1E_58;
M2Eb = (1-b)*M2E_18 + (b)*M2E_58;
M1Ib = (1-b)*M1I_18 + (b)*M1I_58;
M2Ib = (1-b)*M2I_18 + (b)*M2I_58;

VEb = M2Eb - M1Eb.^2;
VIb = M2Ib - M1Ib.^2; 

NEb = VEb./M1Eb;  % (sigma^2/mu)  ER
NIb = VIb./M1Ib; % (sigma^2/mu)  IR



if vout == 1  
m = real(log2(M1Ib./M1Eb)); mp = m(m>0); mm = m(m<0);
hist(mm,linspace(xmin,xmax,bins)); hold on;
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r');
  hist(mp,linspace(xmin,xmax,bins)); xlim([xmin,xmax]);
xlabel('log_{2}(\eta_{IR}) -log_{2}(\eta_{ER})');  set(gca,'YTickLabel',' ');
%ylabel('frequency'); 
title(['b = ',num2str(b), '  Mean expression time, \mu_{\tau}'],'FontSize',F); 
legend('\mu_{IR} < \mu_{ER}','\mu_{IR} > \mu_{ER}');
end

if vout == 2  
v = real(log2(VIb./VEb)); vp = v(v>0); vm = v(v<0);
hist(vm,linspace(xmin,xmax,bins)); hold on;
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r');
  hist(vp,linspace(xmin,xmax,bins)); xlim([xmin,xmax]);
xlabel('log_{2}(\eta_{IR}) -log_{2}(\eta_{ER})');  set(gca,'YTickLabel',' ');
%ylabel('frequency'); 
title(['b = ',num2str(b), '  Variability in expression time, \sigma^2_{\tau}'],'FontSize',F); 
legend('\sigma^2_{IR} < \sigma^2_{ER}','\sigma^2_{IR} > \sigma^2_{ER}');
end


if vout == 3;
n = real(log2(NIb./NEb)); np = n(n>0); nm = n(n<0);
hist(nm,linspace(xmin,xmax,bins)); hold on;
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r');
  hist(np,linspace(xmin,xmax,bins)); xlim([xmin,xmax]);
xlabel('log_{2}(\eta_{IR}) -log_{2}(\eta_{ER})');  set(gca,'YTickLabel',' ');
%ylabel('frequency'); 
title(['b = ',num2str(b), '  Noise in transcript number, \eta'],'FontSize',F); 
legend('\eta_{IR} < \eta_{ER}','\eta_{IR} > \eta_{ER}');

end