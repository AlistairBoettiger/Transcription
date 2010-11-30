
function plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins)

F =10;

M1Eb = b*M1E_18 + (1-b)*M1E_58;
M2Eb = b*M2E_18 + (1-b)*M2E_58;
M1Ib = b*M1I_18 + (1-b)*M1I_58;
M2Ib = b*M2I_18 + (1-b)*M2I_58;

NEb = (M2Eb - M1Eb.^2)./M1Eb;  % (sigma^2/mu)  ER
NIb = (M2Ib - M1Ib.^2)./M1Ib; % (sigma^2/mu)  IR

n = real(log2(NIb./NEb)); np = n(n>0); nm = n(n<0);
hist(nm,linspace(xmin,xmax,bins)); hold on;
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','none','FaceColor','r');
  hist(np,linspace(xmin,xmax,bins)); xlim([xmin,xmax]);
xlabel('log_{2}(\eta_{IR}) -log_{2}(\eta_{ER})');  set(gca,'YTickLabel',' ');
%ylabel('frequency'); 
title(['b = ',num2str(b), '  Noise in transcript number, \eta'],'FontSize',F); 
legend('\eta_{IR} < \eta_{ER}','\eta_{IR} > \eta_{ER}');