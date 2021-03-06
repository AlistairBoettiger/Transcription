
%%                               Plot_complexmodel_results
% Alistair Boettiger
% Levine Lab

%% Required Function  
% plotbdist

clear all;


load multiloop_pre_b_data; 
% the correct way to do this, with function calls.  
%%

xmin = -1; xmax = 2; bins = 60;

hout = figure(2); clf; 
subplot(1,4,1); b = 1; 
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins,3)

figure(2); subplot(1,4,2); b = .9;
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins,3)

figure(2); subplot(1,4,3);  b = .3;
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins,3)

figure(2); subplot(1,4,4);  b = 0;
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins,3)

%%

xmin = -1; xmax = 3; bins = 60;

figure(1); clf;  subplot(1,3,1);  b = 0;
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins,1)

figure(1); subplot(1,3,2);  b = 0;
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins,2)

figure(1); subplot(1,3,3);  b = 0;
plotbdist(M1E_18,M2E_18,M1E_58,M2E_58,M1I_18,M2I_18,M1I_58,M2I_58,b,xmin,xmax,bins,3)