function IRvER_models

Tmax = 100; % max time at which a single run has a substantial distribution
pts = 1000;
t = linspace(.01,Tmax,pts)'; % time (start >0, finish, num points)

% Rate Constants
kab=1;
kba=1;
k12=1;
k21=1;
k23=1;

% Look up MAPLE computed laplace transforms (MGRI, MGRE), numerically
% compute inverse laplace transform for designated parameter values.  
 psi_I = invlap('MGRI',t,0,1e-4,kab,kba,k12,k21,k23);
psi_E = invlap('MGRE',t,0,1e-4,kab,kba,k12,k21,k23);

%Check normalization
sum(psi_I)*Tmax/pts
sum(psi_E)*Tmax/pts

save test;

load test;
% Plot distributions 
figure(1); clf;  colordef black; set(gcf,'color','k');
 plot(psi_E,'w','linewidth',3);  hold on;
 plot(psi_I,'b','linewidth',3);
 legend('Elongation Regulated','Initiation Regulated');
 set(gca,'FontSize',14);
 xlabel('time t','FontSize',16); 
 ylabel('probability of first transcription','FontSize',16);

 
% use computed distributions to simulate multiple runs
T = 1000;
N = 50000; % samples
cdf_I = cumsum(psi_I*Tmax/pts); 
cdf_E = cumsum(psi_E*Tmax/pts);
figure(3); clf; plot(cdf_I);
mRNAsE = zeros(N,1);
mRNAsI = zeros(N,1);

for n=1:N
TimeI = 0; TauI = 0; NI = 0;
    while TimeI < T   
        % draw a random number 0-1, trace over cdf, find corresponding x.  
        [jnk, TauI] = min((cdf_I-rand).^2);
        TimeI = TimeI+TauI ;
        NI = NI+1;
    end
    TimeE = 0; TauE = 0; NE = 0;
    while TimeE < T 
        % draw a random number 0-1, trace over cdf, find corresponding x.  
        [jnk, TauE] = min((cdf_E-rand).^2);
        TimeE = TimeE+TauE ;
        NE = NE+1;
    end
    mRNAsI(n) = NI;
    mRNAsE(n) = NE;
end
figure(2); clf; colordef white; set(gcf,'color','w'); %colordef black; 
subplot(2,1,1); 
hist(mRNAsI,14);  title('Initiation regulated','FontSize',20); 
xlabel('number of mRNA transcripts synthesized','FontSize',16);
ylabel('frequency','FontSize',16);
subplot(2,1,2); 
hist(mRNAsE,14); title('Elongation regulated','FontSize',20); 
xlabel('number of mRNA transcripts synthesized','FontSize',16);
ylabel('frequency','FontSize',16);
% set(gcf,'color','k'); colordef black; 
mI = mean(mRNAsI)
sI = std(mRNAsI)
mE = mean(mRNAsE)
sE = std(mRNAsE)
nI = sI^2/mI^2
nE = sE^2/mE^2




 