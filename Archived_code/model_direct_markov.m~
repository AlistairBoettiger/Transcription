% Basic Fitzimmons Laplace solution for Finite Markov Chains

% Alistair Boettiger                                   Date Begun: 09/22/09
% Levine Lab                                     Functional Since: /09
%                                                   Last Modified: /09
% 
% 

clear all;

% % % parameters
k12 = sym('k12','real');
% k13 = sym('k13','real');
k21 = sym('k21','real');
% k31 = sym('k31','real');
k24 = sym('k24','real'); 
% k34 = sym('k34','real'); 
k41 = sym('k41','real'); 
k45 = sym('k45','real'); 

% % % parameters
% kab = 1;
% kba = 1;
% k12 = 1;
% k21 = 1;
% k23 = 1; 

% Simple shadow model
Gs{1} = zeros(3);
Gs{2} =[[-k12-k12,k12,k12,0,0];
       [k21,-k21-k24,0,k24,0];
       [k21,0,-k21-k24,k24,0];
       [k41,0,0,-k41-k45,k45];
       [0,0,0,0,0]];
% Gs{2} =[[-k12-k13,k12,k13,0,0];
%        [k21,-k21-k24,0,k24,0];
%        [k31,0,-k31-k34,k34,0];
%        [k41,0,0,-k41-k45,k45];
%        [0,0,0,0,0]];
   
  % [m1s, m2s] = SeriesDecomp(Gs); 

% Simple Eloop model
GE{1} = zeros(3);
GE{2} = [[-k12,k12,0,0];
       [k21,-k21-k24,k24,0];
       [k41,0,-k41-k45,k45];
       [0,0,0,0]];
       
   
   % [m1E, m2E] = SeriesDecomp(GE); 
   

   %vars = [k12,k13,k21,k31,k24,k34,k41,k45];
  % vals = [1,1,1,1,1,1,1,1];

   
% grad_mE = jacobian(m1E,vars);
%  grad_sE = jacobian(m2E-m1E^2,vars);
 % grad_nE = jacobian((m2E-m1E^2)/m1E^2,vars);
%   grad_mI = jacobian(m1I,vars);
 % grad_sI = jacobian(m2I-m1I^2,vars);  
%  grad_nI = jacobian((m2I-m1I^2)/m1I^2,vars);
%  grad_dm = jacobian(m1E-m1I,vars);
% grad_ds = jacobian(m2E-m1E^2-(m2I-m1I^2),vars);
% grad_dn = jacobian(((m2E-m1E^2)/m1E^2-(m2I-m1I^2)/m1I^2),vars);
% clear m1E m2E m1I m2I; 


 % grad_dm2 = jacobian(m2E-m1E^2-(m2I-m1I^2),vars);

 %   subs(grad_dn,vars,vals);
   
% % ~~~~~~~~~~~~~~ direct method ~~~~~~~~~~~~~~~~~~~~% 
   syms lambda
     fi=length(Gs{2});
  G= lambda*inv(lambda*eye(fi) - Gs{2});
  phi = G(1,fi);  
  phit = diff(phi,lambda);
  ms_dir = -subs(phit,lambda,0);
  m2s_dir = subs(diff(phit,lambda),lambda,0);
  
 %   
  fe=length(GE{2});
  G = lambda*inv(lambda*eye(fe) - GE{2});
  phi = G(1,fe);  
  phit = diff(phi,lambda);
  mE_dir = -subs(phit,lambda,0);
  m2E_dir = subs(diff(phit,lambda),lambda,0);

%%     


dmean = mE_dir - ms_dir;

Evar = (m2E_dir - mE_dir.^2);
svar = (m2s_dir - ms_dir.^2);
dvar = Evar - svar;
dnoise = Evar/mE_dir^2 - svar/ms_dir^2; 


dMean = zeros(100,1); dVar = zeros(100,1); dNoise = zeros(100,1); 
for j=1:100
       vars = [k12,k21,k24,k41,k45];
       vals = rand(1,5); 
       dMean(j) = subs(dmean,vars,vals);
       dVar(j) =  subs(dvar,vars,vals);
       dNoise(j) = subs(dnoise,vars,vals);
end;

figure(1); clf; hist(dMean);    

figure(1); clf; hist(dVar);    
figure(1); clf; hist(dNoise);    


% k41,k45 << k12, k24;
2/(k12*k24*k45) *[  1/k45 - (k12 + k21 + k24)]

dmean=   
((k21 + k24)*(k41 + k45))/(2*k12*k24*k45)
((k41 + k45))/(2*k12*k45) % as k24--> infty

((k21 + k24)*(k41 + k45))/(2*k12*k24*k45) % as k45 --> infty

