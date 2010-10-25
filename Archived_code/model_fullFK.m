
% Direct computation on Full matrix


clear all;

% Specific Model Constants
kab = sym('kab','real');
kba = sym('kba','real');
k01 = sym('k01','real'); 
k10 = sym('k10','real'); 
k12 = sym('k12','real'); 
k21 = sym('k21','real'); 
k23 = sym('k23','real'); 
k24 = sym('k24','real'); 
k32 = sym('k32','real'); 
k35 = sym('k35','real'); 
k42 = sym('k42','real'); 
k45 = sym('k45','real'); 
k53 = sym('k53','real'); % why isn't this in other model
k54 = sym('k54','real'); % ditto
k56 = sym('k56','real'); 
k65 = sym('k65','real'); 
k67 = sym('k67','real'); 
k78 = sym('k78','real'); 
k89 = sym('k89','real');
k81 = sym('k81','real');
k96 = sym('k96','real');
k98 = sym('k98','real');

G = [[-k12,k12,0,0,0,0,0,0,0];
    [k21,-k21-k23-k24,k23,k24,0,0,0,0,0];
    [0,k32,-k32-k35,0,k35,0,0,0,0];
    [0,k42,0,-k42-k45,k45,0,0,0,0];
    [0,0,k53,k54,-k53-k54-k56,k56,0,0,0];
    [0,0,0,0,k65,-k65-k67,k67,0,0];
    [0,0,0,0,0,0,-k78,k78,0];
    [k81,0,0,0,0,0,0,-k81-k89,k89];
    [0,0,0,0,0,k96,0,k98,-k96-k98]];
    
    
  GE = [[G         , kab*eye(9)];
        [kba*eye(9), G         ]];
  
  
  GI = [[-kab,kab,zeros(1,8)];
       [zeros(9,1)  G]];
  GI(2,1:2) = [kba,-k12-kba];  
    
 % Compute the star matrices in the FK formula
  % end state for elongation reg
  fe = 16;  GEs = GE(1:fe,1:fe);
  GEs(16,:) = zeros(1,length(GEs)); 
  
  % end state for initiation reg
  fi = 8; GIs = GI(1:fi,1:fi);
  GIs(8,:) = zeros(1,length(GIs));
  
    
   % direct method
syms lambda
  Gt = lambda*inv(lambda*eye(fi) - GIs);
  phi = Gt(1,fi);  
  phit = diff(phi,lambda);
  m1I = -subs(phit,lambda,0);
  m2I = subs(diff(phit,lambda),lambda,0);
  
  Gt = lambda*inv(lambda*eye(fe) - GEs);
  phi = Gt(1,fe);  
  phit = diff(phi,lambda);
  m1E = -subs(phit,lambda,0);
  m2E = subs(diff(phit,lambda),lambda,0);