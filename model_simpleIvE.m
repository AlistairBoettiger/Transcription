% Simple model as a control
%

% Alistair Boettiger                                      Date Begun: 04/09
% Levine Lab                                     Functional Since: 06/22/09
%                                               Last Modified: 08/07/09
% 
% Modified 08/07/09 to add explicit representation of ``transcribed''
% state.  
% 

clear all;

% % % parameters
kab = sym('kab','real');
kba = sym('kba','real');
k12 = sym('k12','real'); 
k21 = sym('k21','real'); 
k23 = sym('k23','real'); 
k34 = sym('k34','real'); 


% % % parameters
% kab = 1;
% kba = 1;
% k12 = 1;
% k21 = 1;
% k23 = 1; 

% Simple model
GI{1} = zeros(3);
GI{2} =[[-kab,kab,0,0,0];
        [0,-k12,k12,0,0];
        [0,k21,-k21-k23,k23,0];
        [0,0,0,-k34,k34];
        [0,0,0,0,0 ]];
   
   [m1I, m2I] = SeriesDecomp(GI); 

% Elongation model
GE{1} = zeros(3);
%            1A 2A 3A  1B  2B  3B   4 
GE{2} = [[-k12-kab,k12,kab,0,0,0,0];       % 1A
        [k21,-k21-kab-k23,0,kab,k23,0,0];  % 2A
        [0,0,-kab,0,0,kab,0 ];             % 3A
        [kba,0,0,-kba-k12,k12,0,0];        % 1B
        [0,kba,0,k21,-kba-k21-k23,k23,0];  % 2B
        [0,0,0,0,0,-k34,k34] ;             % 3B
        [0,0,0,0,0,0,0]];                  % 4   
    
  
% GE{2} = [[-k12-kab,k12,kab,0,0,0];
%         [k21,-k21-kab-k23,0,kab,k23,0];
%         [kba,0,-kba-k12,k12,0,0];
%         [0,kba,k21,-kba-k21-k23,0,k23];
%         [0,0,0,0,-kab,kab];
%         [0,0,0,0,0,0]];
       
   [m1E, m2E] = SeriesDecomp(GE); 
   

   vars = [k12,k21,k23,k34,kab,kba];
   vals = [1,1,1,1,1,1];

  % save model_simpIvE m1E m2E m1I m2I; 
   % load model_simpIvE
   
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
   
%       % direct method
%       % % initiation regulation model 
%    syms lambda
%      fi=length(GI{2});
%   G= lambda*inv(lambda*eye(fi) - GI{2});
%   phi = G(1,fi);  
%   phit = diff(phi,lambda);
%   mI_dir = -subs(phit,lambda,0)
%   m2I_dir = subs(diff(phit,lambda),lambda,0);
%   
%   save test2;
%   
%   diff_mean = m1I - mI_dir
%   diff_m2 = m2I - m2I_dir
%     dif_mean = m1b - mI_dir
%   dif_m2 = m2b - m2I_dir
%    
%   subs(dif_mean,[kab,kba,k12,k21,k23],rand(1,5))
%   subs(dif_m2,[kab,kba,k12,k21,k23],rand(1,5))
  
%    
% 
% 
% 
%   
%   fe=length(GE{2});
%   G = lambda*inv(lambda*eye(fe) - GE{2});
%   phi = G(1,fe);  
%   phit = diff(phi,lambda);
%   mE_dir = -subs(phit,lambda,0)
%   m2E_dir = subs(diff(phit,lambda),lambda,0)
%   iG = inv(GE{2}(1:end-1,1:end-1))^2; 
%   iG(1,:)*GE{2}(1:end-1,end);
%   
% 

