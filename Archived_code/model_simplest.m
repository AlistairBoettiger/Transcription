% Simple model as a control
function model_simplest
clear all;



% Elongation model
GE{1} = zeros(3);
GE{2} = [[-1,1,0];
        [1,-2,1];
        [0,1, -1]];
       
   [m1E, m2E] = SeriesDecomp(GE); 

   T = [[-2,2];
        [1,-1]];
    
    
  [m1r,m2r]=FK(GE{2})
  
   % [m1r,m2r]=FK(T)
    
    
function [m1r,m2r]=FK(T) 
syms lambda;
  f = length(T); 
 T(end,:) = zeros(1,f); 
  G = lambda*inv(lambda*eye(f) - T);
  phi = G(1,f);  
  phit = diff(phi,lambda);
  m1r = -subs(phit,lambda,0);
  m2r = subs(diff(phit,lambda),lambda,0);
