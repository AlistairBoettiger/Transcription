
%
%-------------------------------------------------------------------------%
%          Series Decomposition of High Degree Markov Chains              %
%-------------------------------------------------------------------------%
%    
%
% Alistair Boettiger                                   Date Begun: 04/06/09
% Levine Lab                                     Functional Since: 04/13/09
% Functionally Complete                             Last Modified: 04/13/09
%
%
% Notes: 
%   SeriesDecomp.m custom matlab function
%   Take pinchpoint decomposed markov assembly chains and compute first
%   passage times through the chain.  Return first and second moments 
%
% This function is called by:
%   model_CompFull.m 
% 
%
% ------------- Compute LaPlace Transform of Transition Time ------------ % 
function [m1b,m2b] = SeriesDecomp(G)  

n = length(G)-1; % first box of G

% initialize variables
 x=cell(n+1,1); a = sym(zeros(n,1)); b = a; c= a; d=a; Ig=x;
 
% compute eigenvectors of the G_{-sf}^k sub-matrices indexed by k
% x_i = sum_{i~=s,f} x(s)_s (G_{-sf})^{-1}_{ij} G_{jf},  eq (41)
for k = 2:n+1 % k=1
    Ig{k} = inv(G{k}(2:end-1,2:end-1)); % invert G_{-sf}
    x{k} = [-(Ig{k}(:,:)*G{k}(2:end-1,end))',1]; % do sum over i as a matrix product. 
    % note that x^k_f = 1 by definition.  We've shortcut to that here
    %   This gives the vector x^k = {x_i} 
end
x{1} = [1,1]; % x{n+2} = [1,1]; % not needed 
%
% Compute transition probabilities to:   (eval equations (44)-(47))
%     (a) move forward and hit p_{k+1}  (w/o returning to p_k)
%     (b) drop back and hit p_{k-1}  (w/o returning to p_k)
%     (c) move forward but return to p_k (before hitting p_k+1) (excursion into X^{k+1})   
%     (d) drop back into X^{k} but return to p_k (before hitting p_k-1) (excursion into X^k}).                   
%
% note: p_k are the pinch points such that p_0 = s_{k=1} and p_n+1 = f_{k=n}
% We define P as the transition between pinch points so we want to stop the
% chain at each pinch point and ask which of the possabilities happens
% next, (a)-(d).

% sum only over i~=s,f, since x
for k=1:n  % k=0:n-1
    a(k) = -sum( G{k+1}(1,2:end)./( G{k+1}(1,1)+G{k}(end,end)).*x{k+1} );     % (39)
    b(k) = -sum( G{k}(end,1:end-1)./( G{k+1}(1,1)+G{k}(end,end)).*(1-x{k}) ); % (40)
    c(k) = -sum( G{k+1}(1,2:end)./( G{k+1}(1,1)+G{k}(end,end)).*(1-x{k+1}) ); % (41)
    d(k) = -sum( G{k}(end,1:end-1)./( G{k+1}(1,1)+G{k}(end,end)).*x{k} );     % (42)
end

save test2;

% compute matrix P from a and b.
% substitute equations (39) and (40) into (26)
P = sym(zeros(n+1,n+1));
if n>1
P = P + diag(a,1); 
P = P + diag(b,-1);
end
P = P + diag(1-[a;0]-[0;b(1:end)]);
P(end,end-1:end) = [0,1];     % (26) 


% ================================================================ %
% Now compute transition time moments 
%  

%
% E[exp(-lambda tau_{01})]= E[exp(-lambda tau_{->k})] + E[exp(-lambda S^0)]   
% E[exp(-lambda tau_{->k})] = (lambda (lambda I - G^k)^{-1})_{sf}
% E[exp(-lambda S^0)] = -G^k_{ss}/(lambda - G^k_{ss})  
%
% where
% E[exp(-lambda tau_{kk})]= 1/(c_k+d_k)*(c_k*E[exp(-lambda tau_cw^{k+1})]...
%    +d_k*E[exp(-lambda tau_ccw^{k})]   
% E[exp(-lambda tau_cw^{k+1})] = 
%




% %  Laplace Transform approach 
% initialize new variables and data arrays
% syms lambda;
% phi = sym(zeros(1,n)); phiL=phi; phiB=phi; phiF=phi;
% expTf=phi; expTb = phi; expS = phi; expT_cw = phi; expT_ccw =phi;  
% M = sym(zeros(n,n+1)); S = M; Z=M; M2=M;S2=M; 
% 
% 
% % preinvert matrices 
% Gi = cell(1,n+1); 
% for k=1:n+1
%         Gs = G{k}; % define G**
%         f = length(Gs); 
%         Gs(1,:) = sym(zeros(1,f));   
%         Gs(end,:) = sym(zeros(1,f)); 
%         Gi{k} = inv(lambda*eye(f)-Gs);
% end
% 
% 
% for k=1:n % exp[-\lambda \tau_{k,k+1}]   
%     % Additional time to leave the state in the first place    
%         expS(k) = -(G{k+1}(1,1)+G{k}(end,end))/(lambda - G{k+1}(1,1)-G{k}(end,end));                                         
%    %  compute component Laplace transforms for forward jumps
%        % if k<n % computations only valid for upper diagnol k<n                   
%         expTf(k) = -lambda/G{k+1}(1,1)*G{k+1}(1,2:end-1)*(Gi{k+1}(2:end-1,end)./x{k+1}(1:end-1)'); %(43)
% save test;  
%         phiF(k) = expTf(k)*expS(k); % Laplace transform of tau_{k,k+1}
%         phip = diff(phiF(k),lambda);
%         % now we can differentiate laplace transforms to get moments             
%          M(k,k+1) = -P(k,k+1)*subs(phip,lambda,0); % mean of tau_{k,k+1}
%          S(k,k+1) = P(k,k+1)*subs(diff(phip,lambda),lambda,0); % moment2 of tau_{k,k+1}
%          Z(k,k+1) = P(k,k+1)*phiF(k);                                                   
% save test;
%               
%        % computing \tau_{k,k} 
%   % Loop excursions E[\exp(-\lambda \tau_{i,i})]
%       %   clockwise / forward looping
%          expT_cw(k) = -lambda/G{k+1}(1,1)*G{k+1}(1,2:end-1)*(Gi{k+1}(2:end-1,1)./(1-x{k+1}(1:end-1)')); % eq (45)
%        % recall G is indexed from 2 on up to account for G_0-->G_1.  Chain of x is the chain of G
%        % Here we are using matrix product to compute sum as usual 
%        % counterclockwise / backward looping 
%         
%        if k~=1 % zero for dummy matrix at k=1;  
%        expT_ccw(k) = -lambda/G{k}(end,end)*G{k}(end,2:end-1)*(Gi{k}(2:end-1,end)./x{k}(1:end-1)'); % eq (46) % corrected 
%        end
%        phiL(k) = 1/(c(k)+d(k))*(c(k)*expT_cw(k)+d(k)*expT_ccw(k));  % eq (28)
%        phiL(k) = phiL(k)*expS(k); % following defintion in Theorem 3.1
%        % now we can differentiate laplace transforms to get moments
%        phip = diff(phiL(k),lambda);
%        M(k,k) = -P(k,k)*subs(phip,lambda,0); % mean of tau_{0k}
%        S(k,k) =  P(k,k)*subs(diff(phip,lambda),lambda,0); % moment2 of tau_{0k}
%        Z(k,k) = phiL(k)*P(k,k);  
%        % can compute with sum rather than differentiate
% save test       
%        % Backward Excursions
%         if k>1 % computations only valid for lower diagnol 
%             expTb(k) = -lambda/G{k}(end,end)*G{k}(end,2:end-1)*(Gi{k}(2:end-1,1)./(1-x{k}(1:end-1)'));
%             phiB(k) = expTb(k)*expS(k); % Laplace transform of tau_{k,k+1}
%          %  compute component Laplace transforms for backward jumps  
%         
%             phip = diff(phiB(k-1),lambda);
%         % now we can differentiate laplace transforms to get moments
%           M(k,k-1) = -P(k,k-1)*subs(phip,lambda,0); % mean of tau_{k,k+1}
%           S(k,k-1) =  P(k,k-1)*subs(diff(phip,lambda),lambda,0);
%           Z(k,k-1) = phiB(k)*P(k,k-1);
%         end
% end
% save test; 


% % Compute Laplace Transform of tau  
% iZ = inv(eye(n) - Z(:,1:n));
% v = iZ*Z(1:n,n+1);  % laplace transform of tau 
% temp = -diff(v,lambda);
% muT_e = subs(temp,lambda,0);
% sigmaT_e = -subs(diff(temp,lambda),lambda,0);
%
% 
% %  Compute transition moments from moment matrices M and Sigma
%   iP = inv(eye(n)-P(1:end-1,1:end-1));
%    R = [iP;zeros(1,n)];
%   m1T = R*M*ones(n+1,1); 
%   m2T = R*S*ones(n+1,1) + 2*(R*M)^2*ones(n+1,1); 
%   
%   ml1 = m1T(1)
%   ml2 = m2T(1) 
%   save test;


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%  
% Compute Moments for Excursions Between Subchains
%   
 % forward excursions
ET_f=sym(zeros(n,1)); ET2_f = ET_f; ET_cw = ET_f; ET2_cw = ET_f; ET_ccw = ET_f;
ET2_ccw = ET_f; ET_b = ET_f; ET2_b = ET_f; ET_o = ET_f; ET2_o = ET_f; 
ES = ET_f; ES2 = ET_f; M2 = sym(zeros(n,n+1)); S2=M2; 

for k=1:n
 % precompute inversions and modified matrix G**  (Gss)  
    iGb = inv(G{k+1}(2:end-1,2:end-1));
    iG1 = iGb^2;
    iG2 = iGb^3; 
    
 
        % Corollary 3.2 / equation (48)  
 nx = [0,x{k+1}]'; % nx(nx==0) = 1; %denominator
 inx = (nx~=0); inx(1)=0; inx(end) = 0; 
 inxs = inx(2:end-1); 
  
 mx = 1-[0,x{k+1}]';
 imx = (mx~=0); imx(1)=0; imx(end) = 0; 
 imxs = imx(2:end-1); 
  xs = [0,x{k+1}];
  
  save test;
 
    

  
  if isempty( G{k+1}(1,inx)) == 0
     ET_f(k) = -( G{k+1}(1,inx) ./ (G{k+1}(1,1)*xs(inx)) ) * iG1(inxs,inxs)*G{k+1}(inx,end) ; 
     ET2_f(k) = 2*( G{k+1}(1,inx) ./ (xs(inx)*G{k+1}(1,1)) ) * iG2(inxs,inxs)*(G{k+1}(inx,end)) ;
  end
  if isempty( G{k+1}(1,imx)) == 0
     ET_cw(k) =  -( G{k+1}(1,imx) ./ ((1-xs(imx))*G{k+1}(1,1)) ) * iG1(imxs,imxs)*(G{k+1}(imx,1)) ; 
     ET2_cw(k) =  2*( G{k+1}(1,imx) ./ ((1-xs(imx))*G{k+1}(1,1)) ) * iG2(imxs,imxs)*(G{k+1}(imx,1) );
  end
  


 save test; % load test
   i2x = G{k+1}(end,:)'~=0; % don't add terms where the numerator is zero (even if the denom is also zero)    
 icnx = i2x&inx; % combine excluded terms
 icnxs = icnx(2:end-1); % 
  icmx = i2x&imx; % combine excluded terms
 icmxs = icmx(2:end-1); 


if k<n && (isempty( G{k+1}(end,icnx)) == 0)
 % all E[\tau_] are indexed by the index of the starting pinchpoint k+1.  
     ET_ccw(k) = -( G{k+1}(end,icnx) ./ (xs(icnx)*G{k+1}(end,end)) ) * iG1(icnxs,icnxs)*(G{k+1}(icnx,end)) ; 
     ET2_ccw(k) = 2*( G{k+1}(end,icnx) ./ (xs(icnx)*G{k+1}(end,end)) ) *  iG2(icnxs,icnxs)*(G{k+1}(icnx,end)) ; 
end
if k<n && (isempty( G{k+1}(end,icmx)) == 0)
     ET_b(k) =  -( G{k+1}(end,icmx) ./ ((1-xs(icmx))*G{k+1}(end,end)) ) *iG1(icmxs,icmxs)*(G{k+1}(icmx,1)) ; 
     ET2_b(k) =  2*( G{k+1}(end,icmx) ./ ((1-xs(icmx))*G{k+1}(end,end)) ) *iG2(icmxs,icmxs)*(G{k+1}(icmx,1)) ; 
 end
     
 if k~=1
    ET_o(k) = 1./(c(k)+d(k))*(c(k)*ET_cw(k)+d(k)*ET_ccw(k-1));
    ET2_o(k) = 1./(c(k)+d(k))*(c(k)*ET2_cw(k)+d(k)*ET2_ccw(k-1));
 else
    ET_o(k) = 1./(c(k)+d(k))*(c(k)*ET_cw(k)); % ET_ccw not defined for state 1  
    ET2_o(k) = 1./(c(k)+d(k))*(c(k)*ET2_cw(k));
 end

 if (c(k) == 0) && (d(k) == 0)
    ET_o(k) = 0;
    ET2_o(k) = 0;
 end
     
 ES(k) = -1/(G{k+1}(1,1)+G{k}(end,end));  % mean of exponential distribution
 ES2(k) = 2/(G{k+1}(1,1)+G{k}(end,end))^2; % 2nd moment of exponential distribution 
 

 M2(k,k+1) = P(k,k+1)*(ET_f(k)+ES(k)); % mean of tau_{k,k+1}
 S2(k,k+1) = P(k,k+1)*(ET2_f(k)+ES2(k)+2*ET_f(k)*ES(k));
 M2(k,k) = P(k,k)*(ET_o(k)+ES(k)); % mean of tau_{0k}
 S2(k,k) = P(k,k)*(ET2_o(k)+ES2(k)+ 2*ET_o(k)*ES(k) );
 if k>1
    M2(k,k-1) = P(k,k-1)*(ET_b(k-1) + ES(k-1));
    S2(k,k-1) = P(k,k-1)*(ET2_b(k-1)+ES2(k-1)+2*ET_b(k-1)*ES(k-1));   
 end
end
 % Compute momments for tau
 
   iP = inv(eye(n)-P(1:end-1,1:end-1));
   R = [iP;zeros(1,n)];
  m1T = R*M2*ones(n+1,1); 
  m2T = R*S2*ones(n+1,1) + 2*(R*M2)^2*ones(n+1,1);  
  

  
  m1b = m1T(1);
  m2b = m2T(1);
  save test;
 
  
  
%   %
%   ET_f(k) = G{k+1}(1,2:end-1)/G{k+1}(1,1)*iG1*(Gss(2:end-1,end)./x{k+1}(1:end-1)');
%  ET2_f(k) = -2*G{k+1}(1,2:end-1)/G{k+1}(1,1)*iG2*(Gss(2:end-1,end)./x{k+1}(1:end-1)');
%  ET_cw(k) = G{k+1}(1,2:end-1)/G{k+1}(1,1)*iG1*(Gss(2:end-1,1)./nx);
%  ET2_cw(k) = -2*G{k+1}(1,2:end-1)/G{k+1}(1,1)*iG2*(Gss(2:end-1,1)./nx);
%  if G{k+1}(end,end) ~= 0;
%  ET_ccw(k) = G{k+1}(end,2:end-1)/G{k+1}(end,end)*iG1*(Gss(2:end-1,end)./x{k+1}(1:end-1)');
%  ET2_ccw(k) =-2*G{k+1}(end,2:end-1)/G{k+1}(end,end)*iG2*(Gss(2:end-1,end)./x{k+1}(1:end-1)');
%  ET_b(k) = G{k+1}(end,2:end-1)/G{k+1}(end,end)*iG1*(Gss(2:end-1,end)./nx);
%  ET2_b(k) = -2*G{k+1}(end,2:end-1)/G{k+1}(end,end)*iG2*(Gss(2:end-1,end)./nx);
%  end
%  
%  end