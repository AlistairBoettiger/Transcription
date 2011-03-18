%%                          FK_direct_method
% Alistair Boettiger                                   Date Begun: 10/12/10
% Levine Lab                                        Last Modified: 10/12/10
%
% Description
% Direct Feynman-Kac substitution to solve looping from 

clear all;

% % parameters
kab = sym('kab','real');
kba = sym('kba','real');
k12 = sym('k12','real'); 
k21 = sym('k21','real'); 
k23 = sym('k23','real'); 
k34 = sym('k34','real'); 
lambda = sym('lambda','real'); 



G=[[-k56, k56, 0, 0];
    [k65, -k65-k67, k67, 0];
    [0, 0, -k78, k78];
    [0,0,0,0]]; % kab kba integrated below
    
    
 % Compute moments using Feynman-Kac formula
 PhiE =    lambda.*inv(lambda*eye(7) - (G)); 
 dPhiE = diff(PhiE(1,7),lambda);
mE = -subs(dPhiE,lambda,0);


