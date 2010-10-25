%%                          FK_direct_method
% Alistair Boettiger                                   Date Begun: 10/12/10
% Levine Lab                                        Last Modified: 10/12/10
%
% Description
% Direct Feynman-Kac substitution to solve simple model

clear all;

% % parameters
kab = sym('kab','real');
kba = sym('kba','real');
k12 = sym('k12','real'); 
k21 = sym('k21','real'); 
k23 = sym('k23','real'); 
k34 = sym('k34','real'); 
lambda = sym('lambda','real'); 


% Initiation Regulated Model
%         1A  1B 2B 3B 4
GI =    [[-kab,kab,0,0,0];
        [kba,-k12-kba,k12,0,0];
        [0,k21,-k21-k23,k23,0];
        [0,0,0,-k34,k34];
        [0,0,0,0,0 ]];
  

% Elongation Regulated Model
%            1A 2A 3A  1B  2B  3B   4 
    GE = [[-k12-kab,k12,0,kab,0,0,0];       % 1A
        [k21,-k21-kab-k23,k23,0,kab,0,0];  % 2A
        [0,0,-kab,0,0,kab,0 ];             % 3A
        [kba,0,0,-kba-k12,k12,0,0];        % 1B
        [0,kba,0,k21,-kba-k21-k23,k23,0];  % 2B
        [0,0,kba,0,0,-k34-kba,k34] ;       % 3B
        [0,0,0,0,0,0,0]];                  % 4   
    
    
 % Compute moments using Feynman-Kac formula
 PhiE =    lambda.*inv(lambda*eye(7) - (GE)); 
 dPhiE = diff(PhiE(1,7),lambda);
mE = -subs(dPhiE,lambda,0);


 PhiI =    lambda.*inv(lambda*eye(5) - (GI)); 
 dPhiI = diff(PhiI(1,5),lambda);
mI = -subs(dPhiI,lambda,0);