
%                               model_CompFull_pdist.m
%
% Alistair Boettiger                                   Date Begun: 04/06/09
% Levine Lab                                     Functional Since: 10/09/09
% Functional computational code                     Last Modified: 11/07/09

% Notes:
% compute moment generating function of first passage times for arbitrary
% pinch point decomposed (series) Markov Chains using Matlab's symbolic
% Algebra routines
% Compute laplace transform of complete probability distribution
% Modified 11/07/09 to fix bug k53 and k54 missing (had used k35 and k45)


clear all;
% 
% Specific Model Constants
kab = sym('kab','real');
kba = sym('kba','real');
k12 = sym('k12','real'); 
k21 = sym('k21','real'); 
k23 = sym('k23','real'); 
k24 = sym('k24','real'); 
k32 = sym('k32','real'); 
k35 = sym('k35','real');
k53 = sym('k53','real');
k54 = sym('k54','real');
k42 = sym('k42','real'); 
k45 = sym('k45','real'); 
k56 = sym('k56','real'); 
k65 = sym('k65','real'); 
k67 = sym('k67','real'); 
k78 = sym('k78','real'); 


% Submatrices for initiation regulated system
GI{1} = zeros(3); 

GI{2}=[[-kab,kab,0];
    [kba, -kba-k12, k12];
    [0, k21, -k21]];

GI{3}=[[-k23-k24, k23, k24, 0];
    [k32, -k32-k35, 0, k35];
    [k42, 0, -k42-k45, k45];
    [0, k53, k54, -k53-k54]];

GI{4} =[[-k56, k56, 0, 0];
    [k65, -k65-k67, k67, 0];
    [0, 0, -k78, k78];
    [0, 0, 0, 0]];

[vI] = SeriesDecomp_pdist(GI);  

% Submatrices for elongation regulated system 
p = kba/(kab+kba);

GE{1} = zeros(3); 

GE{2}=[[-k12, k12];
       [k21, -k21]];

GE{3}=[[-k23-k24, k23, k24, 0];
    [k32, -k32-k35, 0, k35];
    [k42, 0, -k42-k45, k45];
    [0, k53, k54, -k53-k54]];


%            5   6  7A  7B  8
GE{4} = [[-k56, k56, 0,  0  ,0];
        [k65,-k65-p*k67-(1-p)*k67, p*k67, (1-p)*k67, 0];
        [0,0,-kab,kab,0 ];  
        [0,0,0,-k78,k78 ];
        [0,0,0,0,0]];

        % either the enhancer has already enabled the promoter to fire, and all
    % gates are open, or the enhancer opens the gate after the chain
    % reaches the paused state.  
    

 [vE] = SeriesDecomp_pdist(GE);  



%   %      1    2    3  4    5  6   7  8    9   10  11  12 13  14  15  16
% vars = [k12,k21,k23,k24,k32,k35,k53,k54,k42,k45,k56,k65,k67,k78,kab,kba];
 
 
 save Model2b_pdist