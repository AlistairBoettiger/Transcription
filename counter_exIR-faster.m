
clear all;

load modelFull2b_solns;

pars =   [  k12,      k21,   k23,      k24,     k32,    k35,    k53,    k54,     k42,       k45,     k56,     k65,    k67,       k78,   kab,   kba];         
 vals =  [  .07,     .05,     ones(1,10), 1E-7,    .05,    .07,  1E-10];
    %       k12    k21                  k23    k34

muI_val = subs(m1I,pars,vals); 

muE_val = subs(m1E,pars,vals); 


log(muI_val / muE_val)

