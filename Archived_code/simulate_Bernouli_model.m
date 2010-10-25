

N = 20;
K = 5000;
FAC1 = zeros(1,K); 
FAC2 = zeros(1,K); 

for k = 1:K
    S1 = true(1,N);
    F1 = false(1,N);

     S2 = true(1,N);
    F2 = false(1,N);
    
    for t=1:2
        jump = rand(1,N)>mean(A);
        S1(jump) = false;
        F1(jump) = true;

        fall = rand(1,N)<0;
        F1(fall) = false;
        S1(fall) =  true;
        
        jump = rand(1,N)>mean(B);
        S2(jump) = false;
        F2(jump) = true;

        fall = rand(1,N)<0;
        F2(fall) = false;
        S2(fall) =  true;
    end

    FAC1(k) = sum(S1)/N;
    FAC2(k) = sum(S2)/N;
end

figure(1); clf;
hist(FAC1); 
hold on;
h =  findobj(gca,'type','patch');
set(h,'FaceColor','r','EdgeColor','none'); 
hist(FAC2) 
alpha(.6);
xlim([0,1]);