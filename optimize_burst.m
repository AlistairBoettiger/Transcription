
%%  Extension function from multiloop_pdist


Tmax = 600; pts = 100;
B = linspace(0,1,pts);
noiseE = zeros(pts,1);
noiseI = zeros(pts,1);

for i = 1:pts;
    b = B(i); 


    % simulate draws under the ER Model
    NE = zeros(1,Ntot); 
    for k = 1:Ntot % k cycles of Ntot different cells
        
        T = 0; N=0;
        % first draw is always under the 1->8 distribution
            j = 1+round(rand*(length(FE18_dist)-1)); 
            t_mRNA = FE18_dist(j);
            N = N+1;
            T = t_mRNA + T;
        % ~~~~~~~~
        
        while T<Tmax  % during expression interval
            % choose 5->8 or 1->8 distribution based on bernuolli draw b. 
            if rand < b
                bdist = FE18_dist;
            else
                bdist = FE58_dist;
            end
            
            j = 1+round(rand*(length(bdist)-1)); 
            t_mRNA = bdist(j);
            N = N+1;
            T = t_mRNA + T;
        end    
    NE(k) = N;  % number of transcripts made by cell k.  
    end
   noiseE(i) = std(NE)/mean(NE);
 
    
    
    % simulate draws under the IR Model
    NI = zeros(1,Ntot); 
    for k = 1:Ntot % k cycles of Ntot different cells
        
        T = 0; N=0;
        % first draw is always under the 1->8 distribution
            j = 1+round(rand*(length(FI18_dist)-1)); 
            t_mRNA = FI18_dist(j);
            N = N+1;
            T = t_mRNA + T;
        % ~~~~~~~~
        
        while T<Tmax  % during expression interval
            % choose 5->8 or 1->8 distribution based on bernuolli draw b. 
            if rand < b
                bdist = FI18_dist;
            else
                bdist = FI58_dist;
            end
            
            j = 1+round(rand*(length(bdist)-1)); 
            t_mRNA = bdist(j);
            N = N+1;
            T = t_mRNA + T;
        end    
    NI(k) = N;  % number of transcripts made by cell k.  
    end
   noiseI(i) = std(NI)/mean(NI);
   
end

figure(3); clf; plot(B,noiseI,'b'); hold on; plot(B,noiseE,'r');
 plot(B,noiseI./noiseE,'k'); legend('noise_I','noise_E','ratio');
 set(gcf,'color','w'); ylabel('noise'); xlabel('b');