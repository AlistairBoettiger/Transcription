

T=500;  
tf = zeros(T,2);
tf2 = zeros(T,2);
p = [1,0];
p2 = [1,0];



for t=1:T-1
    figure(2); clf;
    tf(t+1,:) = tf(t,:) + .3*(rand(1,2)-.5);
    subplot(1,2,1); 
    if tf(t+1,1) > 2; tf(t+1,1) = 2; end; 
      if tf(t+1,2) > 2; tf(t+1,2) = 2; end; 
          if tf(t+1,1) < 0; tf(t+1,1) = 0; end; 
      if tf(t+1,2) < 0; tf(t+1,2) = 0; end; 
    
    plot(tf(t+1,1), tf(t+1,2),'b.','MarkerSize',30);
    hold on; plot(p(1),p(2),'r.','MarkerSize',30);
    xlim([0,2]); ylim([0,2]);
    
   
     tf2(t+1,:) = tf2(t,:) + .3*(rand(1,2)-.5);
    subplot(1,2,2); 
    if tf2(t+1,1) > 2; tf2(t+1,1) = 2; end; 
      if tf2(t+1,2) > 2; tf2(t+1,2) = 2; end; 
          if tf2(t+1,1) < 0; tf2(t+1,1) = 0; end; 
      if tf2(t+1,2) < 0; tf2(t+1,2) = 0; end; 
    
    plot(tf2(t+1,1), tf2(t+1,2),'b.','MarkerSize',30);
    hold on; plot(p(1),p(2),'r.','MarkerSize',30);
    xlim([0,2]); ylim([0,2]);
    
    d1 = sqrt( (tf(t+1,1)-p(1)).^2 + (tf(t+1,2)-p(2)).^2 );
    if d1<.1; 
    subplot(1,2,1); plot(tf(t+1,1), tf(t+1,2),'y*','MarkerSize',30);
    tf(t+1,:) = [NaN NaN]; 
    end;
    
        d2 = sqrt( (tf2(t+1,1)-p(1)).^2 + (tf2(t+1,2)-p(2)).^2 );
    if d2<.1; 
    subplot(1,2,2); plot(tf2(t+1,1), tf2(t+1,2),'y*','MarkerSize',30);
    tf2(t+1,:) = [NaN,NaN]; 
    end;
    
    
    pause(.001); 
    
end
