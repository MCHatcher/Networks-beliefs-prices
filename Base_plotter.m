%Baseline plots of beliefs, prices, distributions

count = 1:n;

figure(1)
subplot(1,2,1), hold on, plot(Periods, ptild), 
title('Price deviation'), xlabel('Time')
subplot(1,2,2), hold on, plot(Periods, price)
title('Market price'), xlabel('Time')
  
figure(2)
hold on, subplot(1,2,1), plot(Periods, gstack(i,1:end)), 
title('individual chartist weights'), xlabel('Time')
hold on, subplot(1,2,2), plot(Periods, Demands(i,1:end)),
title('Individual demands'), xlabel('Time')

figure(3)
hold on, subplot(1,2,1), plot(count, g_init),
hold on, plot(mean(g_init), 'o', 'MarkerSize',8,'MarkerEdgeColor','red')
title('Initial distribution of g')
hold on, subplot(1,2,2), plot(count, g(:,end))
hold on, plot(mean(g(:,end)), 'o', 'MarkerSize',8,'MarkerEdgeColor','blue')
title('Final distribution of g')