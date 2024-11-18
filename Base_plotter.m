%Baseline plots of beliefs, prices, distributions

count = 1:n;
set(0,'DefaultLineLineWidth',1)

figure(1)
subplot(1,2,1), hold on, plot(Periods, ptild), 
title('Price deviation'), xlabel('Time')
subplot(1,2,2), hold on, plot(Periods, price)
title('Market price'), xlabel('Time')

for i=1:n
    figure(2)
    hold on, subplot(1,2,1), plot(Periods, gstack(i,1:end)), 
    title('Individual types'), xlabel('Time')   
    hold on, subplot(1,2,2), plot(Periods, Demands(i,1:end)),
    title('Individual demands'), xlabel('Time')
end

figure(3)
hold on, subplot(1,2,1), plot(count, g_init),
hold on, plot(mean(g_init), 'o', 'MarkerSize',8,'MarkerEdgeColor','red')
title('Initial distribution of g')
hold on, subplot(1,2,2), plot(count, g(:,end))
hold on, plot(mean(g(:,end)), 'o', 'MarkerSize',8,'MarkerEdgeColor','blue')
title('Final distribution of g')

T_plot = T;
T_plot = T_plot*(T_plot<=T) + T*(T_plot>T);

for i=1:n
    figure(4)
    hold on, subplot(1,2,1), plot(Periods(1:T_plot+1), gstack(i,1:T_plot+1)), 
    title('Individual types'), xlabel('Time')   
end
%yline(gmean(1))
subplot(1,2,2), hold on, plot(Periods(1:T+1), ptild(1:T+1)), hold on, 
title('Price deviation'), xlabel('Time')

