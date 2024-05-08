%Bipartite_plotter3

set(0,'DefaultLineLineWidth',1), hold on,
subplot(1,2,1), plot(Periods, gmean, 'k'),
title('Average type'), xlabel('Time'), hold on,
subplot(1,2,2), plot(Periods, ptild, 'k'),  
title('Price deviation'), xlabel('Time')