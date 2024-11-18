%Bipartite_plotter3

set(0,'DefaultLineLineWidth',1), hold on,
subplot(2,2,1), plot(Periods, gmean, '--k'),
title('Average type'), xlabel('Time'), hold on,
subplot(2,2,2), plot(Periods, ptild, '--k'), hold on, 
title('Price deviation'), xlabel('Time')
subplot(2,2,3), 