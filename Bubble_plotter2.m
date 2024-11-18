%Bubble_plotter2

set(0,'DefaultLineLineWidth',0.75)
figure(1)
subplot(1,2,1), plot(Periods, ptild,'--k'),  hold on,
title('Price deviation'), xlabel('Time'), axis([-inf,40,0,inf]), hold on,
subplot(1,2,2), plot(Periods, gmean, '--k'),  hold on,
title('Average type'), xlabel('Time'), axis([-inf,40,-inf,inf])

figure(2)
%subplot(2,3,1), plot(Periods, ptild,'--k'),  hold on,
title('Price deviation'), xlabel('Time'), axis([-inf,inf,0,inf]), hold on,
%subplot(2,3,4), plot(Periods, gmean, '--k'),  hold on,
title('Average type'), xlabel('Time'), axis([-inf,inf,0,inf])