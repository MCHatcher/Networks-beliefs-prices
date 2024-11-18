%Bipartite_plotter2
%For bipartite example, left panel

set(0,'DefaultLineLineWidth',1)
%subplot(2,3,1)
%for i=[1,6]
%    grey = i/6;
%    plot(Periods,gstack(:,i),'Color',[grey grey grey]), hold on,
%end
%title('Individual types'), xlabel('Periods')
subplot(1,2,1), plot(Periods,gmean,'--k'), hold on,
axis([-inf,inf,-inf,inf])
title('Average type'), xlabel('Time')