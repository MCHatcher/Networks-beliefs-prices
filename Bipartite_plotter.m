%Bipartite example plotter

set(0,'DefaultLineLineWidth',1)
figure(1)
for i=1:5:6
    grey = (i-1)/6;
hold on, subplot(2,3,3), plot(Periods, gstack(1:end,i), 'Color',[grey grey grey]), 
%title('Individual types'), xlabel('Time')
%hold on, subplot(2,3,3), plot(Periods, Demands(i,1:end), 'b'),
%title('Individual demands'), xlabel('Time'), hold on,
if i>5
hold on, subplot(2,3,3), plot(Periods, gstack(1:end,i), 'Color',[grey grey grey]), 
title('Individual types'), xlabel('Time')
%hold on, subplot(3,2,4), plot(Periods, Demands(i,1:end), 'r'),
%title('Individual demands'), xlabel('Time'), hold on,
end

end
set(0,'DefaultLineLineWidth',1)
subplot(2,3,6), plot(Periods, gmean, 'k'),
title('Average type'), xlabel('Time'), hold on, axis([-inf,inf,0,0.5])
%subplot(2,3,4), plot(Periods, ptild, 'k'),  hold on,
%title('Price deviation'), xlabel('Time')