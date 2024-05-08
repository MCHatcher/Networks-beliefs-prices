%Wheel example plotter

set(0,'DefaultLineLineWidth',0.75)
figure(1)
for i=1:n 
    count(i) = 1;
%hold on, subplot(3,2,1), plot(Periods, gstack(1:end,i)), 
%title('Individual types'), xlabel('Time')
%hold on, subplot(3,2,3), plot(Periods, Demands(1:end,i)),
%title('Individual demands'), xlabel('Time'), hold on,
%hold on, subplot(3,2,2), plot(Periods, gstack(1:end,i)), 
%title('Individual types'), xlabel('Time')
%hold on, subplot(3,2,4), plot(Periods, Demands(1:end,i)),
%title('Individual demands'), xlabel('Time'), hold on,
end

set(0,'DefaultLineLineWidth',0.75)
hold on, subplot(2,2,2), plot(Periods, gmean, 'b'),
title('Average type'), xlabel('Time'), hold on,
axis([-inf,inf,0,0.5001])
hold on, subplot(2,2,4), plot(Periods, ptild, '--b'), hold on, 
plot(Periods, pcrit, 'r')
title('Price deviation'), xlabel('Time')

figure(2) 
hold on, plot(Periods, ptild, 'b'), plot(Periods, pcrit, 'r')
title('Price deviation'), xlabel('Time')
axis([-1,8,-inf,inf])