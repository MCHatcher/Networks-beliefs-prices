%Watts_Strogatz_plotter
%We omit self-loops for the plots here

clear
close all

deg = 2;
rho = 0.2;
n = 100;

rng(15)
h = WattsStrogatz(n,deg,rho);
h1 = WattsStrogatz(n,2*deg,rho);
h2 = WattsStrogatz(n,2*deg,2.5*rho);

figure(1)
subplot(2,2,1)
plot(h,'NodeColor','k','EdgeAlpha',0.5,'EdgeColor',[0.5 0.5 0.5],'Layout','circle');
title('Graph: n = 100, K = 2, \rho = 0.20')

subplot(2,2,2)
A = adjacency(h);
A = full(A);
%for i=1:n A(i,i)=1; end
spy(A,'k',7)
title('Adjacency Matrix: n = 100, K = 2, \rho = 0.20')


subplot(2,2,3)
%histogram(degree(h),'BinMethod','integers','FaceAlpha',0.9);
histogram(degree(h),'FaceColor',[0.7,0.7,0.7])
hold on
%A = adjacency(h1); A = full(A);
%for i=1:n A(i,i)=1; end, h1 = graph(A);
%histogram(degree(h1),'BinMethod','integers','FaceAlpha',0.9);
histogram(degree(h1),'FaceColor',[0,0,0])
title('Degree distribution: n = 100, \rho = 0.2')
xlabel('Degree of node')
ylabel('Number of nodes')
legend('K = 2','K = 4','Location','NorthWest')

subplot(2,2,4)
%histogram(degree(h),'BinMethod','integers','FaceAlpha',0.9);
hold on
%histogram(degree(h2),'BinMethod','integers','FaceAlpha',0.9);
%A = adjacency(h2); A = full(A);
%for i=1:n A(i,i)=1; end h2 = graph(A);
histogram(degree(h),'FaceColor',[0.7,0.7,0.7])
histogram(degree(h2),'FaceColor',[0,0,0])
title('Degree distribution: n = 100','Interpreter','latex')
xlabel('Degree of node')
ylabel('Number of nodes')
legend('K = 2, \rho = 0.2','K = 4, \rho = 0.5','Location','NorthWest')



