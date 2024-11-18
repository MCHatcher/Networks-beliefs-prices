% Matlab code to compute network diameter

G = graph(A);
d=distances(G);
diam = max(d)

%Based on MIT Matlab Tools for Network Analysis (2006-2011)
%http://strategic.mit.edu/downloads.php?page=matlab_networks
adj = A; diam=0;
diam_check = [];
for i=1:size(A,1)
    d=simple_dijkstra(A,i);
    diam = max([max(d),diam]);
    diam_check(i) = diam;
end

diam_check

