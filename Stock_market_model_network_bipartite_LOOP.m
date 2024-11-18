%Simulation of bipartite example with loop. Last updated: March 27 2024. 
%Written by M. Hatcher (m.c.hatcher@soton.ac.uk)

clear; clc;

%Parameter values
xbar = 0; %Supply per person
r = 0.04; phi = 0.5; deltta = 1/phi; gama = 20; dbar = 0.1;  %bipartite example  
%r = 0.04; phi = 1; deltta = 1/phi; gama = 0; dbar = 0.02;  %bubble example  %gama = 1.33,1,0
%r = 0.04; phi = 0.4; deltta = 1/phi; gama = 2; dbar = 0.5;  %wheel example  

T = 100;  %no. of periods
pf = ( dbar - xbar/deltta )/ r;  %Steady state fundamental price
n = 20;   %no. of agents

Ngama = 200;
gama_stack = linspace(0,100,Ngama);
consensus = NaN(Ngama,1); variance =  consensus; max_var = consensus;

for m=1:Ngama
    
    gama = gama_stack(m);

%--------------------------------
%Initial matrices (for storage)
%--------------------------------
U_net = NaN(n,n); sum_U = NaN(n,1); rat = sum_U; 
Beliefs = NaN(n,T); X = Beliefs; U = Beliefs; U_tild = U; g = U;
p = NaN(T,1); dev = p; p_crit = p; gap = p; gbar = p; cap_gain = p; check = p;

%----------------------------
%Initialization of network
%----------------------------
%A = rand(n,n)>.0.5; %Random network structure 
%run bipartite; g_init = [1,1,1,1,1,0,0,0,0,0];  %Initial distribution of g
K = 16; % No. of chartists
run bipartite_K; g_init = [ones(1,K),zeros(1,20-K)];  %Initial distribution of g

%run wheel %run star  
%A(1,n)=0; A(1,2)=0; A(n,n-1)=0; A(n,1)=0;

%----------------
%Initial values
%----------------
gbar_init = sum(g_init)/n;
g0 = g_init';

p0 = pf + 0.25;
plag2 = pf + ((1+r)/gbar_init)^2*(p0-pf);
plag1 = (  dbar  +  (1-sum(g_init)/n)*pf + sum(g_init)/n*plag2  - xbar/deltta )  /(1+r);
plag1_crit = sum(g_init)/n*(xbar/deltta)/((1+r)^2 - sum(g_init)/n); 
gap_lag1 = (plag1 - pf) - plag1_crit;

p0_check = (  dbar  +  (1-sum(g0)/n)*pf + sum(g0)/n*plag1  - xbar/deltta )  /(1+r);
p0_crit = sum(g0)/n*(xbar/deltta)/((1+r)^2- sum(g0)/n);
ptild0 = p0-pf;
gap0 = (p0-pf) - p0_crit;
Beliefs_lag1 = (1-g_init)*pf + g_init*plag2;
Xlag = deltta*(Beliefs_lag1 + dbar - (1+r)*plag1);

%----------------------------------------------
%Computation of demands and fitness (period 0)
%----------------------------------------------
Beliefs0 = (1-g0)*pf + g0*plag1;
X0 = deltta*(Beliefs0 + dbar - (1+r)*p0);
U0 = (p0 + dbar - (1+r)*plag1)*Xlag;
U_tild0 = exp(gama*U0); 

    for i=1:n
        for j=1:n 
            U_net(i,j) = A(i,j)*U_tild0(j);
        end 
    end
    
    cap_gain0 = p0 + dbar - (1+r)*plag1;
    cap_gain_lag1 = plag1 + dbar - (1+r)*plag2;
    gbar0 = mean(g0);
    
%------------------
% Simulation
%------------------ 
for t=1:T

    if t==1 
   
        for i=1:n
      
        sum_U(i) = sum(U_net(i,1:n));

            for j=1:n
                rat(j) = (U_net(i,j)/ sum_U(i) )*g0(j);   %Rel. fitness of rule j for agent i
            end 
  
            g(i,1) = sum(rat);
        end
 
p(1) = (  dbar  +  (1-sum(g(1:n,1))/n)*pf + sum(g(1:n,1))/n*p0 - xbar/deltta )  /(1+r); 
dev(1) = p(1) - pf;
p_crit(1) = sum(g(1:n,1))/n*(xbar/deltta)/((1+r)^2 - sum(g(1:n,1))/n);
gap(1) = dev(1) - p_crit(1);

for i=1:n
    Beliefs(i,1) = (1-g(i,1))*pf + g(i,1)*p0;
    X(i,1) = deltta*(Beliefs(i,1) + dbar - (1+r)*p(1));
    U(i,1) = (p(1) + dbar - (1+r)*p0)*X0(i);
    U_tild(i,1) = exp(gama*U(i,1));
end
    
    cap_gain(1) = p(1) + dbar - (1+r)*p0;
    gbar(1) = sum(g(1:n,1))/n;
   
for i=1:n
    for j=1:n 
        U_net(i,j) = A(i,j)*U_tild(j,1);     
    end 
end

    elseif t>=2 

%----------------------
% Dates t>=2
%----------------------  
     
    for i=1:n
        sum_U(i) = sum(U_net(i,1:n));

        for j=1:n
            rat(j) = (U_net(i,j)/ sum_U(i) )*g(j,t-1);   %Rel. fitness of rule j for agent i
        end 
  
    g(i,t) = sum(rat);     

    end

p(t) = (  dbar  +  (1-sum(g(1:n,t))/n)*pf + sum(g(1:n,t))/n*p(t-1) - xbar/deltta )  /(1+r);
dev(t) = p(t) - pf;
p_crit(t) = sum(g(1:n,t))/n*(xbar/deltta)/((1+r)^2 - sum(g(1:n,t))/n);
gap(t) = dev(t) - p_crit(t);

%Computation of indvidual demands and fitness
Beliefs(1:n,t) = (1-g(1:n,t))*pf + g(1:n,t)*p(t-1);
X(1:n,t) = deltta*( Beliefs(1:n,t) + dbar - (1+r)*p(t) );
Xweighted = (1/n)*X(1:n,t);
U(1:n,t) = (p(t) + dbar - (1+r)*p(t-1))*X(1:n,t-1);
U_tild(1:n,t) = exp(gama*U(1:n,t));


    for i=1:n
        for j=1:n 
            U_net(i,j) = A(i,j)*U_tild(j,t);     
        end 
    end

gbar(t) = sum(g(1:n,t))/n;
cap_gain(t) = p(t) + dbar - (1+r)*p(t-1);
check(t) = sum(Xweighted)-xbar;   %Market clearing

    end

end

gbar_end = gbar(end);

Period = 1:T;
Periods = [0; Period'];
ptild = [p0-pf; dev];  %ptild = [plag2-pf; plag1-pf; p0-pf; dev'];
price = [p0; p];
pcrit = [p0_crit; p_crit]; 
diff = [gap0; gap]; 

gstack = [g0'; g'];
Belief = [Beliefs0'; Beliefs'];
Demands = [X0'; X'];
gmean = [gbar0; gbar];
Cap_gain = [cap_gain0; cap_gain];
zero = zeros(1,length(Periods));

consensus(m) = gbar(end);
variance(m) = var(g(1:n,end));
max_var(m) = max(variance(m));

end

%run Base_plotter %run Bipartite_plotter
%run Bubble_plotter %run Wheel_plotter

var_max = max(max_var)
grey = 0;

%Loop plotter
set(0,'DefaultLineLineWidth',1)
figure(1)
subplot(1,2,2), plot(gama_stack, consensus, '','Color',[grey grey grey]),  hold on,
title('Consensus type'), xlabel('\gamma'), axis([-inf,max(gama_stack),0,0.625]), hold on,
