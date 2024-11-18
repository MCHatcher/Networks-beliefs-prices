%Simulation for gamma --> infty: stochastic simulations and time to consensus.
%Simulation of star network 'running example'. Last updated: Nov 2024.
%Written by M. Hatcher (m.c.hatcher@soton.ac.uk)

clear; clc; %close all; 

%--------------
%Calibration
%--------------
xbar = 0; %Supply per person
r = 0.04;  phi = 0.4; deltta = 1/phi; dbar = 0.5; %star example
T = 6;  %no. of periods
n = 5;   %no. of agents
pf = ( dbar - xbar/deltta )/ r;  %Steady state fundamental price

%Initialization of network
run star, g_init = [0,0.8,1.1,2,2.8];   %%g_init = [0,1.3,1.6,1.8,2];

%-----------------
%Dividend shocks
%-----------------
%Truncated normal distribution
rng(3)
sigma_d = 0.0005;  
pd = makedist('Normal','mu',0,'sigma',sigma_d);
pd_t = truncate(pd,-dbar,dbar);
shock0 = random(pd_t,1,1);
shock = random(pd_t,T,1);

shock0 = -0.5 + 1e-15;

%----------------
%Initial values
%----------------
gbar_init = sum(g_init)/n;
g0 = g_init';
p0 = pf + 1.25;
plag2 = pf + ((1+r)/gbar_init)^2*(p0-pf);
plag1 = pf + ((1+r)/gbar_init)*(p0-pf); 

%Initial matrices (for storage)
Xweighted = NaN(n,1); Beliefs = NaN(n,T); X = Beliefs; g = Beliefs; p = NaN(T,1); 
dev = p; p_crit = p; gap = p; gbar = p; cap_gain = p; check = p;

%Computation of demands and returns (period 0)
Beliefs_lag1 = (1-g_init)*pf + g_init*plag2;
Xlag = deltta*(Beliefs_lag1+ dbar - (1+r)*plag1);
Beliefs0 = (1-g0)*pf + g0*plag1;
X0 = deltta*(Beliefs0 + dbar - (1+r)*p0);
    
cap_gain0 = p0 + dbar + shock0 - (1+r)*plag1;
gbar0 = mean(g0);
dev0 = p0 - pf;  
p0_crit = gbar_init*(xbar/deltta + shock0)/((1+r)^2 - gbar_init);
gap0 = p0 - pf - p0_crit;
    
%------------------
% Simulation
%------------------ 

for t=1:T  

    if t==1

    for i=1:n

        %For gamma --> infty    
        g_init_adj = g0;
        A_vec = A(i,:)';
        g_init_adj(A_vec==0) = NaN;

        %Find maximal agents
        if cap_gain0*dev0 > 0 
            [row,col] = find(g_init_adj==max(g_init_adj));
        elseif cap_gain0*dev0 < 0
            [row,col] = find(g_init_adj==min(g_init_adj));
        end

         nstar = length(row);   %Number of maximal agents
         V = zeros(1,n);
            for j = row
                 V(j) = 1/nstar;
             end 
            g(i,1) = V*g0;
        
    end

gbar(1) = sum(g(1:n,1))/n;
p(1) = (  dbar  +  (1-gbar(1))*pf + gbar(1)*p0 - xbar/deltta )  /(1+r); 
dev(1) = p(1) - pf;
p_crit(1) = gbar(1)*( (xbar/deltta) + shock(1) ) /((1+r)^2 - gbar(1));

%Computation of demands and returns (period 1)
Beliefs(1:n,1) = (1-g(1:n,1))*pf + g(1:n,1)*p0;
X(1:n,1) = deltta*(Beliefs(1:n,1) + dbar - (1+r)*p(1));
cap_gain(1) = p(1) + dbar + shock(1) - (1+r)*p0;
gap(1) = p(1) - pf - p_crit(1);

    elseif t>=2

%----------------------
% Dates t>=2
%----------------------

        for i=1:n
       
            if t==2
                g_adj = g0;
            elseif t > 2
                g_adj = g(1:n,t-2);
            end
            A_vec = A(i,:)';
            g_adj(A_vec==0) = NaN;

            %Find maximal agents
            if cap_gain(t-1)*dev(t-1) > 0
                [row,col] = find(g_adj==max(g_adj));
            elseif cap_gain(t-1)*dev(t-1) < 0
                [row,col] = find(g_adj==min(g_adj));
            end
                 
            nstar = length(row);   %Number of maximal agents
            V = zeros(1,n);
                for j = row
                    V(j) = 1/nstar;
                end 
            g(i,t) = V*g(1:n,t-1);

        end
     
gbar(t) = sum(g(1:n,t))/n;    
p(t) = (  dbar  +  (1-gbar(t))*pf + gbar(t)*p(t-1) - xbar/deltta )  /(1+r);
dev(t) = p(t) - pf;
p_crit(t) = gbar(t)*( (xbar/deltta) + shock(t) )/((1+r)^2 - gbar(t));
gap(t) = dev(t) - p_crit(t);

%Computation of demands and returns (period t)
Beliefs(1:n,t) = (1-g(1:n,t))*pf + g(1:n,t)*p(t-1);
X(1:n,t) = deltta*( Beliefs(1:n,t) + dbar - (1+r)*p(t) );
Xweighted = (1/n)*X(1:n,t);

cap_gain(t) = p(t) + dbar + shock(t) - (1+r)*p(t-1);
check(t) = sum(Xweighted)-xbar;   %Market clearing

    end

end

%Variables for plotting
Period = 1:T;  Periods = [0; Period'];
ptild = [p0-pf; dev];   %ptild = [plag2-pf; plag1-pf; p0-pf; dev'];
price = [p0; p];  pcrit = [p0_crit; p_crit]; diff = [gap0; gap]; 

gstack = [g0 g]; gmean = [gbar0; gbar];
Belief = [Beliefs0 Beliefs]; Demands = [X0 X];
Cap_gain = [cap_gain0; cap_gain];
zero = zeros(1,T);

run Base_plotter
%run Bipartite_plotter
%run Bipartite_plotter2
%run Wheel_plotter






