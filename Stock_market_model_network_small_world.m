%Simulation of baseline model: nests various cases. Last updated: Aug 22, 2023. 
%Written by M. Hatcher (m.c.hatcher@soton.ac.uk)

clear; clc;

%------------------
%Parameter values
%------------------
xbar = 0; %Supply per person  %r=0.04, 0.0215, 0.030769230769231
r = 0.04; phi = 0.4; deltta = 1/phi; gama = 100; dbar = 0.5;  %star example  
T = 60;  %no. of periods
pf = ( dbar - xbar/deltta )/ r;  %Steady state fundamental price
n = 100;   %no. of agents

%--------------------------------
%Initial matrices (for storage)
%--------------------------------
U_net = NaN(n,n); sum_U = NaN(n,1); rat = sum_U; 
Beliefs = NaN(n,T); X = Beliefs; U = Beliefs; U_tild = U; g= U;
p = NaN(T,1); dev = p; p_crit = p; gap = p; gbar = p; cap_gain = p; check = p;

%-----------------
%Dividend shocks
%-----------------
%Truncated normal distribution
rng(3)
sigma_d = 0.075;  
pd = makedist('Normal','mu',0,'sigma',sigma_d);
pd_t = truncate(pd,-dbar,dbar);
shock0 = random(pd_t,1,1);
shocks = random(pd_t,T,1);
%shocks(1:50) = [0.0410986399143831	-0.0413057671521609	0.00203580857073450	0.0931765273083455	0.0945531516959126	-0.0860633011379160	-0.0612018856272910	-0.122309332290207	-0.0111687217192075	-0.141196202592010	-0.00813113012443586	0.0287257949851132	-0.0440505911904066	0.0342938718450832	0.0172323767287788	-0.148326699700773	0.0111048414404568	-0.0484238836844265	-0.0160831498135564	-0.0429300943759335	0.0378573470121357	-0.0112364234192069	-0.0755561107944348	0.00841148735381180	0.0579942531337106	-0.0379638450142672	-0.0574198206395994	-0.0213457924790365	0.114382938157929	0.148296553687402	0.0334878597850676	0.0973402996534312	0.0763783651126819	-0.0233065065363435	-0.0995419310393515	0.0295909274017404	0.0109122758525702	-0.0265709528163859	-0.0566424977259354	-0.0177379711932591	-0.00584505529507396	-0.0461344904624994	-0.0411116172991344	-0.00796980898177939	0.0812040928770656	0.0163437355619785	-0.0429383372816691	-0.0441645200406739	-0.00854939849898497	-0.0616837688788280]; 

%----------------------------
%Initialization of network
%----------------------------
K = 4; rho = 0.1;
rng(5)
h = WattsStrogatz(n,K,rho);
A = adjacency(h);
A = full(A);
g_init = linspace(0,2,n); rng(5), g_init = 2.06*rand(1,n);
plot(h,'NodeColor','k','EdgeAlpha',0.1);
title('Watts-Strogatz Graph with $n = 100$ nodes, $K = 8$, and $\rho = 0.25$', ...
    'Interpreter','latex')
for i=1:n
A(i,i)=1;
end

%----------------
%Initial values
%----------------
gbar_init = sum(g_init)/n;
g0 = g_init';
p0 = pf + 0.1;
plag2 = pf + ((1+r)/gbar_init)^2*(p0-pf);
plag1 = (  dbar  +  (1-gbar_init)*pf + sum(g_init)/n*plag2  - xbar/deltta )  /(1+r);
plag1_crit = gbar_init*(xbar/deltta)/((1+r)^2 - gbar_init); 
gap_lag1 = (plag1 - pf) - plag1_crit;
p0_crit = sum(g0)/n*( xbar/deltta + shock0 )/((1+r)^2- sum(g0)/n);
ptild0 = p0-pf; gap0 = (p0-pf) - p0_crit;

%----------------------------------------------
%Computation of demands and fitness (period 0)
%----------------------------------------------
Beliefs_lag1 = (1-g_init)*pf + g_init*plag2; 
Xlag =  deltta*(Beliefs_lag1 + dbar - (1+r)*plag1);
Beliefs0 = (1-g0)*pf + g0*plag1;
X0 = deltta*(Beliefs0 + dbar - (1+r)*p0);
U0 = (p0 + dbar + shock0 - (1+r)*plag1)*Xlag;
U_tild0 = exp(gama*U0); 

    for i=1:n
        for j=1:n 
            U_net(i,j) = A(i,j)*U_tild0(j);
        end 
    end
    
cap_gain0 = p0 + dbar + shock0 - (1+r)*plag1;
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
            rat(j) = (U_net(i,j)/ sum_U(i) )*g0(j);  %Rel. fitness of rule j for agent i
            %rat(j) = vpa( (U_net(i,j)/ sum_U(i) )*g0(j), 50);  %For accuracy
        end 
  
            g(i,1) = sum(rat);
    
        end
 
        p(1) = (  dbar  +  (1-sum(g(1:n,1))/n)*pf + sum(g(1:n,1))/n*p0 - xbar/deltta )  /(1+r); 
        dev(1) = p(1) - pf; 
        p_crit(1) = sum(g(1:n,1))/n*( xbar/deltta + shocks(1) )/((1+r)^2 - sum(g(1:n,1))/n);
        gap(1) = dev(1) - p_crit(1);

        %Computation of demands and fitness
        Beliefs(1:n,1) = (1-g(1:n,1))*pf + g(1:n,1)*p0;
        X(1:n,1) = deltta*(Beliefs(1:n,1) + dbar - (1+r)*p(1));
        U(1:n,1) = (p(1) + dbar + shocks(1) - (1+r)*p0)*X0;
        U_tild(1:n,1) = exp(gama*U(1:n,1));
    
        cap_gain(1) = p(1) + dbar + shocks(1) - (1+r)*p0;
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
            %rat(j) = vpa( (U_net(i,j)/ sum_U(i) )*g(j,t-1), 50);   %For accuracy
        end 
  
g(i,t) = sum(rat);     

     end

p(t) = (  dbar  +  (1-sum(g(1:n,t))/n)*pf + sum(g(1:n,t))/n*p(t-1) - xbar/deltta )  /(1+r);
dev(t) = p(t) - pf;
p_crit(t) = sum(g(1:n,t))/n*( xbar/deltta + shocks(t) )/((1+r)^2 - sum(g(1:n,t))/n);
gap(t) = dev(t) - p_crit(t);

%-------------------------------------
%Computation of demands and fitness
%-------------------------------------
Beliefs = (1-g(1:n,t))*pf + g(1:n,t)*p(t-1);
X(1:n,t) = deltta*( Beliefs + dbar - (1+r)*p(t) );
Xweighted = (1/n)*X(1:n,t);
U(1:n,t) = (p(t) + dbar + shocks(t) - (1+r)*p(t-1))*X(1:n,t-1);
U_tild(1:n,t) = exp(gama*U(1:n,t));

for i=1:n
    for j=1:n 
        U_net(i,j) = A(i,j)*U_tild(j,t);     
    end 
end

gbar(t) = sum(g(1:n,t))/n;
cap_gain(t) = p(t) + dbar + shocks(t) - (1+r)*p(t-1);
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

