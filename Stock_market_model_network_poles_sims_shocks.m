%Simulation of application with fluctutating price and types. Last updated: Aug 22, 2023. 
%Written by M. Hatcher (m.c.hatcher@soton.ac.uk)

%clear; clc;

%Calibration
xbar = 0; %Supply per person
r = 0.04; phi = 0.5; deltta = 1/phi; 
gama = 1e-3; dbar = 0.5;  %star poles example  %gama = 0,5,10,45,75, T=1100
T = 5001;  %no. of periods
samp = T; samp2 = 5; %for plots
pf = ( dbar - xbar/deltta )/ r;  %Steady state fundamental price
n = 3;   %no. of agents

%--------------------------------
%Initial matrices (for storage)
%--------------------------------
U_net = NaN(n,n); sum_U = NaN(n,1); rat = sum_U; 
Beliefs = NaN(n,T); X = Beliefs; U = Beliefs; U_tild = U; g= U; 
p = NaN(T,1); dev = p; p_crit = p; gap = p; gbar = p; cap_gain = p; check = p; 
dev_lag = dev; gbar_lag = dev;

%-----------------
%Dividend shocks
%-----------------
%Truncated normal distribution
%rng(3)
%sigma_d = 0.01;  
%pd = makedist('Normal','mu',0,'sigma',sigma_d);
%pd_t = truncate(pd,-dbar,dbar);
%shocks = random(pd_t,T,1);
shocks = zeros(T,1);

%----------------------------
%Initial types and network
%----------------------------
g_init = zeros(1,n);
eps = 0.001;
g_init(1,1) = 1; g_init(1,2) = 1+r; g_init(1,3) = (1+r+eps)^2;
A = [1 0 0; 1 1 1; 0 0 1]; 

%----------------
%Initial values
%----------------
gbar_init = sum(g_init)/n;
p0 = pf + 1; %1
plag2 = pf + ((1+r)/gbar_init)^2*(p0-pf);
plag1 = (  dbar  +  (1-sum(g_init)/n)*pf + sum(g_init)/n*plag2  - xbar/deltta )  /(1+r);
plag1_crit = sum(g_init)/n*( (xbar/deltta) )/((1+r)^2 - sum(g_init)/n); 
gap_lag1 = (plag1 - pf) - plag1_crit;

g0 = g_init';
p0_check = (  dbar  +  (1-sum(g0)/n)*pf + sum(g0)/n*plag1  - xbar/deltta )  /(1+r);
p0_crit = sum(g0)/n*(xbar/deltta)/((1+r)^2- sum(g0)/n);
ptild0 = p0-pf;  gap0 = (p0-pf) - p0_crit;
       
%----------------------------------------------
%Computation of demands and fitness (period 0)
%----------------------------------------------
Beliefs_lag1 = (1-g_init)*pf + g_init*plag2;
Xlag = deltta*(Beliefs_lag1 + dbar - (1+r)*plag1);
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
p_crit(1) = sum(g(1:n,1))/n*( (xbar/deltta) + shocks(1) ) /((1+r)^2 - sum(g(1:n,1))/n);
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
        end 
  
    g(i,t) = sum(rat);        

    end

p(t) = (  dbar  +  (1-sum(g(1:n,t))/n)*pf + sum(g(1:n,t))/n*p(t-1) - xbar/deltta )  /(1+r);
dev(t) = p(t) - pf;
dev_lag(t) = dev(t-1);
p_crit(t) = sum(g(1:n,t))/n*( (xbar/deltta) + shocks(t) ) /((1+r)^2 - sum(g(1:n,t))/n);
gap(t) = dev(t) - p_crit(t);

%Computation of indvidual demands and fitness
Beliefs(1:n,t) = (1-g(1:n,t))*pf + g(1:n,t)*p(t-1);
X(1:n,t) = deltta*( Beliefs(1:n,t) + dbar - (1+r)*p(t) );
Xweighted = (1/n)*X(1:n,t);
U(1:n,t) = (p(t) + dbar + shocks(t) - (1+r)*p(t-1))*X(1:n,t-1);
U_tild(1:n,t) = exp(gama*U(1:n,t));

    for i=1:n
        for j=1:n 
            U_net(i,j) = A(i,j)*U_tild(j,t);     
        end 
    end

gbar(t) = sum(g(1:n,t))/n;
gbar_lag(t) = gbar(t-1);
cap_gain(t) = p(t) + dbar + shocks(t) - (1+r)*p(t-1);
check(t) = sum(Xweighted)-xbar;   %Market clearing

    end

end

gbar_end = gbar(end);

%Variables for plotting
Period = 1:T;
Periods = [0; Period'];
%ptild = [plag2-pf; plag1-pf; p0-pf; dev'];
ptild = [p0-pf; dev];
price = [p0; p];
pcrit = [p0_crit; p_crit]; 
diff = [gap0; gap]; 

gstack = [ g0 g];
Belief = [Beliefs0 Beliefs];
Demands = [X0 X];
gmean = [gbar0; gbar];
Cap_gain = [cap_gain0; cap_gain];
zero = zeros(1,length(Periods));

set(0,'DefaultLineLineWidth',1)

figure(1)
plot(Periods(1:samp),ptild(1:samp),'--k'), title('Price deviation'), xlabel('Time'), hold on,
axis([-inf,inf,0,3000])

%figure(2)
%subplot(1,2,1), plot(Periods(1:samp),ptild(1:samp),'--k'), title('Price deviation'), xlabel('Time'), hold on,
%axis([-inf,inf,0,inf])
%subplot(1,2,2), plot(Periods(1:samp2), ptild(1:samp2),'--k'), title('Price deviation (zoomed in)'), xlabel('Time'), hold on

%figure(3)
%subplot(1,2,2), plot(ptild(1:samp) - pcrit(1:samp),ptild(1:samp)), title('Price deviation vs critical price (\gamma=2.5)')

%figure(4)
%subplot(1,2,1), plot(Periods(1:samp),pcrit(1:samp)), hold on, plot(Periods(1:samp),ptild(1:samp)), title('Price deviation vs critical price (\gamma=2.5)')
