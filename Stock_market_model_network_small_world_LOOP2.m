%Simulation of small world network: final application. Loops for consensus or rewiring prob. Last updated: Nov 2024. 
%Written by M. Hatcher (m.c.hatcher@soton.ac.uk)

clear; clc;

%------------------
%Parameter values
%------------------
xbar = 0; %Supply per person  %r=0.04
r = 0.04; phi = 0.4; deltta = 1/phi; gama = 180; dbar = 1;  %star example  
T = 800;  %no. of periods
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
sigma_d = 0.075;  %1/3 %0.075 %0.15   
pd = makedist('Normal','mu',0,'sigma',sigma_d);
pd_t = truncate(pd,-dbar,dbar);
shock0 = random(pd_t,1,1);
shocks = random(pd_t,T,1);

num_K = 17;  %maximum mean degree
max_var = NaN(num_K,1);
dummy_K = zeros(num_K,1); 
deg_K = zeros(num_K,1);

N_rand_net = 30;
consensus = NaN(N_rand_net,num_K);
rho = 0.2; K = 2;

for k=1:num_K

    dummy1 = zeros(N_rand_net,1);
    var_g = NaN(N_rand_net,1); 

    %K = k+1;
    %deg_K(k) = K;
    rho = 0.1 + 0.05*(k-1);
    deg_K(k) = rho;

%----------------------------
%Initialization of network
%----------------------------

for z=1:N_rand_net

    rng(z)
    h = WattsStrogatz(n,K,rho);
    A = adjacency(h);
    A = full(A);
    g_init = linspace(0,1.95,n);  %%rng(5), g_init = 1.96*rand(1,n); 
    for i=1:n
        A(i,i)=1;
    end

    %----------------
    %Initial values
    %----------------
    gbar_init = sum(g_init)/n;
    g0 = g_init';
    p0 = pf + 0.10;
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
    g = NaN(n,T);
    for t=1:T
        gbar = NaN(T,1);

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

     if isreal(g(1:n,t))
        var_g(z) = var(g(1:n,t));
     end

    if isreal(g(1:n,t)) && var_g(z) < 1e-8 
        dummy1(z) = 1;
        break
    end


    end

    gbar(isnan(gbar)) = [];

    if dummy1(z)==1
        consensus(z,k) = gbar(end);
    end


end

    if min(dummy1)==1

        dummy_K(k) = 1;

    end

max_var(k) = max(var_g);

end

max_var_all = max(max_var)
min_dummy = min(dummy_K)

%Plot results
grey = 0;  %black = 0, half-grey = 0.5  grey = 0.75

for z=1:N_rand_net

    figure(1)
    %plot(deg_K,consensus(z,:),'.', 'MarkerSize',6, 'Color',[grey grey grey]), axis([-inf,inf,-inf,inf]), hold on,
    plot(deg_K,consensus(z,:),'.', 'MarkerSize',7.5, 'Color',[grey grey grey]), axis([-inf,inf,-inf,inf]), hold on,

end

%xlabel('Mean degree / 2  (K)'), ylabel('Consensus')
xlabel('Rewiring probability \rho'), ylabel('Consensus')
yline(1+r,'--k')

