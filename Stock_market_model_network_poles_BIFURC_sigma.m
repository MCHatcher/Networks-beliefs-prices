%Stock market model numerical bifurcation-type diagram

%clear; clc; close all
n_sigma = 100;
sigma_min = 1e-8; sigma_max = 0.4;
sigma_vec = linspace(sigma_min,sigma_max,n_sigma);
n_sims = 30;
n_samp = 50;
T = 2100;  %no. of periods

samp0 = NaN(n_samp,n_sims); samp1 = samp0; 
x_plot = NaN(n_samp*n_sims,n_sigma); percent = NaN(n_sigma,1); 


for m = 1:n_sigma
    
    sigma_d = sigma_vec(m);
    brk = zeros(n_sims,1);

    %Truncated normal distribution
    pd = makedist('Normal','mu',0,'sigma',sigma_d);
    dbar = 0.5;
    pd_t = truncate(pd,-dbar,dbar);

    for v=1:n_sims

        rng(500+v);
        %Stock_market_model_network_poles_sims_shocks_insert 
        Stock_market_model_network_poles_sims_gama_infty_shocks_insert
        samp0(:,v) = dev(end+1-n_samp:end);
        samp1(:,v) = gbar(end+1-n_samp:end);

        %Check for boundedness
        r1 = 1 - isreal(dev(end)); r2 = isnan(dev(end));  r3 = isinf(dev(end));

        if r1 + r2 + r3 > 0
               brk(v) = 1;
        end

    end

    %Sims with no attractor
    percent(m) = 100*sum(brk)/n_sims;

    %Prepare to plot
    x_plot(:,m) = reshape(samp0,1,[]);

end 

grey = 0.1;
figure(1)
subplot(1,2,2), plot(sigma_vec,x_plot(:,1:end),'o','MarkerSize',2,'Color',[grey grey grey])
%title('Bifurcation diagram'), xlabel('\gamma'), ylabel('Price deviation'), hold on, axes()
%subplot(2,3,3), plot(sigma_vec(m),samp0(:,m),'o','MarkerSize',2.5), title('Bifurcation diagram (with detail)'), 
xlabel('\sigma_d'), ylabel('Price deviation'), hold on,
axis([-inf,inf,-inf,inf]), title('\gamma \rightarrow \infty')
%subplot(2,3,4), plot(sigma_vec(m),samp1(:,m),'o','MarkerSize',2.5), xlabel('\gamma'), 
%ylabel('Average type'), title('Bifurcation diagram'), hold on,
%subplot(1,2,2), plot(sigma_vec(1:m),samp1(:,1:m),'o','MarkerSize',2,'Color',[grey grey grey]), %title('Bifurcation diagram (with detail)')
%xlabel('\gamma'), ylabel('Average type'), hold on,

