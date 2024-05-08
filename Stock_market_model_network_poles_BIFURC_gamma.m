%Stock market model numerical bifurcation-type diagram

%clear; clc; close all

n_gama = 80;
gama_min = 10; gama_max = 1000;
gama_vec = linspace(gama_min,gama_max,n_gama);
n_sims = 8;
n_samp = 50;
T = 1000;  %no. of periods

samp0 = NaN(n_samp,n_sims); samp1 = samp0; 
x_plot = NaN(n_samp*n_sims,n_gama); percent = NaN(n_gama,1); 


for m = 1:n_gama
    
    gama = gama_vec(m);
    brk = zeros(n_sims,1);

    for v=1:n_sims

        rng(v);
        %Stock_market_model_network_poles_sims_gama_infty_shocks_insert 
        Stock_market_model_network_poles_sims_shocks_insert 
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
subplot(1,2,1), plot(gama_vec,x_plot(:,1:end),'o','MarkerSize',2,'Color',[grey grey grey])
%title('Bifurcation diagram'), xlabel('\gamma'), ylabel('Price deviation'), hold on, axes()
%subplot(2,3,3), plot(sigma_vec(m),samp0(:,m),'o','MarkerSize',2.5), title('Bifurcation diagram (with detail)'), 
xlabel('\gamma'), ylabel('Price deviation'), hold on,
axis([-inf,inf,0,inf])
%subplot(2,3,4), plot(sigma_vec(m),samp1(:,m),'o','MarkerSize',2.5), xlabel('\gamma'), 
%ylabel('Average type'), title('Bifurcation diagram'), hold on,
%subplot(1,2,2), plot(sigma_vec(1:m),samp1(:,1:m),'o','MarkerSize',2,'Color',[grey grey grey]), %title('Bifurcation diagram (with detail)')
%xlabel('\gamma'), ylabel('Average type'), hold on,

