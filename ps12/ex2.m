clc
close all
clear all

%Parameters
sigma = 0.2; %volatility
s0 = 100; %initial stock price
r = 0.01; %risk-free rate
T = 1; %maturity
K = [70, 90, 110, 130]; %strike price
m = [2,3,4,5];

%Initialization
%price_bs = zeros(1,4);
price_bs = zeros(4,4);
for k = 1:4
    price_bs(k,:) = BS(sigma, s0, r, K(k));
end
%price_bs(1,:) = BS(sigma, s0, r, K(1));
p_exact_sim = zeros(4, 4);
p_euler_sch = zeros(4, 4);
p_sim_euler_sch = zeros(4, 4);
p_sec_ord_sch = zeros(4, 4);
p_exact_sim_cv = zeros(4, 4);
p_euler_sch_cv = zeros(4, 4, 4);

%Computation
for i = 1:4
    M = 1000*(2^(2*m(i))); %number of simulations
    N = 10*(2^m(i)); %number of time steps
    dt = 1/N; %time step
    
    %price_bs(i) = BS(M, N, dt, sigma, s0, r, K(1));
    [p_exact_sim(i,1), p_exact_sim(i,2), p_exact_sim(i,3), p_exact_sim(i,4)] = exact_sim(M, N, dt, sigma, s0, r, K(1));
    [p_euler_sch(i,1), p_euler_sch(i,2), p_euler_sch(i,3), p_euler_sch(i,4)]  = euler_sch(M, N, dt, sigma, s0, r, K(1));
    [p_sim_euler_sch(i,1), p_sim_euler_sch(i,2), p_sim_euler_sch(i,3), p_sim_euler_sch(i,4)] = simp_euler_sch(M, N, dt, sigma, s0, r, K(1));
    [p_sec_ord_sch(i,1), p_sec_ord_sch(i,2),p_sec_ord_sch(i,3),p_sec_ord_sch(i,4)] = sec_ord_sch(M, N, dt, sigma, s0, r, K(1));
    [p_exact_sim_cv(i,1), p_exact_sim_cv(i,2),p_exact_sim_cv(i,3), p_exact_sim_cv(i,4)] = exact_sim_cv(M, N, dt, sigma, s0, r, K(1));
    for k = 1:4
        [p_euler_sch_cv(i,k,1), p_euler_sch_cv(i,k,2), p_euler_sch_cv(i,k,3), p_euler_sch_cv(i,k,4)] = euler_sch_cv(M, N, dt, sigma, s0, r, K(k));
    end
end

%Plots
figure()
t = tiledlayout(2,3);
title(t, 'Results for parts i) to vi)')
% i) Exact simulation
nexttile
plot(m,price_bs(1,:),'r',m,p_exact_sim(:,1),'g',m,p_exact_sim(:,3),'b',m,p_exact_sim(:,4),'b')
xlabel('m')
ylabel('Option Price')
legend('BS price', 'Estimation', '95% Confidence interval')
title('i) Exact simulation')
% ii) Euler scheme
nexttile
plot(m,price_bs(1,:),'r',m,p_euler_sch(:,1),'g',m,p_euler_sch(:,3),'b',m,p_euler_sch(:,4),'b')
xlabel('m')
ylabel('Option Price')
legend('BS price', 'Estimation', '95% Confidence interval')
title('ii) Euler scheme')
% iii) Simplified Euler scheme
nexttile
plot(m,price_bs(1,:),'r',m,p_sim_euler_sch(:,1),'g',m,p_sim_euler_sch(:,3),'b',m,p_sim_euler_sch(:,4),'b')
xlabel('m')
ylabel('Option Price')
legend('BS price', 'Estimation', '95% Confidence interval')
title('iii) Simplified Euler scheme')
% iv) Second order scheme
nexttile
plot(m,price_bs(1,:),'r',m,p_sec_ord_sch(:,1),'g',m,p_sec_ord_sch(:,3),'b',m,p_sec_ord_sch(:,4),'b')
xlabel('m')
ylabel('Option Price')
legend('BS price', 'Estimation', '95% Confidence interval')
title('iv) Second order scheme')
% v) Exact simulation with control variate
nexttile
plot(m,price_bs(1,:),'r',m,p_exact_sim_cv(:,1),'g',m,p_exact_sim_cv(:,3),'b',m,p_exact_sim_cv(:,4),'b')
xlabel('m')
ylabel('Option Price')
legend('BS price', 'Estimation', '95% Confidence interval')
title('v) Exact simulation using control variate estimator')
% vi) Euler scheme with control variate
nexttile
plot(m,price_bs(1,:),'r',m,p_euler_sch_cv(:,1,1),'g',m,p_euler_sch_cv(:,1,3),'b',m,p_euler_sch_cv(:,1,4),'b')
xlabel('m')
ylabel('Option Price')
legend('BS price', 'Estimation', '95% Confidence interval')
title('vi) Euler scheme using control variate estimator')

%Euler scheme with control variate for all K's
figure()
t = tiledlayout(2,2);
title(t, 'Results for vii), Euler scheme using control variate estimator for K = 70, 90, 110, 130')
for k = 1:4
    nexttile
    plot(m,price_bs(k,:),'r',m,p_euler_sch_cv(:,k,1),'g',m,p_euler_sch_cv(:,k,3),'b',m,p_euler_sch_cv(:,k,4),'b')
    xlabel('m')
    ylabel('Option Price')
    legend('BS price', 'Estimation', '95% Confidence interval')
    title(sprintf('K = %d',K(k)))
end

%FUNCTIONS
% BS formula put option
function [price,stdev,conflow,confup] = BS(sigma, s0, r, K)
% Computes price of a European option
% M = number of simulations
% N = number of time steps

    % BS formula put option
    d1 = (log(s0/K)+r+sigma^2/2)/sigma;
    d2 = (log(s0/K)+r-sigma^2/2)/sigma;
    price = K*exp(-r)*normcdf(-d2)-s0*normcdf(-d1); %Put option
    %price = s0*normcdf(d1)-K*exp(-r)*normcdf(d2); -> Call option
    stdev = 0; % dummy
    conflow=price; % dummy
    confup=price; % dummy
end

% i) Exact simulation
function [price,stdev,conflow,confup] = exact_sim(M, N, dt, sigma, s0, r, K)
    Z = randn(1,M); % generate normal rv
    S = s0*exp((r-sigma^2/2)+sigma*Z);
        % final stock price
    [price,stdev,conflow,confup] = emp(S,r,K,M);
        % estimated option price
end

% ii) Euler scheme
function [price,stdev,conflow,confup] = euler_sch(M, N, dt, sigma, s0, r, K)
    S=s0*ones(1,M); % initialize stock price
    for i=1:N
        Z = randn(1,M); % generate normal rv
        S = S + r*S*dt + sigma*sqrt(dt)*S.*Z; % Euler recursion
    end
    [price,stdev,conflow,confup] = emp(S,r,K,M);
end

% iii) Simplified Euler scheme
function [price,stdev,conflow,confup] = simp_euler_sch(M, N, dt, sigma, s0, r, K)
    S=s0*ones(1,M); % initialize stock price
    for i=1:N
        Z = 2*randi([0,1],1,M)-1; % generate Bernoulli rv (randint is an obsolete function)
        S = S + r*S*dt + sigma*sqrt(dt)*S.*Z; % Euler recursion
    end
    [price,stdev,conflow,confup] = emp(S,r,K,M);
end

% iv) Second order scheme
function [price,stdev,conflow,confup] = sec_ord_sch(M, N, dt, sigma, s0, r, K)
    S=s0*ones(1,M); % initialize stock price
    for i=1:N
        Z = randn(1,M); % generate normal rv
        S = S + r*S*dt + sigma*sqrt(dt)*S.*Z + r*sigma*(dt)^(3/2)*S.*Z+ 1/2*sigma^2*dt*S.*(Z.^2-1)+ 1/2*r^2*dt^2*S;
                % second order recursion
    end
    [price,stdev,conflow,confup] = emp(S,r,K,M);
end

% v) Exact simulation and use e^{-rT}S_T as a control variable
function [price,stdev,conflow,confup] = exact_sim_cv(M, N, dt, sigma, s0, r, K)
    Z = randn(1,M); % generate normal rv
    S = s0*exp((r-sigma^2/2)+sigma*Z);
        % final stock price
    Y = exp(-r)*S;
    X = exp(-r)*F(S,K);
    me_X = mean(exp(-r)*F(S,K));
    b = (X-me_X).*(Y-s0)/((Y-s0).^2);
    X_b = X - b*(Y-s0);
    
    price = mean(X_b);
    stdev = std(X_b);
    conflow = price-1.96*stdev/sqrt(M);
    confup = price+1.96*stdev/sqrt(M);
    %[price,stdev,conflow,confup] = emp(S, r);
end

% vi) Euler scheme and use e^{-rT}S_T as a control variable
function [price,stdev,conflow,confup] = euler_sch_cv(M, N, dt, sigma, s0, r, K)
    S=s0*ones(1,M); % initialize stock price
    for i=1:N
        Z = randn(1,M); % generate normal rv
        S = S + r*S*dt + sigma*sqrt(dt)*S.*Z; % Euler recursion
    end
    Y = exp(-r)*S;
    X = exp(-r)*F(S,K);
    me_X = mean(exp(-r)*F(S,K));
    b = (X-me_X).*(Y-s0)/((Y-s0).^2);
    X_b = X - b*(Y-s0);
    
    price = mean(X_b);
    stdev = std(X_b);
    conflow = price-1.96*stdev/sqrt(M);
    confup = price+1.96*stdev/sqrt(M);
    %[price,stdev,conflow,confup] = emp(S, r);
end

% Auxiliary subfunction computing empirical mean, standard deviation, and
% 95% confidence interval of sample S
function [me,st,colow,coup] = emp(S, r, K, M)
    me = mean(exp(-r)*F(S,K));
        % estimated option price
    st = std(exp(-r)*F(S,K));
        % empirical standard deviation
    colow = me-1.96*st/sqrt(M);
        % 95% confidence interval lower bound
    coup = me+1.96*st/sqrt(M);
        % 95% confidence interval upper bound
end

%Payoff of a put option
function F = F(x,K)
    F = max(0,K-x);
    % F = x;
end