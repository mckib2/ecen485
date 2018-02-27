%% Lab 5 - PLL
% Nicholas McKibben
% ECEn 485
% 2018-02-21

clear;
close all;

%% First order PLL
% Params
order = 2;
w0 = 2*pi*.0105;
w0dds = 2*pi*.01;
theta = pi;
zeta = 1/sqrt(2);
BT = .002;
N = 2500;
k0 = 1; kp = 1;
plot_speed = N - 1; % controls how often plot is refreshed

% design a loop filter
[ lf_b,lf_a ] = LF(order,zeta,BT,1,k0,kp);
lf_zi = zeros(1,max(length(lf_a),length(lf_b))-1);

% design dds filter
[ dds_b,dds_a ] = DDS(k0);
dds_zi = zeros(1,max(length(dds_a),length(dds_b))-1);

% create our input signal
in = exp(1j*(w0*(1:N) + theta));

% start looping through each sample
out = zeros(1,N); e = out;
for n = 2:N-1
    e(n) = angle(in(n)*conj(out(n))); % Find the phase error
    [ v,lf_zi ] = filter(lf_b,lf_a,e(n),lf_zi); % Run through the loop filter
    [ theta_hat,dds_zi ] = filter(dds_b,dds_a,v + w0dds,dds_zi); % Run through the DDS
    out(n+1) = exp(1j*(theta_hat)); % construct our output
    
    % Show us what's going on
    if ~mod(n,plot_speed)
        show(n,N,e,in,out,theta_hat);
    end
end