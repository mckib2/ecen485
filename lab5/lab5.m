%% Lab 5 - PLL
% Nicholas McKibben
% ECEn 485
% 2018-02-21

clear;
close all;

%% First order PLL
% Params
w0 = 2*pi*.01;
theta = pi/2;
zeta = 100;%1/sqrt(2);
BT = .05;
N = 1000;
k0 = .02;
plot_speed = 15; % controls how often plot is refreshed

% design a loop filter
[ lf_b,lf_a ] = LF(zeta,BT);
lf_zi = zeros(1,max(length(lf_a),length(lf_b))-1);

% design dds filter
[ dds_b,dds_a ] = DDS(k0,w0);
dds_zi = zeros(1,max(length(dds_a),length(dds_b))-1);

% create our input signal
in = exp(1j*(w0*(1:N) + theta));

% start looping through each sample
out = zeros(1,N); e = out;
for n = 2:N
    e(n) = angle(in(n)*conj(out(n-1))); % Find the phase error
    [ v,lf_zi ] = filter(lf_b,lf_a,e(n),lf_zi); % Run through the loop filter
    [ theta_hat,dds_zi ] = filter(dds_b,dds_a,v,dds_zi); % Run through the DDS
    out(n) = exp(1j*(w0*n + theta_hat)); % construct our output
    
    % Show us what's going on
    if ~mod(n,plot_speed)
        show(n,N,e,in,out,theta_hat);
    end
end

function [ b,a ] = DDS(k0,w0)
    % Give filter parameters for the DDS
    b = [ k0 w0 ];
    a = [ 1 -1 ];
end

function [ b,a ] = LF(zeta,BT,N)
    % Generate loop filter coefficients.
    % There are two ways we can call this function:
    %    > BnT
    %    > BnTs
    % If we get an N, then we know we have the BnTs case. If not, then we
    % can set N = 1 and everything degrades gracefully.
    if nargin < 3
        N = 1;
    end
        
    den = 1 + (2*zeta/N)*(BT/(zeta + 1/(4*zeta))) + ...
        (BT/(N*(zeta + 1/(4*zeta))))^2;

    K1_num = (4*zeta/N)*(BT/(zeta + 1/(4*zeta)));
    K1 = K1_num/den;
    
    K2_num = (4/N^2)*(BT/(zeta + 1/(4*zeta)))^2;
    K2 = K2_num/den;
    
    % from example
%     K1 = .1479;
%     K2 = .0059;
    
    b = [ (K1+K2) -K1 ];
    a = [ 1 -1 ];
end

function show(n,N,e,in,out,theta_hat)
    % Plot phase error and signals
    figure(1);
    subplot(2,1,1);
    plot(e(2:n),'k-');
    xlim([ 0 N ]);
    title('Phase Error');
    xlabel(sprintf('$ e(n) $ = %g',e(n)),'Interpreter','latex');
    ylabel('$ \theta_e(n) $','Interpreter','latex');
    
    subplot(2,1,2);
    plot(real(in(2:n)),'k-');
    hold on;
    plot(real(out(2:n)),'k--');
    xlabel([ '$ \hat{\theta} $ = ' num2str(theta_hat) ],'Interpreter','latex');
    xlim([ 0 N ]);
    ylim([ -1 1 ]);
    ylabel('$ Re\{e^{j(\cdot)}\} $','Interpreter','latex');
    hold off;
    drawnow;
end