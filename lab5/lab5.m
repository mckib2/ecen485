%% Lab 5 - PLL
% Nicholas McKibben
% ECEn 485
% 2018-02-21

clear;
close all;

%% First order PLL

w0 = 2*pi*.01;
theta = pi/2;
zeta = 1/sqrt(2);
BT = .05;
N = 500;
plot_speed = 8; % controls how often plot is refreshed

% design a loop filter
[ lf_b,lf_a ] = LF(zeta,BT,N);
lf_zi = zeros(1,max(length(lf_a),length(lf_b))-1);

% design dds filter
[ dds_b,dds_a ] = DDS(1,w0);
dds_zi = zeros(1,max(length(dds_a),length(dds_b))-1);

out = zeros(1,N); in = out; e = out;
for n = 2:N
    in(n) = exp(1j*(w0*n + theta));
    [ theta_hat,e(n),lf_zi,dds_zi ] = pll(in(n)*conj(out(n-1)),lf_b,lf_a,lf_zi,dds_b,dds_a,dds_zi);
    out(n) = exp(1j*(w0*n + theta_hat));
    
    % Show us what's going on
    if ~mod(n,plot_speed)
        show(n,N,e,in,out,theta_hat);
    end
end

function [ theta_hat,e,lf_zf,dds_zf ] = pll(in,lf_b,lf_a,lf_zi,dds_b,dds_a,dds_zi)
    e = angle(in);
    [ v,lf_zf ] = filter(lf_b,lf_a,e,lf_zi);
    [ theta_hat,dds_zf ] = filter(dds_b,dds_a,v,dds_zi);
end

function [ b,a ] = DDS(k0,w0)
    b = [ k0 w0 ];
    a = [ 1 -1 ];
end

function [ b,a ] = LF(zeta,BT,N)
    den = 1 + (2*zeta/N)*(BT/(zeta + 1/(4*zeta))) + ...
        (BT/(N*(zeta + 1/(4*zeta))))^2;

    K1_num = (4*zeta/N)*(BT/(zeta + 1/(4*zeta)));
    K1 = K1_num/den;
    
    K2_num = (4/N^2)*(BT/(zeta + 1/(4*zeta)))^2;
    K2 = K2_num/den;
    
    % from example
    K1 = .1479;
    K2 = .0059;
    
    b = [ (K1+K2) -K1 ];
    a = [ 1 -1 ];
end

function show(n,N,e,in,out,theta_hat)
    figure(1);
    subplot(2,1,1);
    plot(e(1:n),'k-');
    xlim([ 0 N ]);
    title('Phase Error');
    xlabel(sprintf('e(n) = %g',e(n)));
    
    subplot(2,1,2);
    plot(real(in(1:n)),'k-');
    hold on;
    plot(real(out(1:n)),'k--');
    xlabel([ '$ \hat{\theta} $ = ' num2str(theta_hat) ],'Interpreter','latex');
    xlim([ 0 N ]);
    ylim([ -1 1 ])
    hold off;
    drawnow;
end