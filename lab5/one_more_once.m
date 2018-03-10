%% Lab 5 - PLL
% Nicholas McKibben
% ECEn 485
% 2018-03-10

clear;
close all;

% Get the loop filter constants
order = 2;
w0 = 2*pi*.0105;
w0dds = 2*pi*.01;
theta = pi;
zeta = 1/sqrt(2);
BT = .002;
N = 2500;
k0 = 1; kp = 1;
[ ~,~,K1,K2 ] = LF(order,zeta,BT,1,k0,kp);

in = exp(1j*((1:N)*w0 + theta));
out = zeros(size(in));

%% Let's try the state space formulation way
% Initial conditions
sep = 0;
sip = 0;
sl = 1;

for ii = 1:N
    sa = in(ii)*sl;
    sb = angle(sa);
    sc = sb*K1;
    sd = sb*K2;
    sf = sep; % previous se
    se = sd + sf;
    sg = sc + se;
    sh = sg*k0;
    sj = sip; % previous si
    si = sh + w0dds + sj;
    sk = exp(1j*sj);
    sl = conj(sk);
    
    % Update previous vals
    sep = se;
    sip = si;
    
    % Build the output
    out(ii) = sk;
end

%% Look at it
figure(1);
plot(1:N,real(in),'r',1:N,real(out),'b');