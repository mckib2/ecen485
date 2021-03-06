%% Lab 6b - BPSK Differential Encoding
% Nicholas McKibben
% ECEn 485
% 2018-03-04

clear;
close all;

% This is to get m2ascii function
if ~exist('lab0','dir')
    fprintf('Adding lab0 to path...\n');
    addpath('../lab0','-end');
end
% This is to get PLL functions
if ~exist('lab5','dir')
    fprintf('Adding lab5 to path...\n');
    addpath('../lab5','-end');
end

% Params
M = 2;
N = 8;
E = 1;
w0 = 2*pi*.2;
beta = .5;
span = 12;
SYNC = [ 0 0 0 1 0 1 1 0 ];
DataL = 280;
repeats = 4; % The packet is repeated to allow your PLL to acquire and lock
A = sqrt(E);
R = 1;
Fs = R*N;

% Grab the filter coefficients
b = rcosdesign(beta,span,N);

% PLL
order = 2;
k0 = 1; kp = 1;
BT = .01;
zeta = 1/sqrt(2);
w0dds = 0;

% Loop filter
[ ~,~,K1,K2 ] = LF(order,zeta,BT,1,k0,kp);

% Design the detector
load('bpskcrdedata.mat');
r_t  = bpskcrdedata(2,:);
to = bpskcrdedata(1,:);

% Make a local oscillator with which to mix
LOx = @(x) sqrt(2)*cos(w0*x);
LOy = @(x) -1*sqrt(2)*sin(w0*x);

% Now do the mixin'
Ir_t = r_t.*LOx(to);
Qr_t = r_t.*LOy(to);

% Now do the filterin'
x_t = conv(b,Ir_t);
y_t = conv(b,Qr_t);

% Downsamplin'
xk = x_t(1:N:end);
yk = y_t(1:N:end);

% Make some input
in = xk + 1j*yk;
out = zeros(size(in));
a = zeros(size(in));
e = zeros(size(in));

%% Let's try the state space formulation way
% Initial conditions
sep = 0;
sip = 0;
sl = 1;

for ii = 1:numel(in)
    % rotation
    sa = in(ii)*sl;
    
    % decision
    a(ii) = sign(real(sa));
    
    % mix before loop filter
    sb = a(ii)*imag(sa);
    
    % save the error
    e(ii) = sb;
    
    sc = sb*K1;
    sd = sb*K2;
    sf = sep; % previous se
    se = sd + sf;
    sg = sc + se;
    sh = sg*k0;
    sj = sip; % previous si
    si = sh + w0dds + sj;
    sk = cos(sj) + 1j*sin(sj);
    sl = conj(sk);
    
    % Update previous vals
    sep = se;
    sip = si;
    
    % Build the output
    out(ii) = sk;
end

%% Look at it
figure(1);
plot(e);

%% Differential Encoding
d_hat = zeros(size(a));
d_hat(a > 0) = 1;

for ii = 2:numel(d_hat)
    if d_hat(ii) == 0 && d_hat(ii-1) == 0
        b(ii) = 1;
    elseif d_hat(ii) == 0 && d_hat(ii-1) == 1
        b(ii) = 0;
    elseif d_hat(ii) == 1 && d_hat(ii-1) == 0
        b(ii) = 0;
    elseif d_hat(ii) == 1 && d_hat(ii-1) == 1
        b(ii) = 1;
    end
end

%% Check for the SYNC
m = b;
m(m < 0) = 0;
mstr = strjoin(string(m),'');
idx = strfind(mstr,strjoin(string(SYNC),''));

if isempty(idx)
    fprintf('Checking for pi phase..\n');
    m = -a;
    m(m < 0) = 0;
    mstr = strjoin(string(m),'');
    idx = strfind(mstr,strjoin(string(SYNC),''));
end

% try each possibility
for ii = 1:numel(idx)
    start = idx(end-ii+1) + numel(SYNC);
    o = 0;
    try
        m2 = m(start+o:start+DataL-1+o);
        fprintf('Message %d: %s\n',ii,m2ascii(m2,M));
    catch
    end
end