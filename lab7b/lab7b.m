%% Lab 7b - Carrier Phase Synchronization for QPSK, Differential Encoding
% Nicholas McKibben
% ECEN 485
% 2018-03-10

clear;
close all;

% This is to get PLL functions
if ~exist('lab5','dir')
    fprintf('Adding lab5 to path...\n');
    addpath('../lab5','-end');
end
% This is to get pam helper functions
if ~exist('pam','dir')
    fprintf('Adding PAM to path...\n');
    addpath('../pam','-end');
end
% This is to get m2ascii function
if ~exist('lab0','dir')
    fprintf('Adding lab0 to path...\n');
    addpath('../lab0','-end');
end

%% Params
N = 8;
M = 4;
w0 = 2*pi*.3;
E = 2;
A = sqrt((3*E)/(2*(M - 1)));

% Matched filter design
beta = .5;
span = 12;
b = rcosdesign(beta,span,N);

% PLL
order = 2;
k0 = 1; kp = 1;
BT = .01;
zeta = 1/sqrt(2);
w0dds = 0;

% Loop filter
[ ~,~,K1,K2 ] = LF(order,zeta,BT,1,k0,kp);

%% Detector time
load('qpskcrdedata.mat');
n = qpskcrdedata(1,:);
r_t = qpskcrdedata(2,:);

% Make a local oscillator with which to mix
LOx = @(x) sqrt(2)*cos(w0*x);
LOy = @(x) -1*sqrt(2)*sin(w0*x);

% Now do the mixin'
Ir_t = r_t.*LOx(n);
Qr_t = r_t.*LOy(n);

% Now do the filterin'
x_t = conv(b,Ir_t);
y_t = conv(b,Qr_t);

% Downsamplin'
xk = x_t(1:N:end);
yk = y_t(1:N:end);

% Make some input
in = xk + 1j*yk;
out = zeros(size(in));
a0 = zeros(size(in));
a1 = zeros(size(in));
e = zeros(size(in));

%% Let's try the state space formulation way
% Initial conditions
sr = 1;
skp = 0;
sop = 0;

for ii = 1:numel(in)
    % rotation
    sa = in(ii)*sr;
    
    % decisions
    a0(ii) = A*sign(real(sa));
    a1(ii) = A*sign(imag(sa));
    
    sf = a0(ii)*imag(sa);
    sg = a1(ii)*real(sa);
    sh = sf - sg;
    
    % save the error
    e(ii) = sh;
    
    % loop filter
    si = sh*K1;
    sj = sh*K2;
    sl = skp; % previous sk
    sk = sj + sl;
    sm = si + sk;
    
    % DDS
    sn = sm*k0;
    sp = sop; % previous so
    so = sn + w0dds + sp;
    sq = cos(sp) + 1j*sin(sp);
    sr = conj(sq);
    
    % Update previous vals
    skp = sk;
    sop = so;
    
    % Build the output
    out(ii) = sq;
end

%% Look at it
figure(1);
plot(e);
title('Error');
xlabel('sample');
ylabel('e(n)');

%% Differential Encoding
a0 = a0; % Choose sign to get the right answer
a1 = -a1;  % Choose sign to get the right answer

d_hat0 = zeros(size(a0));
d_hat1 = zeros(size(a1));
d_hat0(a0 > 0) = 1;
d_hat1(a1 > 0) = 1;

d_hat = [ d_hat0; d_hat1 ];
d_hat = char(strjoin(string(d_hat(:).'),''));

% assume 00 at beginning
d_hat = strcat('00',d_hat);

idx = 0;
for ii = (1+log2(M)):log2(M):numel(d_hat)
    idx = idx + 1;
    cur = strcat(d_hat(ii),d_hat(ii+1));
    prev =  strcat(d_hat(ii-2),d_hat(ii-1));
    
    if strcmp(prev,'00') && strcmp(cur,'00')
        b(idx) = 0;
    elseif strcmp(prev,'00') && strcmp(cur,'01')
        b(idx) = 1;
    elseif strcmp(prev,'00') && strcmp(cur,'10')
        b(idx) = 2;
    elseif strcmp(prev,'00') && strcmp(cur,'11')
        b(idx) = 3;
        
    elseif strcmp(prev,'01') && strcmp(cur,'00')
        b(idx) = 2;
    elseif strcmp(prev,'01') && strcmp(cur,'01')
        b(idx) = 0;
    elseif strcmp(prev,'01') && strcmp(cur,'10')
        b(idx) = 3;
    elseif strcmp(prev,'01') && strcmp(cur,'11')
        b(idx) = 1;
        
        
    elseif strcmp(prev,'10') && strcmp(cur,'00')
        b(idx) = 1;
    elseif strcmp(prev,'10') && strcmp(cur,'01')
        b(idx) = 3;
    elseif strcmp(prev,'10') && strcmp(cur,'10')
        b(idx) = 0;
    elseif strcmp(prev,'10') && strcmp(cur,'11')
        b(idx) = 2;
        
    elseif strcmp(prev,'11') && strcmp(cur,'00')
        b(idx) = 3;
    elseif strcmp(prev,'11') && strcmp(cur,'01')
        b(idx) = 2;
    elseif strcmp(prev,'11') && strcmp(cur,'10')
        b(idx) = 1;
    elseif strcmp(prev,'11') && strcmp(cur,'11')
        b(idx) = 0;
    end
end

%% Find Unique Word
uw = [ 0 0 0 0 0 0 0 0 ];
idx = strfind(strjoin(string(b),''),strjoin(string(uw),''));

DataL = 175;

% try each possibility
for ii = 1:numel(idx)
    start = idx(end-ii+1) + numel(uw);
    o = 0;
    try
        s_hat = b(start+o:start+DataL-1+o);
        fprintf('Message %d: %s\n',ii,m2ascii(s_hat,M));
    catch
    end
end

%% Eye diagrams
eyediagram(x_t,N);
title('x(t)');
eyediagram(y_t,N);
title('y(t)');

%% Constellation
scatterplot(a0 + 1j*a1);
title('Constellation');