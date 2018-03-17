%% Lab 7a - Carrier Phase Synchronization for QPSK Using the Unique Word
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
load('qpskcruwdata.mat');
n = qpskcruwdata(1,:);
r_t = qpskcruwdata(2,:);

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

%% Find Unique Word
uw = [ 1   1   1   1   1   0   0   0;
       1   0   0   1   0   1   0   0 ];
m0 = -a0; % adjust sign until you find message
m1 = a1;  % adjust sign until you find message

% convert to bits
m0(m0 < 0) = 0; m1(m1 < 0) = 0;

% let's try to find the uw
uw0str = strjoin(string(uw(1,:)),'');
uw1str = strjoin(string(uw(2,:)),'');
m0str  = strjoin(string(m0),'');
m1str  = strjoin(string(m1),'');

id0 = strfind(m0str,uw0str);
id1 = strfind(m1str,uw1str);

% find matching indicies
idx = intersect(id0,id1);

keys = { '00', '01', '10', '11' };
vals = { 0,1,2,3 };
lut = LUT(keys,vals);

DataL = 84;

% try each possibility
for ii = 1:numel(idx)
    start = idx(end-ii+1) + numel(uw);
    o = -N;
    try
        m00 = m0(start+o:start+DataL-1+o);
        m11 = m1(start+o:start+DataL-1+o);
        m = [ m00; m11 ];
        
        s_hat = zeros(1,size(m,2));
        for nn = 1:size(m,2)
            s_hat(nn) = lut.forward(char(strjoin(string(m(:,nn).'),'')));
        end
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
scatterplot(xk + 1j*yk);
title('Constellation');