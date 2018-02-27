%% Lab 6a - BPSK Unique Word
% Nicholas McKibben
% ECEn 485
% 2018-02-27

clear;
close all;

% This is to get m2ascii function
if ~exist('lab0','dir')
    fprintf('Adding lab0 to path...\n');
    addpath('../lab0','-end');
end
% This is to get PLL function
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
uw = [ 0 0 0 1 0 1 1 0 ];
DataL = 224;
repeats = 4; % The packet is repeated 4 times to allow your PLL to acquire and lock.
A = E;
R = 1;
Fs = R*N;

% Grab the filter coefficients
b = rcosdesign(beta,span,N);

% PLL
order = 2;
k0 = 1; kp = 1;
BT = .01;
zeta = 1/sqrt(2);

% Loop filter
[ lf_b,lf_a ] = LF(order,zeta,BT,N,k0,kp);
lf_zi = zeros(1,max(length(lf_a),length(lf_b))-1);

% DDS
[ dds_b,dds_a ] = DDS(k0);
dds_zi = zeros(1,max(length(dds_a),length(dds_b))-1);

% Design the detector
load('bpskcruwdata.mat');
r_t  = bpskcruwdata(2,:);
to = bpskcruwdata(1,:);

% Make a local oscillator with which to mix
LOx = @(x) sqrt(2)*cos(w0*x);
LOy = @(x) -1*sqrt(2)*sin(w0*x);

% Now do the mixin'
Ir_t = r_t.*LOx(to);
Qr_t = r_t.*LOx(to);

% Now do the filterin'
fdelay = span/(2*R);
x_t = filter(fliplr(b),1,[ Ir_t zeros(1,fdelay*Fs) ]);
y_t = filter(fliplr(b),1,[ Qr_t zeros(1,fdelay*Fs) ]);

% Compensate for group delay
x_t = x_t(fdelay*Fs+1:end);
y_t = y_t(fdelay*Fs+1:end);

% Make sure everything's looking hunky dory
% eyediagram(x_t,N,1);
% title('x(t)');
% eyediagram(y_t,N,1);
% title('y(t)');

% Do the downsampling
xk = downsample(x_t,N);
yk = downsample(y_t,N);

% Feed into PLL
R = @(v,t) [ cos(t) -sin(t); sin(t) cos(t) ]*v;
e = zeros(1,numel(xk));
out = e; a = e; in = e;
theta_hat = 0;
for k = 1:(numel(xk)-1)
    % Make the decision and mix
%     a(k) = sign(real(x_t(k)*conj(out(k))));
%     in(k) = a(k)*yk(k)*conj(out(k));
    v = R([ xk(k) yk(k) ].',-pi/4);
%     v = sqrt(xk(k)^2+yk(k)^2)*exp(atan(yk(k)/xk(k)))*conj(out(k));
    figure(1);
    hold on;
    plot(v(1),v(2),'o');
    xlim([ -1.5 1.5 ]);
    ylim([ -1.5 1.5 ]);
    drawnow;
    
    a(k) = sign(v(1));
    in(k) = a(k)*v(2);
    
%     e(k) = angle(in(k)); % Find the phase error
    e(k) = in(k);
    [ v,lf_zi ] = filter(lf_b,lf_a,e(k),lf_zi); % Run through the loop filter
    [ theta_hat,dds_zi ] = filter(dds_b,dds_a,v + w0,dds_zi); % Run through the DDS
    out(k+1) = exp(1j*theta_hat); % construct our output
end

% For fun, let's look at the PLL plots
figure(2);
show(numel(xk)-1,numel(xk)-1,e,in,out,theta_hat);

%% Check for the unique word
m = -a;
m(m < 0) = 0;
mstr = strjoin(string(m),'');
idx = strfind(mstr,strjoin(string(uw),''));
start = idx(end) + numel(uw);
o = 0;
m2 = m(start+o:start+DataL-1+o);
disp(m2ascii(m2,M));


% % Grab the message
% message = m2ascii(m(end-DataL+1:end),M);
% disp(message);

% test = sign(xk);
% test(test < 0) = 0;
% a = 14;
% message = m2ascii(test(end-DataL+1-a:end-a),M);
% disp(message);