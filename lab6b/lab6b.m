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

% Loop filter
[ lf_b,lf_a ] = LF(order,zeta,BT,N,k0,kp);
lf_zi = zeros(1,max(length(lf_a),length(lf_b))-1);

% DDS
[ dds_b,dds_a ] = DDS(k0);
dds_zi = zeros(1,max(length(dds_a),length(dds_b))-1);

% Design the detector
load('bpskcrdedata.mat');
r_t  = bpskcrdedata(2,:);
to = bpskcrdedata(1,:);

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

% Do the downsampling
xk = downsample(x_t,N);
yk = downsample(y_t,N);

% Now's time for the PLL fun...