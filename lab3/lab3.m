%% Lab 3a - Binary Phase Shift Keying
% Nicholas McKibben
% ECEn 485
% 2018-02-05

clear;
close all;

% This is to get m2ascii function
if ~exist('lab0','dir')
    fprintf('Adding lab0 to path...\n');
    addpath('../lab0','-end');
end
% This is to get pam helper functions
if ~exist('pam','dir')
    fprintf('Adding PAM to path...\n');
    addpath('../pam','-end');
end

% s = [ 1 0 0 1 1 1 0 1 ];
N = 8;
M = 2;
E = 1; A = amp(E,M);
beta = 1;
span = 12;
sps = 8;
b = rcosdesign(beta,span,sps);
% b = ones(1,N)/sqrt(N);
Omega0 = 2*pi*1/4;

s = randi([ 0 (M-1) ],1,200);

keys = {  0,1 };
vals = { -A,A };
LUT2 = LUT(keys,vals);
n = 0:(N*length(s) - 1);

r = modulator(s,LUT2,b,N);
[ s_hat,~,x,xk ] = demodulator(r,b,LUT2,N);

s_hat = s_hat(12:end);

% eyediagram(x,N,1,1);
% scatterplot(xk,N);

% Check to make sure everything is hunky dory
if ~isequal(s(1:numel(s_hat)),s_hat)
    fprintf('Something has gone terribly wrong...\n');
end

%% Exercise
load('bpskdata.mat');
r = bpskdata(2,:);
n = bpskdata(1,:);
[ s,~,x,xk ] = demodulator(r.*cos(Omega0*n),b,LUT2,N);

message = m2ascii(s(end-154:end-N),M);
disp(message);

% eyediagram(x,N,1,1);
% scatterplot(xk,N);