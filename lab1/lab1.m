%% Lab 1 - Binary-PAM
% Nicholas McKibben
% ECEn 485
% 2018-02-03

clear;
close all;

if ~exist('lab0','dir')
    fprintf('Adding lab0 to path...\n');
    addpath('../lab0','-end');
end
if ~exist('pam','dir')
    fprintf('Adding PAM to path...\n');
    addpath('../pam','-end');
end

%% Test with binary-PAM
sig = [ 1 0 0 1 ];
b = .25*ones(1,16);
N = 16;
M = 2;
keys = { 0,1};
vals = {-1,1};
LUT2 = LUT(keys,vals);

r = modulator(sig,LUT2,b,N);
[ ~,~,x ] = demodulator(r,b,LUT2,N);

figure(1);
plot(r); hold on; plot(x);
title('demodulator input and the matched filter output');

load('bb2data.mat');
r = bb2data(2,:);
[ s,~,x ] = demodulator(r,b,LUT2,N);

figure(2);
plot(r); hold on; plot(x);

message2 = m2ascii(s,M);
disp(message2);