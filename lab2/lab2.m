%% Lab 2 - 8-ary PAM
% Nicholas McKibben
% ECEn 485
% 2018-02-03

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

% Test out the detector design
sig = [ 0 2 5 6 ];
N = 16;
b = ones(1,N)/sqrt(N);
M = 8;
E = 189;
A = amp(E,M);
keys = {  0,   1,   3,   2, 6, 7,  5,  4 };
vals = { -7*A,-5*A,-3*A,-A, A, 3*A,5*A,7*A };
LUT8 = LUT(keys,vals);

r = modulator(sig,LUT8,b,N);
[ ~,~,x ] = demodulator(r,b,LUT8,N);
figure(3);
plot(r); hold on; plot(x);
title('demodulator input and the matched filter output');

% Now with the example data to get the message
load('bb8data.mat');
r = bb8data(2,:);
[ s,~,x ] = demodulator(r,b,LUT8,N);

figure(4);
plot(r); hold on; plot(x);

message8 = m2ascii(s,M);
disp(message8);
