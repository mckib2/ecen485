%% Lab 2 - 8-ary PAM
% Nicholas McKibben
% ECEn 485
% 2018-02-03

clear;
close all;

%% Test with binary-PAM
sig = [ 1 0 0 1 ];
b = .25*ones(1,16);
N = 16;
M = 2;

r = modulator(sig,@binary_lut,b,N);
[ ~,~,x ] = demodulator(r,b,M,N);

figure(1);
plot(r); hold on; plot(x);
title('demodulator input and the matched filter output');

load('bb2data.mat');
r = bb2data(2,:);
[ s,~,x ] = demodulator(r,b,M,N);

figure(2);
plot(r); hold on; plot(x);

addpath('../lab0/','-end');
m2ascii(s,M);

%% Now try out 8-ary PAM
sig = [ 1 0 0 1 ];
b = .25*ones(1,16);
N = 16;
M = 2;

%% 8-ary look up table
function [ out ] = eight_lut(in,~)
    if nargin < 2
        % symbol to A
    else
        % A to symbol
    end
end

%% Binary look up table
function [ out ] = binary_lut(in,~)
    if nargin < 2
        % symbol to A
        out = in;
        out(out == 0) = -1;
    else
        % A to symbol
        out = double(in >= 0);
    end
end

%% Modulator
function [ r ] = modulator(s,lut,b,N)
    x = upsample(lut(s),N);
    r = filter(b,1,x);
end

%% Demodulator
function [ s,a,x ] = demodulator(r,b,M,N)

    % Matched filter
    x = filter(b,1,[0 r]);

    % Downsample
    xk = downsample(x,N);
    xk = xk(2:end);
    
    % Make a decision
    d = -(M-1):2:M;
    for ii = 1:numel(xk)
        test_sig = d - xk(ii);
        [ ~,k ] = min(abs(test_sig));
        a(ii) = d(k);
        s(ii) = binary_lut(a(ii),1);
    end
end