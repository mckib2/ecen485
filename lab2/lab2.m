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

r = modulater(sig,@binary_lut,b,N);
[ a,x ] = demodulater(r,b,2,16);

figure(1);
plot(r); hold on; plot(x);
title('demodulator input and the matched filter output');




%% 8-ary look up table
function [ out ] = eight_lut(in)
    out = in;
    out
end

%% Binary look up table
function [ out ] = binary_lut(in)
    out = in;
    out(out == 0) = -1;
end

%% Modulator
function [ r ] = modulater(s,lut,b,N)
    x = upsample(lut(s),N);
    r = filter(b,1,x);
end

%% Demodulator
function [ a,x ] = demodulater(r,b,M,N)

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
    end
end