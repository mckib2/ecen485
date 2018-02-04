%% Lab 2 - 8-ary PAM
% Nicholas McKibben
% ECEn 485
% 2018-02-03

clear;
close all;

%% Test with binary-PAM
% sig = [ 1 0 0 1 ];
% b = .25*ones(1,16);
% N = 16;
% M = 2;
% LUT2 = containers.Map({0,1},{-1,1});
% rLUT2 = containers.Map({-1,1},{0,1});
% 
% r = modulator(sig,LUT2,b,N);
% [ ~,~,x ] = demodulator(r,b,rLUT2,M,N);
% 
% figure(1);
% plot(r); hold on; plot(x);
% title('demodulator input and the matched filter output');
% 
% load('bb2data.mat');
% r = bb2data(2,:);
% [ s,~,x ] = demodulator(r,b,rLUT2,M,N);
% 
% figure(2);
% plot(r); hold on; plot(x);
% 
% addpath('../lab0/','-end');
% message2 = m2ascii(s,M);
% disp(message2);

%% Now try out 8-ary PAM
% Test out the detector design
sig = [ 0 2 5 6 ];
N = 16;
b = ones(1,N)/sqrt(N);
M = 8;
E = 189;
A = amp(E,M);
keys = {  0, 1, 3, 2,6,7,5,4 };
vals = { -7*A,-5*A,-3*A,-A,A,3*A,5*A,7*A };
LUT8 = containers.Map(keys,vals);
rLUT8 = containers.Map(vals,keys);

r = modulator(sig,LUT8,b,M,N);
[ ~,~,x ] = demodulator(r,b,rLUT8,M,N);

figure(3);
plot(r); hold on; plot(x);
title('demodulator input and the matched filter output');

%% 
load('bb8data.mat');
r = bb8data(2,:);
[ s,~,x ] = demodulator(r,b,rLUT8,M,N);

figure(4);
plot(r); hold on; plot(x);

message8 = m2ascii(s,M);
disp(message8);

function [ A ] = amp(E,M)
    A = sqrt(3*E/(M^2-1));
end

%% Modulator
function [ r ] = modulator(s,lut,b,M,N)
    keys = cell(size(s));
    for ii = 1:numel(s)
        keys{ii} = s(ii);
    end
    x = upsample(cell2mat(values(lut,keys)),N);
    r = filter(b,1,x);
end

%% Demodulator
function [ s,a,x ] = demodulator(r,b,rlut,M,N)

    % Matched filter
    x = filter(b,1,[0 r]);

    % Downsample
    xk = downsample(x,N);
    xk = xk(2:end);
    
    % Make a decision
    d = cell2mat(keys(rlut));
    a = zeros(1,numel(xk));
    s = a;
    for ii = 1:numel(xk)
        test_sig = d - xk(ii);
        [ ~,k ] = min(abs(test_sig));
        a(ii) = d(k);
        s(ii) = rlut(a(ii));
    end
end