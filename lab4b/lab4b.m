%% Lab 4b - M-ary Quadrature Amplitude Modulation (MQAM)
% Nicholas McKibben
% ECEn 485
% 2018-02-18

close all;
clear;

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


% Params given in lab
M = 16;
N = 8;
w0 = 2*pi*1/4;
E = 40;
beta = .3;
span = 16;
DataL = 322;
A = sqrt(3*E/(2*(M-1)));
R = 1;
Fs = R*N;

% Get our filter coefficients
b = rcosdesign(beta,span,N);
% b = b/max(abs(b)); % normalize filter coefficients

% Define our local oscillators
LOx = @(x) sqrt(2)*cos(w0*x);
LOy = @(x) -1*sqrt(2)*sin(w0*x);


%% Detector

load('qam16data.mat');
r_t = qam16data(2,:);
to = qam16data(1,:);

% Mix it up
Ir_t = r_t.*LOx(to);
Qr_t = r_t.*LOy(to);

% Filter - make sure to flush out all the useful samples at the end
fdelay = span/(2*R);
x_t = filter(fliplr(b),1, [ Ir_t zeros(1,fdelay*Fs) ]);
y_t = filter(fliplr(b),1, [ Qr_t zeros(1,fdelay*Fs) ]);

% Compensate for group delay
x_t = x_t(fdelay*Fs+1:end);
y_t = y_t(fdelay*Fs+1:end);

% Make sure everything's looking hunky dory
eyediagram(x_t,N,1);
title('x(t)');
eyediagram(y_t,N,1);
title('y(t)');

% Downsample
xk = downsample(x_t,N);
yk = downsample(y_t,N);

% And finally tune it up so we get the samples we were looking for
a = N+1;
xk = xk(a:DataL+a-1);
yk = yk(a:DataL+a-1);

%% Let's make some decisions...
kx = num2cell(0:M-1);
vx = num2cell(A*[ -3 -3 -3 -3 -1 -1 -1 -1  3  3  3  3  1  1  1  1 ]);
vy = num2cell(A*[  3  1 -3 -1  3  1 -3 -1  3  1 -3 -1  3  1 -3 -1 ]);
LUTx = LUT(kx,vx);
LUTy = LUT(kx,vy);

dx = cell2mat(values(LUTx.forward));
[ ~,idx ] = min((dx.' - xk).^2);
x_hat = dx(idx);

dy = cell2mat(values(LUTy.forward));
[ ~,idy ] = min((dy.' - yk).^2);
y_hat = dy(idy);

% Look at the constellation
figure(3);
plot(x_hat,y_hat,'o');
title('Constellation, 16-ary QAM');

s_hat = zeros(1,DataL);
for ii = 1:numel(x_hat)
    s_hat(ii) = helper(x_hat(ii)/A,y_hat(ii)/A);
end

% Grab the message
m = m2ascii(s_hat,M);
disp(m);

% There's a more elegant way to do this.  I don't have time right not.
function [ out ] = helper(x,y)
    if x == -3
        if y == 3
            out = 0;
        elseif y == 1
            out = 1;
        elseif y == -3
            out = 2;
        elseif y == -1
            out = 3;
        end
    elseif x == -1
        if y == 3
            out = 4;
        elseif y == 1
            out = 5;
        elseif y == -3
            out = 6;
        elseif y == -1
            out = 7;
        end
    elseif x == 3
        if y == 3
            out = 8;
        elseif y == 1
            out = 9;
        elseif y == -3
            out = 10;
        elseif y == -1
            out = 11;
        end
    elseif x == 1
        if y == 3
            out = 12;
        elseif y == 1
            out = 13;
        elseif y == -3
            out = 14;
        elseif y == -1
            out = 15;
        end
    end
end