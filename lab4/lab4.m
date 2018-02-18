%% Lab 3b
% Nicholas McKibben
% ECEn 485
% 2018-02-10

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

% Params
R = 1; % data rate - just for funsies
N = 8;
Fs = R*N;
beta = 0.5;
span = 12;
b = rcosdesign(beta,span,N);
b = b/max(abs(b)); % normalize filter coefficients
w0 = 2*pi*1/4;

% Find A for 8PSK
E = 9;
M = 8;
A = sqrt(E);

% Define our local oscillators
LOx = @(x) sqrt(2)*cos(w0*x);
LOy = @(x) -1*sqrt(2)*sin(w0*x);

%% Test signal time
kx = { 0,1,         2,        3, 4,         5, 6, 7 };
vx = { A,A/sqrt(2),-A/sqrt(2),0, A/sqrt(2), 0,-A,-A/sqrt(2) };
vy = { 0,A/sqrt(2), A/sqrt(2),A,-A/sqrt(2),-A, 0,-A/sqrt(2) };
LUTx = LUT(kx,vx);
LUTy = LUT(kx,vy);

% sig = [ 3 1 4 1 5 ];
sig = randi([ 0 7 ],[ 1 200 ]);

i_t = cell2mat(values(LUTx.forward,num2cell(sig)));
q_t = cell2mat(values(LUTy.forward,num2cell(sig)));

% i_t = [ A A -A -A ]/sqrt(2);
% q_t = [ A -A A -A ]/sqrt(2);

% For the polyphase filter, append span/2 zeros at the end to flush all
% the useful samples out of the filter
i_t = [ i_t zeros(1,span/2) ];
q_t = [ q_t zeros(1,span/2) ];

% Now do the upsampling
i_t = upsample(i_t,N);
q_t = upsample(q_t,N);

% Now we can filter
I_t = filter(b,1,i_t);
Q_t = filter(b,1,q_t);

% Compensate for the raised cosine filter group delay by delaying the
% input signal
fdelay = span/(2*R);
I_t = I_t(fdelay*Fs+1:end);
Q_t = Q_t(fdelay*Fs+1:end);

% For our test signal, it was 4 long - mix in the LO
to = 1000*(0:numel(sig)*N - 1)/Fs;
I_t = I_t.*LOx(to);
Q_t = Q_t.*LOy(to);

% Now add 'em together and you've got a modulated signal
s_t = I_t + Q_t;

%% Now for the detector part
% r_t = s_t; % test signal
% DataL = numel(sig);

% Load in the no-kidding real data
load('psk8data.mat');
r_t = psk8data(2,:);
DataL = 350; % the message is this long
to = psk8data(1,:); % new time vector

% Mix it up
Ir_t = r_t.*LOx(to);
Qr_t = r_t.*LOy(to);

% Filter - make sure to flush out all the useful samples at the end
x_t = filter(fliplr(b),1, [ Ir_t zeros(1,fdelay*Fs) ]);
y_t = filter(fliplr(b),1, [ Qr_t zeros(1,fdelay*Fs) ]);

% Compensate for group delay
x_t = x_t(fdelay*Fs+1:end);
y_t = y_t(fdelay*Fs+1:end);

% Downsample
xk = downsample(x_t,N);
yk = downsample(y_t,N);

% And finally pick out the samples we were looking for - fiddle around for
% it is the best I can do at this point...
xk = xk(7:DataL+6);
yk = yk(7:DataL+6);


% Now do some deciding - all we need to do is check angle
figure(1);
thetak = atan2(yk,xk);
polarscatter(thetak,sqrt(yk.^2 + xk.^2));

hold on;
x0 = cell2mat(values(LUTx.forward));
y0 = cell2mat(values(LUTy.forward));
theta0 = [ atan2(y0,x0) -pi ];
polarscatter(theta0,A*ones(size(theta0)),'r+');

[ ~,idx ] = min((theta0.' - thetak).^2);
theta_hat = theta0(idx);
% polarscatter(theta_hat,1.5*A*ones(size(theta_hat)),'go');

kt = num2cell(theta0);
LUTtheta = LUT(kt,[ kx 6 ]); % map -pi to pi

s_hat = cell2mat(values(LUTtheta.forward,num2cell(theta_hat)));

% Display the message
m = m2ascii(s_hat,M);
disp(m);

% eyediagram time!
eyediagram(x_t(N:end),N,1,1);
title('x(t)');
eyediagram(y_t(N:end),N,1,1);
title('y(t)');
